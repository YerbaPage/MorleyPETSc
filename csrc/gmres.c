/*
 *  gmres.c
 *
 *  Created by Xuehai Huang on 3/15/2016.
 *  Copyright 2016 WZU. All rights reserved.
 *
 */

/*! \file gmres.c
 *  \brief Generalized Minimum Residual method 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "precond.h"
#include "matvec.h"

/**
 * \fn void givens(int beta, dCSRmat *H, dvector *y)
 * \brief Perform Givens rotations to compute y |beta*e_1- H*y|
 * \param beta the norm of residual r_0
 * \param H (m+1)*m upper Hessenberg dCSRmat matrix 
 * \param y minimizer of |beta*e_1- H*y|
 * \return void
 *
 * Note: This is a private function used by GMRes only.
 */
void givens(double beta, dCSRmat *H, dvector *y)
{
	int Hsize=H->row;
	double *b = (double*)calloc(Hsize, sizeof(double));
	b[0]=beta;
		
	int i, j, istart, idiag, ip1start;
	double h0,h1,r,c,s,tempi,tempip1, sum;
	
	for(i=0;i<Hsize-1;i++)
	{
		istart=H->IA[i];
		ip1start=H->IA[i+1];
		if(i==0)
			idiag=istart;
		else
			idiag=istart+1;
		
		h0=H->val[idiag]; // h0=H[i][i]
		h1=H->val[H->IA[i+1]]; // h1=H[i+1][i]
		r=sqrt(h0*h0+h1*h1);
		c=h0/r;
		s=h1/r;
		
		for(j=idiag;j<ip1start;j++) {
			tempi=H->val[j];
			tempip1=H->val[ip1start+(j-idiag)];
			H->val[j]=c*tempi+s*tempip1;
			H->val[ip1start+(j-idiag)]=c*tempip1-s*tempi;
		}
		
		tempi=c*b[i]+s*b[i+1];
		tempip1=c*b[i+1]-s*b[i];
		
		b[i]=tempi;
		b[i+1]=tempip1;
	}
		
	for(i=Hsize-2;i>=0;i--)
	{
		sum=b[i];
		istart=H->IA[i];
		if(i==0)
			idiag=istart;
		else
			idiag=istart+1;
		
		for(j=Hsize-2;j>i;j--)
			sum-=H->val[idiag+j-i]*y->val[j];
		
		y->val[i]=sum/H->val[idiag];
	}
	
	free(b);
}

/**
 * \fn int gmres(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
 * \brief A preconditioned generalized minimum residual method (with restarts) (GMRES) method for solving Au=b 
 * \param *A pointer to the coefficient matrix
 * \param *b pointer to the dvector of right hand side
 * \param *u pointer to the dvector of dofs
 * \param MaxIt integer, maximal number of iterations
 * \param tol the tolerance for stopage
 * \param *pre pointer to the structure of precondition (precond) 
 * \param print_level how much information to print out
 * \return the number of iterations
 */
int gmres(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
{
	int i,j, index;
	int m=restart;
	int row=A->row;
	double beta, tempe, tempb, error, error_, beta0;
	int iter=0;
	
	if (m<1 || m>row) m=row;
	
	// compute norm^2 for right hand side
	tempb=sqrt(dot_dvector(b,b));

	dvector r;
	create_dvector(row, &r);
	
	// r = b-A*u
	copy_dvector(b, &r);
	sparse_mv(-1.0, A, u->val, r.val);
	tempe = sqrt(dot_dvector(&r, &r));

	error=tempe/tempb;
	if (error<tol) {
		free_dvector(&r); return iter;
	}
	
	dCSRmat H;
	
	// generate the structure of H, i.e. H.IA, H.JA
	H.row=m+1;
	H.col=m;
	H.nnz=m*(m+3)/2;
	H.IA=(int*)calloc(H.row+1, sizeof(int));
	H.JA=(int*)calloc(H.nnz, sizeof(int));
	H.val=(double*)calloc(H.nnz, sizeof(double));
	H.IA[1]=m;
	for(i=2;i<=H.row;i++) H.IA[i]=H.IA[i-1]+m+2-i;

	for(i=0;i<H.row;i++)
	{
		if(i==0)
			index=0;
		else
			index=i-1;
		
		for(j=H.IA[i];j<H.IA[i+1];j++) {
			H.JA[j]=index;
			index++;
		}
	}
	
	dvector v[m+1];
	for (i = 0; i<m + 1; i++)
		create_dvector(row, &v[i]);

	dvector y, z, w;
	create_dvector(m, &y);
	create_dvector(row, &z);
	create_dvector(row, &w);

	double hij;
	
	while(iter<MaxIt)
	{
		// z = B*r
		if (pre == NULL) {
			copy_dvector(&r, &z);  /* No preconditioner, B=I */
		}
		else {
			pre->fct_dvec(&r, &z, pre->data); /* Preconditioning */
		}

		beta=sqrt(dot_dvector(&z, &z));
		if (iter == 0)
			beta0 = beta;

		error = tempe / tempb;  // relative residual
		error_ = beta / beta0;  // preconditioned relative residual
		if (print_level>1)
			printf("Iteration %3d: relative residual = %e,  preconditioned relative residual = %e\n", iter, error, error_);
		if (error<tol) break;

		// v_0=z/beta
		axy_dvector(1. / beta, &z, &v[0]);

		for(j=0;j<m;j++)
		{
			// r=Av_j
			init_dvector(&r, 0.0);
			sparse_mv(1.0,A,v[j].val,r.val);
			
			// w = B*r
			if (pre == NULL) {
				copy_dvector(&r, &w);  /* No preconditioner, B=I */
			}
			else {
				pre->fct_dvec(&r, &w, pre->data); /* Preconditioning */
			}

			for(i=0;i<=j;i++)
			{
				if(i==0)
					index=0;
				else
					index=i-1;
				hij = dot_dvector(&w, &v[i]);
				H.val[H.IA[i] + j - index] = hij;
				axpy_dvector(-hij, &v[i], &w);
			}
			
			hij=sqrt(dot_dvector(&w, &w));  // h_{j+1,j}=\|w\|_2
			H.val[H.IA[j+1]]=hij;
			
			// v_{j+1}=w/h_{j+1,j}
			axy_dvector(1. / hij, &w, &v[j + 1]);
		}
				
		givens(beta, &H, &y);
		
		// u_m=u_0+ V_m*y_m
		for(i=0;i<m;i++)
			axpy_dvector(y.val[i], &v[i], u);
				
		iter++;
		
		// r = b-A*u
		copy_dvector(b,&r);
		sparse_mv(-1.0,A,u->val,r.val);
		tempe=sqrt(dot_dvector(&r,&r));
	}

	if (print_level>0) {
		if (iter >= MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e and preconditioned relative residual %e.\n", MaxIt, error, error_);
		else
			printf("Number of iterations = %d with relative residual %e and preconditioned relative residual %e.\n", iter, error, error_);
	}
	
	free_dvector(&z);
	free_dvector(&w);
	free_dvector(&r);

	for(i=0;i<m+1;i++) free_dvector(&v[i]);

	free_csr_matrix(&H);

	free_dvector(&y);

	return iter;
}

/**
* \fn int gmres2b(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
* \brief A preconditioned generalized minimum residual method (with restarts) (GMRES) method for solving Au=b in 2 blocks
* \param *A pointer to the coefficient matrix
* \param *b pointer to the dvector of right hand side
* \param *u pointer to the dvector of dofs
* \param MaxIt integer, maximal number of iterations
* \param tol the tolerance for stopage
* \param *pre pointer to the structure of precondition (precond)
* \param print_level how much information to print out
* \return the number of iterations
*/
int gmres2b(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
{
	int i, j, index, n[2];
	int m = restart;
	
	n[0] = b[0].row;
	n[1] = b[1].row;
	
	double beta, tempe, tempb, error, error_, beta0;
	int iter = 0;

	if (m<1 || m>n[0] + n[1]) m = n[0] + n[1];

	// compute norm^2 for right hand side
	tempb = sqrt(dot_dvector2b(b, b));

	dvector r[2];

	for (i = 0; i < 2; i++)
	{
		create_dvector(n[i], &r[i]);
	}

	// r = b-A*u
	copy_dvector2b(b, r);
	sparse_mv2b(-1.0, A, u, r);
	tempe = sqrt(dot_dvector2b(r, r));

	error = tempe / tempb;
	if (error<tol) {
		free_dvector(&r[0]); 
		free_dvector(&r[1]);
		return iter;
	}

	dCSRmat H;

	// generate the structure of H, i.e. H.IA, H.JA
	H.row = m + 1;
	H.col = m;
	H.nnz = m*(m + 3) / 2;
	H.IA = (int*)calloc(H.row + 1, sizeof(int));
	H.JA = (int*)calloc(H.nnz, sizeof(int));
	H.val = (double*)calloc(H.nnz, sizeof(double));
	H.IA[1] = m;
	for (i = 2; i <= H.row; i++) H.IA[i] = H.IA[i - 1] + m + 2 - i;

	for (i = 0; i<H.row; i++)
	{
		if (i == 0)
			index = 0;
		else
			index = i - 1;

		for (j = H.IA[i]; j<H.IA[i + 1]; j++) {
			H.JA[j] = index;
			index++;
		}
	}

	dvector v[m + 1][2];
	for (i = 0; i < m + 1; i++)
	{
		create_dvector(n[0], &v[i][0]);
		create_dvector(n[1], &v[i][1]);
	}

	dvector y, z[2], w[2];
	create_dvector(m, &y);
	for (i = 0; i < 2; i++)
	{
		create_dvector(n[i], &z[i]);
		create_dvector(n[i], &w[i]);
	}

	double hij;

	while (iter<MaxIt)
	{
		// z = B*r
		if (pre == NULL) {
			copy_dvector2b(r, z);  /* No preconditioner, B=I */
		}
		else {
			pre->fct_dvec(r, z, pre->data); /* Preconditioning */
		}

		beta = sqrt(dot_dvector2b(z, z));
		if (iter == 0)
			beta0 = beta;

		error = tempe / tempb;  // relative residual
		error_ = beta / beta0;  // preconditioned relative residual
		if (print_level>1)
			printf("Iteration %3d: relative residual = %e,  preconditioned relative residual = %e\n", iter, error, error_);
		if (error<tol) break;

		// v_0=z/beta
		axy_dvector2b(1. / beta, z, v[0]);

		for (j = 0; j<m; j++)
		{
			// r=Av_j
			init_dvector2b(r, 0.0);
			sparse_mv2b(1.0, A, v[j], r);

			// w = B*r
			if (pre == NULL) {
				copy_dvector2b(r, w);  /* No preconditioner, B=I */
			}
			else {
				pre->fct_dvec(r, w, pre->data); /* Preconditioning */
			}

			for (i = 0; i <= j; i++)
			{
				if (i == 0)
					index = 0;
				else
					index = i - 1;
				hij = dot_dvector2b(w, v[i]);
				H.val[H.IA[i] + j - index] = hij;
				axpy_dvector2b(-hij, v[i], w);
			}

			hij = sqrt(dot_dvector2b(w, w));  // h_{j+1,j}=\|w\|_2
			H.val[H.IA[j + 1]] = hij;

			// v_{j+1}=w/h_{j+1,j}
			axy_dvector2b(1. / hij, w, v[j + 1]);
		}

		givens(beta, &H, &y);

		// u_m=u_0+ V_m*y_m
		for (i = 0; i<m; i++)
			axpy_dvector2b(y.val[i], v[i], u);

		iter++;

		// r = b-A*u
		copy_dvector2b(b, r);
		sparse_mv2b(-1.0, A, u, r);
		tempe = sqrt(dot_dvector2b(r, r));
	}

	if (print_level>0) {
		if (iter >= MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e and preconditioned relative residual %e.\n", MaxIt, error, error_);
		else
			printf("Number of iterations = %d with relative residual %e and preconditioned relative residual %e.\n", iter, error, error_);
	}

	for (i = 0; i < 2; i++)
	{
		free_dvector(&z[i]);
		free_dvector(&w[i]);
		free_dvector(&r[i]);
	}

	for (i = 0; i < m + 1; i++)
	{
		free_dvector(&v[i][0]);
		free_dvector(&v[i][1]);
	}

	free_csr_matrix(&H);

	free_dvector(&y);

	return iter;
}

/**
* \fn int fgmres(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
* \brief A preconditioned flexible generalized minimum residual method (with restarts) (GMRES) method for solving Au=b
* \param *A pointer to the coefficient matrix
* \param *b pointer to the dvector of right hand side
* \param *u pointer to the dvector of dofs
* \param MaxIt integer, maximal number of iterations
* \param tol the tolerance for stopage
* \param *pre pointer to the structure of precondition (precond) (right preconditioning)
* \param print_level how much information to print out
* \return the number of iterations
*/
int fgmres(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
{
	int i, j, k, n;
	int m = restart;

	n = b[0].row;

	double beta, absres0, absres, absresOld, relres, delta;
	int converge;
	int iter = 0;

	// restart number
	if (m<1 || m>n) m = n;
	if (m > MaxIt) m = MaxIt;

	dvector r;
	dvector v[m + 1], z[m];
	dvector Hj, g, y, c, s;
	ddenmat R;

	// initial memory
	create_dvector(n, &r);
	for (i = 0; i < m + 1; i++)
	{
		create_dvector(n, &v[i]);
	}
	for (i = 0; i < m; i++)
	{
		create_dvector(n, &z[i]);
	}
	create_dvector(m + 1, &Hj);
	create_dvector(m + 1, &g);
	create_dvector(m, &y);
	create_dvector(m, &c);
	create_dvector(m, &s);
	create_dden_matrix(m, m, &R);


	// r = b-A*u
	copy_dvector(b, &r);
	sparse_mv(-1.0, A, u->val, r.val);
	beta = sqrt(dot_dvector(&r, &r));

	absres0 = beta;
	absres = beta;

	if (print_level > 0)
	{
		printf("It Num | ||r||/||r0|| |    ||r||     | Conv. Factor \n");
		printf("%6d | %12.5e | %12.5e | %f\n", 0, 1.0, absres, 0.0);
	}

	// Main Loop
	while (iter<MaxIt)
	{
		// reset converge
		converge = 0;

		// v_0=r/beta
		axy_dvector(1. / beta, &r, &v[0]);

		// form right hand side for the hessenberg system
		g.val[0] = beta;

		// loop for restart
		for (j = 0; j<m; j++)
		{
			// z_j = B*v_j
			if (pre == NULL) {
				copy_dvector(&v[j], &z[j]);  /* No preconditioner, B=I */
			}
			else {
				pre->fct_dvec(&v[j], &z[j], pre->data); /* Preconditioning */
			}

			// w=Az_j
			init_dvector(&v[j + 1], 0.0);
			sparse_mv(1.0, A, z[j].val, v[j + 1].val);

			// modified Gram-Schmidt
			for (i = 0; i <= j; i++)
			{
				Hj.val[i] = dot_dvector(&v[j + 1], &v[i]);
				axpy_dvector(-Hj.val[i], &v[i], &v[j + 1]);
			}

			// new orthonormal basis
			Hj.val[j + 1] = sqrt(dot_dvector(&v[j + 1], &v[j + 1]));
			axy_dvector(1. / Hj.val[j + 1], &v[j + 1], &v[j + 1]); // becareful small Hj.val[j+1]

																   // Use Givens transformation to get upper triangular system R
			R.val[0][j] = Hj.val[0];

			// apply the previous Givens transformations
			if (j > 0)
			{
				for (i = 1; i <= j; i++)
				{
					R.val[i][j] = -s.val[i - 1] * R.val[i - 1][j] + c.val[i - 1] * Hj.val[i];
					R.val[i - 1][j] = c.val[i - 1] * R.val[i - 1][j] + s.val[i - 1] * Hj.val[i];
					/*					temp = c.val[i - 1] * R.val[i - 1][j] + s.val[i - 1] * Hj.val[i];
					R.val[i][j] = -s.val[i - 1] * R.val[i - 1][j] + c.val[i - 1] * Hj.val[i];
					R.val[i - 1][j] = temp;*/
				}
			}

			// new Givens transformation
			delta = sqrt(R.val[j][j] * R.val[j][j] + Hj.val[j + 1] * Hj.val[j + 1]);
			c.val[j] = R.val[j][j] / delta;
			s.val[j] = Hj.val[j + 1] / delta;

			R.val[j][j] = c.val[j] * R.val[j][j] + s.val[j] * Hj.val[j + 1];

			// apply Givens transformation to Right hand side g
			g.val[j + 1] = -s.val[j] * g.val[j];
			g.val[j] = c.val[j] * g.val[j];

			// count iterations
			iter++;

			// check convergence g[j+1]=||b-Au_i|
			absresOld = absres;
			absres = fabs(g.val[j + 1]);
			relres = absres / absres0;
			if (print_level>0)
				printf("%6d | %12.5e | %12.5e | %f\n", iter, relres, absres, absres / absresOld);
			if (relres < tol)
			{
				converge = 1;
				break;
			}
		} // j

		if (j == m)
			j--;

		// solve the upper trangular matrix
		for (i = j; i >= 0; i--)
		{
			y.val[i] = g.val[i];
			for (k = j; k>i; k--)
				y.val[i] -= R.val[i][k] * y.val[k];
			y.val[i] /= R.val[i][i];
		}

		// solution: u_j=u_0+ Z_j*y_j
		for (i = 0; i <= j; i++)
			axpy_dvector(y.val[i], &z[i], u);

		// update residual: r = b-A*u and restart
		copy_dvector(b, &r);
		sparse_mv(-1.0, A, u->val, r.val);
		beta = sqrt(dot_dvector(&r, &r));
		// check convergence
		if (converge && beta / absres0<tol)
			break;
	} // main loop

	if (print_level>0) {
		if (iter >= MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e and residual %e.\n", MaxIt, relres, absres);
		else
			printf("Number of iterations = %d with relative residual %e and residual %e.\n", iter, relres, absres);
	}

	free_dvector(&r);

	for (i = 0; i < m + 1; i++)
		free_dvector(&v[i]);

	for (i = 0; i < m; i++)
		free_dvector(&z[i]);

	free_dvector(&Hj);
	free_dvector(&g);
	free_dvector(&y);
	free_dvector(&c);
	free_dvector(&s);
	free_dden_matrix(&R);

	return iter;
}

/**
* \fn int fgmres_den(ddenmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
* \brief A preconditioned flexible generalized minimum residual method (with restarts) (GMRES) method for solving Au=b
* \param *A pointer to the coefficient matrix
* \param *b pointer to the dvector of right hand side
* \param *u pointer to the dvector of dofs
* \param MaxIt integer, maximal number of iterations
* \param tol the tolerance for stopage
* \param *pre pointer to the structure of precondition (precond) (right preconditioning)
* \param print_level how much information to print out
* \return the number of iterations
*/
int fgmres_den(ddenmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
{
	int i, j, k, n;
	int m = restart;

	n = b[0].row;

	double beta, absres0, absres, absresOld, relres, delta;
	int converge;
	int iter = 0;

	// restart number
	if (m<1 || m>n) m = n;
	if (m > MaxIt) m = MaxIt;

	dvector r;
	dvector v[m + 1], z[m];
	dvector Hj, g, y, c, s;
	ddenmat R;

	// initial memory
	create_dvector(n, &r);
	for (i = 0; i < m + 1; i++)
	{
		create_dvector(n, &v[i]);
	}
	for (i = 0; i < m; i++)
	{
		create_dvector(n, &z[i]);
	}
	create_dvector(m + 1, &Hj);
	create_dvector(m + 1, &g);
	create_dvector(m, &y);
	create_dvector(m, &c);
	create_dvector(m, &s);
	create_dden_matrix(m, m, &R);


	// r = b-A*u
	copy_dvector(b, &r);
	denmat_mv(-1.0, A, u->val, r.val);
	beta = sqrt(dot_dvector(&r, &r));

	absres0 = beta;
	absres = beta;

	if (print_level > 0)
	{
		printf("It Num | ||r||/||r0|| |    ||r||     | Conv. Factor \n");
		printf("%6d | %12.5e | %12.5e | %f\n", 0, 1.0, absres, 0.0);
	}

	// Main Loop
	while (iter<MaxIt)
	{
		// reset converge
		converge = 0;

		// v_0=r/beta
		axy_dvector(1. / beta, &r, &v[0]);

		// form right hand side for the hessenberg system
		g.val[0] = beta;

		// loop for restart
		for (j = 0; j<m; j++)
		{
			// z_j = B*v_j
			if (pre == NULL) {
				copy_dvector(&v[j], &z[j]);  /* No preconditioner, B=I */
			}
			else {
				pre->fct_dvec(&v[j], &z[j], pre->data); /* Preconditioning */
			}

			// w=Az_j
			init_dvector(&v[j + 1], 0.0);
			denmat_mv(1.0, A, z[j].val, v[j + 1].val);

			// modified Gram-Schmidt
			for (i = 0; i <= j; i++)
			{
				Hj.val[i] = dot_dvector(&v[j + 1], &v[i]);
				axpy_dvector(-Hj.val[i], &v[i], &v[j + 1]);
			}

			// new orthonormal basis
			Hj.val[j + 1] = sqrt(dot_dvector(&v[j + 1], &v[j + 1]));
			axy_dvector(1. / Hj.val[j + 1], &v[j + 1], &v[j + 1]); // becareful small Hj.val[j+1]

																   // Use Givens transformation to get upper triangular system R
			R.val[0][j] = Hj.val[0];

			// apply the previous Givens transformations
			if (j > 0)
			{
				for (i = 1; i <= j; i++)
				{
					R.val[i][j] = -s.val[i - 1] * R.val[i - 1][j] + c.val[i - 1] * Hj.val[i];
					R.val[i - 1][j] = c.val[i - 1] * R.val[i - 1][j] + s.val[i - 1] * Hj.val[i];
					/*					temp = c.val[i - 1] * R.val[i - 1][j] + s.val[i - 1] * Hj.val[i];
					R.val[i][j] = -s.val[i - 1] * R.val[i - 1][j] + c.val[i - 1] * Hj.val[i];
					R.val[i - 1][j] = temp;*/
				}
			}

			// new Givens transformation
			delta = sqrt(R.val[j][j] * R.val[j][j] + Hj.val[j + 1] * Hj.val[j + 1]);
			c.val[j] = R.val[j][j] / delta;
			s.val[j] = Hj.val[j + 1] / delta;

			R.val[j][j] = c.val[j] * R.val[j][j] + s.val[j] * Hj.val[j + 1];

			// apply Givens transformation to Right hand side g
			g.val[j + 1] = -s.val[j] * g.val[j];
			g.val[j] = c.val[j] * g.val[j];

			// count iterations
			iter++;

			// check convergence g[j+1]=||b-Au_i|
			absresOld = absres;
			absres = fabs(g.val[j + 1]);
			relres = absres / absres0;
			if (print_level>0)
				printf("%6d | %12.5e | %12.5e | %f\n", iter, relres, absres, absres / absresOld);
			if (relres < tol)
			{
				converge = 1;
				break;
			}
		} // j

		if (j == m)
			j--;

		// solve the upper trangular matrix
		for (i = j; i >= 0; i--)
		{
			y.val[i] = g.val[i];
			for (k = j; k>i; k--)
				y.val[i] -= R.val[i][k] * y.val[k];
			y.val[i] /= R.val[i][i];
		}

		// solution: u_j=u_0+ Z_j*y_j
		for (i = 0; i <= j; i++)
			axpy_dvector(y.val[i], &z[i], u);

		// update residual: r = b-A*u and restart
		copy_dvector(b, &r);
		denmat_mv(-1.0, A, u->val, r.val);
		beta = sqrt(dot_dvector(&r, &r));
		// check convergence
		if (converge && beta / absres0<tol)
			break;
	} // main loop

	if (print_level>0) {
		if (iter >= MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e and residual %e.\n", MaxIt, relres, absres);
		else
			printf("Number of iterations = %d with relative residual %e and residual %e.\n", iter, relres, absres);
	}

	free_dvector(&r);

	for (i = 0; i < m + 1; i++)
		free_dvector(&v[i]);

	for (i = 0; i < m; i++)
		free_dvector(&z[i]);

	free_dvector(&Hj);
	free_dvector(&g);
	free_dvector(&y);
	free_dvector(&c);
	free_dvector(&s);
	free_dden_matrix(&R);

	return iter;
}

/**
* \fn int fgmres2b(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
* \brief A preconditioned flexible generalized minimum residual method (with restarts) (GMRES) method for solving Au=b in 2 blocks
* \param *A pointer to the coefficient matrix
* \param *b pointer to the dvector of right hand side
* \param *u pointer to the dvector of dofs
* \param MaxIt integer, maximal number of iterations
* \param tol the tolerance for stopage
* \param *pre pointer to the structure of precondition (precond) (right preconditioning)
* \param print_level how much information to print out
* \return the number of iterations
*/
int fgmres2b(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level)
{
	int i, j, k, n[2];
	int m = restart;

	n[0] = b[0].row;
	n[1] = b[1].row;

	double beta, absres0, absres, absresOld, relres, delta;
	int converge;
	int iter = 0;

	// restart number
	if (m<1 || m>n[0] + n[1]) m = n[0] + n[1];
	if (m > MaxIt) m = MaxIt;

	dvector r[2];
	dvector v[m + 1][2], z[m][2];
	dvector Hj, g, y, c, s;
	ddenmat R;

	// initial memory
	for (i = 0; i < 2; i++)
	{
		create_dvector(n[i], &r[i]);
	}
	for (i = 0; i < m + 1; i++)
	{
		create_dvector(n[0], &v[i][0]);
		create_dvector(n[1], &v[i][1]);
	}
	for (i = 0; i < m; i++)
	{
		create_dvector(n[0], &z[i][0]);
		create_dvector(n[1], &z[i][1]);
	}
	create_dvector(m + 1, &Hj);
	create_dvector(m + 1, &g);
	create_dvector(m, &y);
	create_dvector(m, &c);
	create_dvector(m, &s);
	create_dden_matrix(m, m, &R);


	// r = b-A*u
	copy_dvector2b(b, r);
	sparse_mv2b(-1.0, A, u, r);
	beta = sqrt(dot_dvector2b(r, r));

	absres0 = beta;
	absres = beta;

	if (print_level > 0)
	{
		printf("It Num | ||r||/||r0|| |    ||r||     | Conv. Factor \n");
		printf("%6d | %12.5e | %12.5e | %f\n", 0, 1.0, absres, 0.0);
	}

	// Main Loop
	while (iter<MaxIt)
	{
		// reset converge
		converge = 0;

		// v_0=r/beta
		axy_dvector2b(1. / beta, r, v[0]);

		// form right hand side for the hessenberg system
		g.val[0] = beta;

		// loop for restart
		for (j = 0; j<m; j++)
		{
			// z_j = B*v_j
			if (pre == NULL) {
				copy_dvector2b(v[j], z[j]);  /* No preconditioner, B=I */
			}
			else {
				pre->fct_dvec(v[j], z[j], pre->data); /* Preconditioning */
			}

			// w=Az_j
			init_dvector2b(v[j + 1], 0.0);
			sparse_mv2b(1.0, A, z[j], v[j+1]);

			// modified Gram-Schmidt
			for (i = 0; i <= j; i++)
			{
				Hj.val[i] = dot_dvector2b(v[j + 1], v[i]);
				axpy_dvector2b(-Hj.val[i], v[i], v[j + 1]);
			}

			// new orthonormal basis
			Hj.val[j+1] = sqrt(dot_dvector2b(v[j + 1], v[j + 1]));
			axy_dvector2b(1. / Hj.val[j + 1], v[j + 1], v[j + 1]); // becareful small Hj.val[j+1]

			// Use Givens transformation to get upper triangular system R
			R.val[0][j] = Hj.val[0];
			
			// apply the previous Givens transformations
			if (j > 0)
			{
				for (i = 1; i <= j; i++)
				{
					R.val[i][j] = -s.val[i - 1] * R.val[i - 1][j] + c.val[i - 1] * Hj.val[i];
					R.val[i - 1][j] = c.val[i - 1] * R.val[i - 1][j] + s.val[i - 1] * Hj.val[i];
/*					temp = c.val[i - 1] * R.val[i - 1][j] + s.val[i - 1] * Hj.val[i];
					R.val[i][j] = -s.val[i - 1] * R.val[i - 1][j] + c.val[i - 1] * Hj.val[i];
					R.val[i - 1][j] = temp;*/
				}
			}

			// new Givens transformation
			delta = sqrt(R.val[j][j] * R.val[j][j] + Hj.val[j + 1] * Hj.val[j + 1]);
			c.val[j] = R.val[j][j] / delta;
			s.val[j] = Hj.val[j + 1] / delta;

			R.val[j][j] = c.val[j] * R.val[j][j] + s.val[j] * Hj.val[j+1];

			// apply Givens transformation to Right hand side g
			g.val[j + 1] = -s.val[j] * g.val[j];
			g.val[j] = c.val[j] * g.val[j];

			// count iterations
			iter++;

			// check convergence g[j+1]=||b-Au_i|
			absresOld = absres;
			absres = fabs(g.val[j + 1]);
			relres = absres / absres0;
			if (print_level>0)
				printf("%6d | %12.5e | %12.5e | %f\n", iter, relres, absres, absres / absresOld);
			if (relres < tol)
			{
				converge = 1;
				break;
			}
		} // j

		if (j == m)
			j--;

		// solve the upper trangular matrix
		for (i = j; i >= 0; i--)
		{
			y.val[i] = g.val[i];
			for (k = j; k>i; k--)
				y.val[i] -= R.val[i][k] * y.val[k];
			y.val[i] /= R.val[i][i];
		}

		// solution: u_j=u_0+ Z_j*y_j
		for (i = 0; i<=j; i++)
			axpy_dvector2b(y.val[i], z[i], u);

		// update residual: r = b-A*u and restart
		copy_dvector2b(b, r);
		sparse_mv2b(-1.0, A, u, r);
		beta = sqrt(dot_dvector2b(r, r));
		// check convergence
		if (converge && beta / absres0<tol)
			break;
	} // main loop

	if (print_level>0) {
		if (iter >= MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e and residual %e.\n", MaxIt, relres, absres);
		else
			printf("Number of iterations = %d with relative residual %e and residual %e.\n", iter, relres, absres);
	}

	for (i = 0; i < 2; i++)
		free_dvector(&r[i]);

	for (i = 0; i < m + 1; i++)
	{
		free_dvector(&v[i][0]);
		free_dvector(&v[i][1]);
	}
	for (i = 0; i < m; i++)
	{
		free_dvector(&z[i][0]);
		free_dvector(&z[i][1]);
	}

	free_dvector(&Hj);
	free_dvector(&g);
	free_dvector(&y);
	free_dvector(&c);
	free_dvector(&s);
	free_dden_matrix(&R);

	return iter;
}

/**
 *	\fn int gmresBCD(dCSRmat *A, dvector *b, dvector *u, int MaxIt, double tol, precond *pre, int print_level)
 *	\brief A preconditioned generalized minimum residual method (with restarts) (GMRES) method for solving (D-CB)u=b 
 *	\param *A	 pointer to the coefficient matrix
 *	\param *b	 pointer to the dvector of right hand side
 *	\param *u	 pointer to the dvector of dofs
 *	\param MaxIt integer, maximal number of iterations
 *	\param tol double float, the tolerance for stopage
 *	\param *pre pointer to the structure of precondition (precond) 
 *\param print_level how much information to print out
 *	\return the number of iterations
 */
int gmresBCD(dCSRmat *B, dCSRmat *C, dCSRmat *D, dvector *b, dvector *u, int restart, int MaxIt, double tol, int print_level)
{
	int i,j, index;
	int m=restart;
	int row=D->row;
	double *r=(double*)calloc(row, sizeof(double));
	double *ftemp=(double*)calloc(B->row, sizeof(double));
	double beta, tempe, tempb, error;
	int iter=0;
	
	if (m<1 || m>row) m=row;
	
	// compute norm^2 for right hand side
	tempb=sqrt(dot_array(row,b->val,b->val));
	
	// r = b-A*u
	copy_array(row,b->val,r);
	sparse_mv(-1.0,D,u->val,r);
	sparse_mv0(1.0, B, u->val, ftemp);
	sparse_mv(1.0,C,ftemp,r);
	tempe=sqrt(dot_array(row,r,r));
	
	error=tempe/tempb;
	if (error<tol) {
		free(r); free(ftemp); return iter;
	}
	
	dCSRmat H;
	
	// generate the structure of H, i.e. H.IA, H.JA
	H.row=m+1;
	H.col=m;
	H.nnz=m*(m+3)/2;
	H.IA=(int*)calloc(H.row+1, sizeof(int));
	H.JA=(int*)calloc(H.nnz, sizeof(int));
	H.val=(double*)calloc(H.nnz, sizeof(double));
	H.IA[1]=m;
	for(i=2;i<=H.row;i++) H.IA[i]=H.IA[i-1]+m+2-i;

	for(i=0;i<H.row;i++)
	{
		if(i==0)
			index=0;
		else
			index=i-1;
		
		for(j=H.IA[i];j<H.IA[i+1];j++) {
			H.JA[j]=index;
			index++;
		}
	}
	
	dvector v[m+1];
	for(i=0;i<m+1;i++) {
		v[i].row=row;
		v[i].val=(double*)calloc(row, sizeof(double));
	}

	dvector y;
	y.row=m;
	y.val=(double*)calloc(y.row, sizeof(double));
	
	double *z, *w, hij;
	z=(double *)calloc(row,sizeof(double));
	w=(double *)calloc(row,sizeof(double));
	
	while(iter<MaxIt)
	{
		// z = B*r
		copy_array(row,r,z);  /* No preconditioner, B=I */
		
		beta=sqrt(dot_array(row,z,z));
		
		// v_0=z/beta
		for(i=0;i<row;i++) v[0].val[i]=z[i]/beta;
				
		for(j=0;j<m;j++)
		{
			// r=Av_j
			init_array(row, r, 0.0);
			sparse_mv(1.0,D,v[j].val,r);
			sparse_mv0(1.0, B, v[j].val, ftemp);
			sparse_mv(-1.0,C,ftemp,r);
	
			
			// w = B*r
			copy_array(row,r,w);  /* No preconditioner, B=I */
			
			for(i=0;i<=j;i++)
			{
				if(i==0)
					index=0;
				else
					index=i-1;
				hij=dot_array(row,w,v[i].val);
				H.val[H.IA[i]+j-index]=hij;
				axpy_array(row, -hij, v[i].val, w);
			}
			
			hij=sqrt(dot_array(row,w,w));  // h_{j+1,j}=\|w\|_2
			H.val[H.IA[j+1]]=hij;
			
			// v_{j+1}=w/h_{j+1,j}
			for(i=0;i<row;i++) v[j+1].val[i]=w[i]/hij;
		}
				
		givens(beta, &H, &y);
		
		// u_m=u_0+ V_m*y_m
		for(i=0;i<m;i++)
			axpy_array(row, y.val[i], v[i].val, u->val);
				
		iter++;
		
		// r = b-A*u
		copy_array(row,b->val,r);
		sparse_mv(-1.0,D,u->val,r);
		sparse_mv0(1.0, B, u->val, ftemp);
		sparse_mv(1.0,C,ftemp,r);
		tempe=sqrt(dot_array(row,r,r));
		
		error=tempe/tempb;
		if (print_level>1)
			printf("Iteration %3d: relative residual = %e\n",iter,error);
		if(error<tol) break;
	}

	if (print_level>0) {
		if (iter>=MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, error);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, error);
	}
	
	free(z); 
	free(w);
	free(r);
	free(ftemp);

	for(i=0;i<m+1;i++) free(v[i].val);

	free(H.IA);
	free(H.JA);
	free(H.val);

	free(y.val);

	return iter;
}
