/*
 *  cg.c
 *
 *  Created by Chensong Zhang on 03/28/2009.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file cg.c
 *  \brief Preconditioned Conjugate Gradient Method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "precond.h"
#include "matvec.h"

/**
 * \fn int pcg(dCSRmat *A, dvector *b, dvector *u, int MaxIt, double tol, precond *pre, int print_level)
 *	 \brief A preconditioned conjugate gradient (CG) method for solving Au=b 
 *	 \param *A	 pointer to the coefficient matrix
 *	 \param *b	 pointer to the dvector of right hand side
 *	 \param *u	 pointer to the dvector of DOFs
 *	 \param MaxIt integer, maximal number of iterations
 *	 \param tol double float, the tolerance for stopage
 *	 \param *pre pointer to the structure of precondition (precond) 
 * \param print_level how much information to print out
 *	 \return the number of iterations
 */
int pcg(dCSRmat *A, dvector *b, dvector *u, int MaxIt, double tol, precond *pre, int print_level)
{
	int iter = 0, m = A->row;
	double absres0, absres, relres1, relres2, factor;
	double alpha, beta, temp1, temp2, normb2, normu2;
	//	return 0;////////////////////////////////////////////////
	dvector p, z, r, t;
	create_dvector(m, &p);
	create_dvector(m, &z);
	create_dvector(m, &r);
	create_dvector(m, &t);
	//	double p[m], z[m], r[m], t[m];
	//return 0;/////////////////////////////
	// (b,b)
	normb2 = dot_array(m, b->val, b->val);

	// r = b-A*u
	copy_array(m, b->val, r.val);
	sparse_mv(-1.0, A, u->val, r.val);

	temp2 = dot_array(m, r.val, r.val);
	relres1 = sqrt(temp2 / normb2);
	if (relres1<tol) {
		free_dvector(&p);
		free_dvector(&z);
		free_dvector(&r);
		free_dvector(&t);
		return iter;
	}

	// z = B*r
	if (pre == NULL)
		copy_array(m, r.val, z.val); /* No preconditioner, B=I */
	else
		pre->fct(A, r.val, z.val, pre->data); /* Preconditioning */

											  // p = z
	copy_array(m, z.val, p.val);

	// temp1=(z_{k-1},r_{k-1})
	temp1 = dot_array(m, z.val, r.val);

	// initial residual
	absres0 = 1.e8;

	while (iter++ < MaxIt & normb2 > 1e-32)
	{
		// t=A*p
		init_array(m, t.val, 0.0);
		sparse_mv(1.0, A, p.val, t.val);

		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2 = dot_array(m, t.val, p.val);
		alpha = temp1 / temp2;

		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		axpy_array(m, alpha, p.val, u->val);

		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		sparse_mv(-alpha, A, p.val, r.val);

		temp2 = dot_array(m, r.val, r.val);
		if (temp2<1e-30) {
			iter = iter*(-1) - 1; break;
		}

		// relative residual = ||b-Au||_2/||b||_2=||r||_2/||b||_2
		absres = sqrt(temp2);
		normu2 = dot_dvector(u, u);
		normb2 = dot_array(m, b->val, b->val);

		relres1 = absres / sqrt(normb2);
		relres2 = absres / sqrt(normu2); // ||r||/||x|| for stopping crietia
		factor = absres / absres0;

		// output iteration information if needed	
		if (print_level>1) {
			if (iter == 1) {
				printf("It Num |  ||r||/||b|| |  ||r||/||x|| |     ||r||    |  Conv. Factor\n");
				printf("%6d | %12.5e | %12.5e | %12.5e | \n", iter, relres1, relres2, absres);
			}
			else {
				printf("%6d | %12.5e | %12.5e | %12.5e | %10.5f\n", iter, relres1, relres2, absres, factor);
			}
		}

		// update relative residual here
		absres0 = absres;

		// check convergence
		if (relres1<tol) break;

		// z_k = B*r_k
		if (pre == NULL)
			copy_array(m, r.val, z.val);	 /* No preconditioner, B=I */
		else
			pre->fct(A, r.val, z.val, pre->data); /* preconditioning */

												  // compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2 = dot_array(m, z.val, r.val);
		beta = temp2 / temp1;
		temp1 = temp2;

		// compute p_k = z_k + beta_k*p_{k-1}
		axpby_array(m, 1.0, z.val, beta, p.val);
	}

	if (iter<0) {
		iter = iter*(-1) - 1; return iter;
	}

	if (print_level>0) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, relres1);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, relres1);
	}
	free_dvector(&p);
	free_dvector(&z);
	free_dvector(&r);
	free_dvector(&t);
	return iter;
}


/**
 * \fn int den_pcg(ddenmat *A, dvector *b, dvector *u, int MaxIt, double tol, den_precond *pre, int print_level)
 *	 \brief A preconditioned conjugate gradient (CG) method for solving Au=b 
 *	 \param *A	 pointer to the coefficient matrix
 *	 \param *b	 pointer to the dvector of right hand side
 *	 \param *u	 pointer to the dvector of DOFs
 *	 \param MaxIt integer, maximal number of iterations
 *	 \param tol double float, the tolerance for stopage
 *	 \param *pre pointer to the structure of precondition (den_precond) 
 * \param print_level how much information to print out
 *	 \return the number of iterations
 */
int den_pcg(ddenmat *A, dvector *b, dvector *u, int MaxIt, double tol, den_precond *pre, int print_level)
{
	int iter=0,m=A->row;
	double alpha, beta, error, temp1, temp2, tempb;
	double *p, *z, *r, *t;
	
	p=(double *)calloc(m,sizeof(double));
	z=(double *)calloc(m,sizeof(double));
	r=(double *)calloc(m,sizeof(double));
	t=(double *)calloc(m,sizeof(double));
	
	// (b,b)
	tempb=dot_array(m,b->val,b->val);
	
	// r = b-A*u
	copy_array(m,b->val,r);
	denmat_mv(-1.0,A,u->val,r);
	
	temp2=dot_array(m,r,r);
	if(temp2<1e-30) {
		free(p);
		free(r);
		free(z);
		free(t);
		return 0;
	 }
	
	// z = B*r
	if (pre == NULL)
		copy_array(m,r,z); /* No preconditioner, B=I */
	else
		pre->fct(A,r,z,pre->data); /* Preconditioning */
	
	// p = z
	copy_array(m,z,p);
	
	// temp1=(z_{k-1},r_{k-1})
	temp1=dot_array(m,z,r);
	
	while(iter<MaxIt)
	{
		iter++;
		
		// t=A*p
		init_array(m,t,0.0);
		denmat_mv(1.0,A,p,t);
		
		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2=dot_array(m,t,p);
		alpha=temp1/temp2;
		
		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		axpy_array(m,alpha,p,u->val);
		
		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		denmat_mv(-alpha,A,p,r);
		
		temp2=dot_array(m,r,r);
		/*	if(temp2<1e-30) {
		 iter=iter*(-1)-1; break;
		 }*/
		
		// relative residual = ||b-Au||_2/||b||_2=||r||_2/||b||_2
		error=sqrt(temp2/tempb);		
		if (print_level>1)
			printf("Iteration %3d: relative residual = %e\n",iter,error);
		if (error<tol) break;
		
		// z_k = B*r_k
		if (pre == NULL)
			copy_array(m,r,z);	 /* No preconditioner, B=I */
		else
			pre->fct(A,r,z,pre->data); /* preconditioning */
		
		// compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2=dot_array(m,z,r);
		beta=temp2/temp1;
		temp1=temp2;
		
		// compute p_k = z_k + beta_k*p_{k-1}
		axpby_array(m,1.0,z,beta,p);
	}
	
	if(iter<0) {
		iter=iter*(-1)-1;
		free(p);
		free(r);
		free(z);
		free(t);
		return iter;
	}
	
	if (print_level>0) {
		if (iter>MaxIt)
			printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, error);
		else
			printf("Number of iterations = %d with relative residual %e.\n", iter, error);
	}
	
	free(p);
	free(r);
	free(z);
	free(t);
	
	return iter;
}

