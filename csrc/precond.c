/*
 *  precond.c
 *  
 *
 *  Created by Chensong on 3/28/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file precond.c
 *  \brief Preconditioners
 */
 
#include <stdio.h> 
#include <stdlib.h>

#include "precond.h"
#include "matvec.h"
 
/**
 * \fn void precond_null(dCSRmat *A, double *r, double *z, void *data)
 * \brief Do nothing preconditioner z=I*r
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 */
void precond_null(dCSRmat *A, double *r, double *z, void *data)
{
	int m=A->row;	
	copy_array(m,r,z);
}

/**
 * \fn void precond_diag(dCSRmat *A, double *r, double *z, void *data)
 * \brief Diagonal preconditioner z=inv(D)*r
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 */
void precond_diag(dCSRmat *A, double *r, double *z, void *data)
{
	int i, m=A->row;	
	dvector *diag = (dvector *)data;
	for (i=0;i<m;i++) {
		if (abs(diag->val[i])>1e-16) {
			z[i]=r[i]/diag->val[i];
		}
		else {
			z[i]=r[i];
		}
	}	
}

/**
 * \fn void precond_classicAMG(dCSRmat *A, double *r, double *z, void *data)
 * \brief get z from r by classic AMG
 * \param *A pointer to the stiffness matrix
 * \param *r pointer to residual
 * \param *z pointer to preconditioned residual
 * \param *data pointer to precondition data
 */
void precond_classicAMG(dCSRmat *A, double *r, double *z, void *data)
{
	int i, m=A->row;
	dvector rr, zz;
	precond_data *amgdata=data;
	
	int smoother = amgdata->smoother;
	int pre = amgdata->presmooth_iter;
	int post = amgdata->postsmooth_iter;
	int maxit = amgdata->max_iter;
	int max_levels = amgdata->max_levels;
	
	zz.row=rr.row=m;
	rr.val=r; zz.val=z;
	init_dvector(&zz,0.0);
	
	for(i=0;i<maxit;i++)
		multigrid(amgdata->Aarray[0], &rr, &zz, amgdata->Aarray[1], 
						  amgdata->Aarray[2], 0, max_levels, smoother, pre, post, 1);	
}

/**
* \fn void precond_aspLaplaceVec2(dCSRmat *A, double *r, double *z, void *data)
* \brief get z from r by auxiliary space preconditioner of vector Laplacian operator
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_aspLaplaceVec2(dCSRmat *A, double *r, double *z, void *data)
{
	int i, m = A->row;
	dvector rr, zz, z1, z2, r2;
	dvector z2s[2], r2s[2];
	precond_data *aspdata = data;

	int smoother = aspdata->smoother;
	int smooth_iter = aspdata->smooth_iter;
	int mg_smoother = aspdata->mg_smoother;
	int mg_smooth_iter = aspdata->mg_smooth_iter;
	int maxit = aspdata->max_iter;
	int levelNum = aspdata->max_levels;
	dCSRmat *R = aspdata->R;
	dCSRmat *P = aspdata->P;
	dCSRmat *As = aspdata->As;
	dCSRmat *Rs = aspdata->Rs;
	dCSRmat *Ps = aspdata->Ps;
	dOBDmat *swzB = aspdata->swzB;
	
	int precond_type = aspdata->precond_type; // 1 additive; 2 multiplicative  

	int m1[levelNum], m2[levelNum];
	m1[0] = smooth_iter;
	m2[0] = smooth_iter;
	for (i = 1; i < levelNum; i++)
	{
		m1[i] = mg_smooth_iter;
		m2[i] = mg_smooth_iter;
	}

	zz.row = rr.row = m;
	rr.val = r; zz.val = z;

	init_array(m, z, 0);

	create_dvector(m, &z1);
	create_dvector(R->row, &z2);
	create_dvector(R->row, &r2);
	z2s[0].row = R->row / 2;
	z2s[1].row = R->row / 2;
	r2s[0].row = R->row / 2;
	r2s[1].row = R->row / 2;
	z2s[0].val = z2.val;
	z2s[1].val = z2.val + z2s[0].row;
	r2s[0].val = r2.val;
	r2s[1].val = r2.val + r2s[0].row;

	// additive
	if (precond_type == 1)
	{
		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&z1, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&z1, 0, m - 1, 1, A, &rr, smooth_iter);
			gs(&z1, m - 1, 0, -1, A, &rr, smooth_iter);
		}
		else if (smoother == SMSWZ) {
			mulschwarz(&z1, 0, swzB->nb - 1, 1, A, &rr, swzB, smooth_iter);
			mulschwarz(&z1, swzB->nb - 1, 0, -1, A, &rr, swzB, smooth_iter);
		}

		sparse_mv0(1.0, R, r, r2.val);
		for (i = 0; i < maxit; i++)
		{
			multigrid(As, &r2s[0], &z2s[0], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		//	multigridvar(As, &r2s[0], &z2s[0], Rs, Ps, 0, levelNum, mg_smoother, m1, m2, 1);
		}
		for (i = 0; i < maxit; i++)
		{
			multigrid(As, &r2s[1], &z2s[1], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		//	multigridvar(As, &r2s[1], &z2s[1], Rs, Ps, 0, levelNum, mg_smoother, m1, m2, 1);
		}
		sparse_mv0(1.0, P, z2.val, zz.val);

		axpy_dvector(1.0, &z1, &zz);
	}
	else // multiplicative
	{
		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == GS) {
			gs(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
			gs(&zz, m - 1, 0, -1, A, &rr, smooth_iter);
		}
		else if (smoother == MSWZ) {
			mulschwarz(&zz, 0, swzB->nb - 1, 1, A, &rr, swzB, smooth_iter);
		}
		else if (smoother == SMSWZ) {
			mulschwarz(&zz, 0, swzB->nb - 1, 1, A, &rr, swzB, smooth_iter);
			mulschwarz(&zz, swzB->nb - 1, 0, -1, A, &rr, swzB, smooth_iter);
		}

		/** form residual z1 = rr - A zz */
		copy_dvector(&rr, &z1);
		sparse_mv(-1.0, A, zz.val, z1.val);

		sparse_mv0(1.0, R, z1.val, r2.val);

		for (i = 0; i < maxit; i++)
		{
			//	multigridvar(As, &r2s[0], &z2s[0], Rs, Ps, 0, levelNum, mg_smoother, m1, m2, 1);
			multigrid(As, &r2s[0], &z2s[0], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		}
		for (i = 0; i < maxit; i++)
		{
			//	multigridvar(As, &r2s[1], &z2s[1], Rs, Ps, 0, levelNum, mg_smoother, m1, m2, 1);
			multigrid(As, &r2s[1], &z2s[1], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
		}

		sparse_mv(1.0, P, z2.val, zz.val);

		/** smoothing */
		if (smoother == JACOBI) {
			jacobi(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
		}
		else if (smoother == GS) {
			gs(&zz, m - 1, 0, -1, A, &rr, smooth_iter);
		}
		else if (smoother == SGS) {
			gs(&zz, 0, m - 1, 1, A, &rr, smooth_iter);
			gs(&zz, m - 1, 0, -1, A, &rr, smooth_iter);
		}
		else if (smoother == MSWZ) {
			mulschwarz(&zz, swzB->nb - 1, 0, -1, A, &rr, swzB, smooth_iter);
		}
		else if (smoother == SMSWZ) {
			mulschwarz(&zz, 0, swzB->nb - 1, 1, A, &rr, swzB, smooth_iter);
			mulschwarz(&zz, swzB->nb - 1, 0, -1, A, &rr, swzB, smooth_iter);
		}
	}

	free_dvector(&z1);
	free_dvector(&r2);
	free_dvector(&z2);
}

/**
* \fn void precond_TriAsP1Elas(dCSRmat *A, double *r, double *z, void *data)
* \brief get z from r by block triangular preconditioner with auxiliary space method of linear elasticity
* discretized from mixed FEM for linear elasticity
*
*                 |M  B'| |sigma| = |f|
*                 |B -C | |u|     = |g|
*
* Use |D  B'| as the preconditioner in gmres and compute the inverse by
*     |B -C |
*
* the factorization
*
*                 |D  B'| |I Dinv*B'| = |D      0      |
*                 |B -C | |0   -I   |   |B B*Dinv*B'+C |
* More generally, we have
*                 |A11 A12| |I A11inv*A12| = |A11           0         |
*                 |A21 A22| |0     -I    |   |A21  A21*A11inv*A12-A22 |
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_TriAsP1Elas(dvector *r, dvector *z, void *data)
{
	int i, m = r[1].row;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *M = aspdata->precA[0];
	dCSRmat *Bt = aspdata->precA[1];
	dCSRmat *B = aspdata->precA[2];
	dCSRmat *A = aspdata->precA[3];
	double *scale = aspdata->precond_scale;
	int mass_precond_type = aspdata->mass_precond_type;

	// preconditioning r[0] by diagonal preconditioner
	if (aspdata->Minv == NULL)
	{
		for (i = 0; i < z[0].row; i++)
			z[0].val[i] = r[0].val[i] / diag->val[i];

		//	gs(&z[0], 0, z[0].row - 1, 1, M, &r[0], 3);
		//	gs(&z[0], z[0].row - 1, 0, -1, M, &r[0], 3);
	}
	else
		dBDmat_mv0(1.0, aspdata->Minv, &r[0], &z[0]);



	dvector tempVec;
	create_dvector(r[1].row, &tempVec);
	// tempVec = r[1] - B*z[0]
	copy_dvector(&r[1], &tempVec);
	sparse_mv(-1.0, B, z[0].val, tempVec.val);


	if (dot_dvector(&tempVec, &tempVec) < 1e-50)
		init_dvector(&z[1], 0);
	else
	{
		precond_aspLaplaceVec2(A, tempVec.val, z[1].val, data);
		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_aspLaplaceVec2;

		// solver part
		int iter = pcg(A, &tempVec, &z[1], 100, 1e-8, prec, 0);
		//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}

	free_dvector(&tempVec);

	int tri = 0;
	if (tri)
	{
		//		axy_dvector(-1.0, &z[1], &z[1]);/////////////
	}
	else
	{
		// z[0] = z[0] + Dinv*Bt*z[1]
		create_dvector(r[0].row, &tempVec);
		if (aspdata->Minv == NULL)
		{
			sparse_mv(1.0, Bt, z[1].val, tempVec.val);
			for (i = 0; i < tempVec.row; i++)
				tempVec.val[i] = tempVec.val[i] / diag->val[i];
		}
		else
		{
			dvector tempVec1;
			create_dvector(r[0].row, &tempVec1);
			sparse_mv(1.0, Bt, z[1].val, tempVec1.val);
			dBDmat_mv0(1.0, aspdata->Minv, &tempVec1, &tempVec);
			free_dvector(&tempVec1);
		}
		axpy_dvector(1.0, &tempVec, &z[0]);
		free_dvector(&tempVec);

		// z[1] = -z[1]
		axy_dvector(-1.0, &z[1], &z[1]);
	}


	//	printf("%lf, %lf\n", scale[0],scale[1]);
//	axy_dvector(scale[0], &z[0], &z[0]);
//	axy_dvector(scale[1], &z[1], &z[1]);

}

/**
* \fn void precond_DiagAsP1Elas(dCSRmat *A, double *r, double *z, void *data)
* \brief get z from r by block diagonal preconditioner with auxiliary space method of linear elasticity
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_DiagAsP1Elas(dvector *r, dvector *z, void *data)
{
	int i, m = r[1].row;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *M = aspdata->precA[0];
	dCSRmat *A = aspdata->precA[3];
	double *scale = aspdata->precond_scale;

	// preconditioning r[0] by diagonal preconditioner
//	for (i = 0; i < z[0].row; i++)
//		z[0].val[i] = r[0].val[i] / diag->val[i];
	if (aspdata->Minv == NULL)
	{
		for (i = 0; i < z[0].row; i++)
			z[0].val[i] = r[0].val[i] / diag->val[i];

		// preconditioning r[0] by symmetric Gaussian smoother
//		gs(&z[0], 0, z[0].row - 1, 1, M, &r[0], 3);
//		gs(&z[0], z[0].row - 1, 0, -1, M, &r[0], 3);
	}
	else
		dBDmat_mv0(1.0, aspdata->Minv, &r[0], &z[0]);
	
	//	gs(&z[0], 0, z[0].row - 1, 1, M, &r[0], 3);
	//	gs(&z[0], z[0].row - 1, 0, -1, M, &r[0], 3);

	// preconditioning r[0] by ASP
	//	copy_dvector(&r[1], &z[1]);
	//	precond_aspLaplaceVec2(Adg, r[1].val, z[1].val, data);

	//	for (i = 0; i < z[1].row; i++)
	//		z[1].val[i] *= -1;
	//	init_dvector(&z[1], 0);
	if (dot_dvector(&r[1], &r[1])<1e-50)
		init_dvector(&z[1], 0);
	else
	{
		precond_aspLaplaceVec2(A, r[1].val, z[1].val, data);

		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_aspLaplaceVec2;

		// solver part
		int iter = pcg(Adg, &r[1], &z[1], 2, 1e-8, prec, 0);
		//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}
	//	init_dvector(&z[1], 0);

	//	printf("%lf, %lf\n", scale[0],scale[1]);
//	axy_dvector(scale[0], &z[0], &z[0]);
//	axy_dvector(scale[1], &z[1], &z[1]);

	/*printf("r[1].:\n");
	for (i = 0; i < r[1].row; i++)
	printf("%lf, ", r[1].val[i]);
	printf("\n");
	printf("z[1].:\n");
	for (i = 0; i < z[1].row; i++)
	printf("%lf, ", z[1].val[i]);
	printf("\n");
	*/
}

/**
* \fn void precond_AbfpAsP1Stokes(dvector *r, dvector *z, void *data)
* \brief get z from r by approximate block factorization preconditioner with auxiliary space method of Stokes problem
* discretized from mixed FEM for linear elasticity
*
*                 |A  B'| |sigma| = |f|
*                 |B -C | |u|     = |g|
*
* Use |A-B'*Dinv*B   B'| as the preconditioner in gmres and compute the inverse by
*     |    B        -D |
*
* the factorization
*
*                 |A-B'*Dinv*B   B'| = |A  -B'| |   I      0 |
*                 |    B        -D |   |0   D | |Dinv*B   -I |
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_AbfpAsP1Stokes(dvector *r, dvector *z, void *data)
{
	int i, m = r[1].row;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *A = aspdata->precA[0];
	dCSRmat *Bt = aspdata->precA[1];
	dCSRmat *B = aspdata->precA[2];
	dCSRmat *M = aspdata->precA[3];
	double *scale = aspdata->precond_scale;

	// preconditioning r[1] by diagonal preconditioner
	if (aspdata->Minv == NULL)
	{
		for (i = 0; i < z[1].row; i++)
			z[1].val[i] = r[1].val[i] / diag->val[i];

		//	gs(&z[1], 0, z[1].row - 1, 1, M, &r[1], 3);
		//	gs(&z[1], z[1].row - 1, 0, -1, M, &r[1], 3);
	}
	else
		dBDmat_mv0(1.0, aspdata->Minv, &r[1], &z[1]);



	dvector tempVec;
	create_dvector(r[0].row, &tempVec);
	// tempVec = r[0] + Bt*z[1]
	copy_dvector(&r[0], &tempVec);
	sparse_mv(1.0, Bt, z[1].val, tempVec.val);


	if (dot_dvector(&tempVec, &tempVec) < 1e-50)
		init_dvector(&z[0], 0);
	else
	{
		precond_aspLaplaceVec2(A, tempVec.val, z[0].val, data);
		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_aspLaplaceVec2;

		// solver part
		int iter = pcg(A, &tempVec, &z[0], 100, 1e-8, prec, 0);
		//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}

	free_dvector(&tempVec);

	// z[1] = -z[1] + Dinv*B*z[0]
	axy_dvector(-1.0, &z[1], &z[1]);
	create_dvector(r[1].row, &tempVec);
	if (aspdata->Minv == NULL)
	{
		sparse_mv(1.0, B, z[0].val, tempVec.val);
		for (i = 0; i < tempVec.row; i++)
			tempVec.val[i] = tempVec.val[i] / diag->val[i];
	}
	else
	{
		dvector tempVec1;
		create_dvector(r[1].row, &tempVec1);
		sparse_mv(1.0, B, z[0].val, tempVec1.val);
		dBDmat_mv0(1.0, aspdata->Minv, &tempVec1, &tempVec);
		free_dvector(&tempVec1);
	}
	axpy_dvector(1.0, &tempVec, &z[1]);
	free_dvector(&tempVec);
}

/**
* \fn void precond_DiagAsP1Stokes(dvector *r, dvector *z, void *data)
* \brief get z from r by block diagonal preconditioner with auxiliary space method of Stokes equation
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_DiagAsP1Stokes(dvector *r, dvector *z, void *data)
{
	int i;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *A = aspdata->precA[0];
	dCSRmat *M = aspdata->precA[3];
	double *scale = aspdata->precond_scale;

	// preconditioning r[1] by diagonal preconditioner
	if (aspdata->Minv == NULL)
	{
		//	for (i = 0; i < z[1].row; i++)
	//		z[1].val[i] = r[1].val[i] / diag->val[i];

	// preconditioning r[1] by symmetric Gaussian smoother
		gs(&z[1], 0, z[1].row - 1, 1, M, &r[1], 3);
		gs(&z[1], z[1].row - 1, 0, -1, M, &r[1], 3);
	}
	else
		dBDmat_mv0(1.0, aspdata->Minv, &r[1], &z[1]);

	// preconditioning r[0] by ASP
	//	copy_dvector(&r[1], &z[1]);
	//	precond_aspLaplaceVec2(Adg, r[1].val, z[1].val, data);

	//	for (i = 0; i < z[1].row; i++)
	//		z[1].val[i] *= -1;
	//	init_dvector(&z[1], 0);
	if (dot_dvector(&r[0], &r[0])<1e-50)
		init_dvector(&z[0], 0);
	else
	{
		precond_aspLaplaceVec2(A, r[0].val, z[0].val, data);

		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_aspLaplaceVec2;

		// solver part
		int iter = pcg(Adg, &r[1], &z[1], 2, 1e-8, prec, 0);
		//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}
	//	init_dvector(&z[1], 0);

	//	printf("%lf, %lf\n", scale[0],scale[1]);
//	axy_dvector(scale[0], &z[0], &z[0]);
//	axy_dvector(scale[1], &z[1], &z[1]);

	/*printf("r[1].:\n");
	for (i = 0; i < r[1].row; i++)
	printf("%lf, ", r[1].val[i]);
	printf("\n");
	printf("z[1].:\n");
	for (i = 0; i < z[1].row; i++)
	printf("%lf, ", z[1].val[i]);
	printf("\n");
	*/
}

/**
* \fn void precond_DiagAsP1StokesNcP1_P0(dvector *r, dvector *z, void *data)
* \brief get z from r by block diagonal preconditioner with auxiliary space method of linear elasticity
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_DiagAsP1StokesNcP1_P0(dvector *r, dvector *z, void *data)
{
	int i;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *A = aspdata->precA[0];
	dCSRmat *M = aspdata->precA[3];
	double *scale = aspdata->precond_scale;

	// preconditioning r[1] by diagonal preconditioner
	for (i = 0; i < z[1].row; i++)
		z[1].val[i] = r[1].val[i] / diag->val[i];

	//	gs(&z[0], 0, z[0].row - 1, 1, M, &r[0], 3);
	//	gs(&z[0], z[0].row - 1, 0, -1, M, &r[0], 3);

	// preconditioning r[0] by ASP
	//	copy_dvector(&r[1], &z[1]);
	//	precond_aspLaplaceVec2(Adg, r[1].val, z[1].val, data);

	//	for (i = 0; i < z[1].row; i++)
	//		z[1].val[i] *= -1;
	//	init_dvector(&z[1], 0);
	if (dot_dvector(&r[0], &r[0])<1e-50)
		init_dvector(&z[0], 0);
	else
	{
		precond_aspLaplaceVec2(A, r[0].val, z[0].val, data);

		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_aspLaplaceVec2;

		// solver part
		int iter = pcg(A, &r[0], &z[0], 100, 1e-8, prec, 0);
		//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}
	//	init_dvector(&z[1], 0);

	//	printf("%lf, %lf\n", scale[0],scale[1]);
//	axy_dvector(scale[0], &z[0], &z[0]);
//	axy_dvector(scale[1], &z[1], &z[1]);

	/*printf("r[1].:\n");
	for (i = 0; i < r[1].row; i++)
	printf("%lf, ", r[1].val[i]);
	printf("\n");
	printf("z[1].:\n");
	for (i = 0; i < z[1].row; i++)
	printf("%lf, ", z[1].val[i]);
	printf("\n");
	*/
}

/**
* \fn void precond_DiagAMGStokesNcP1_P0(dvector *r, dvector *z, void *data)
* \brief get z from r by block diagonal preconditioner with auxiliary space method of linear elasticity
* \param *A pointer to the stiffness matrix
* \param *r pointer to residual
* \param *z pointer to preconditioned residual
* \param *data pointer to precondition data
*/
void precond_DiagAMGStokesNcP1_P0(dvector *r, dvector *z, void *data)
{
	int i;
	precond_data *aspdata = data;
	dvector *diag = aspdata->diag;
	dCSRmat *A = aspdata->precA[0];
	dCSRmat *M = aspdata->precA[3];
	double *scale = aspdata->precond_scale;

	// preconditioning r[1] by diagonal preconditioner
	for (i = 0; i < z[1].row; i++)
		z[1].val[i] = r[1].val[i] / diag->val[i];

	//	gs(&z[0], 0, z[0].row - 1, 1, M, &r[0], 3);
	//	gs(&z[0], z[0].row - 1, 0, -1, M, &r[0], 3);

	// preconditioning r[0] by ASP
	//	copy_dvector(&r[1], &z[1]);
	//	precond_aspLaplaceVec2(Adg, r[1].val, z[1].val, data);

	//	for (i = 0; i < z[1].row; i++)
	//		z[1].val[i] *= -1;
	//	init_dvector(&z[1], 0);
	if (dot_dvector(&r[0], &r[0])<1e-50)
		init_dvector(&z[0], 0);
	else
	{
		dvector zs[2], rs[2];
		int smoother = aspdata->smoother;
		int smooth_iter = aspdata->smooth_iter;
		int mg_smoother = aspdata->mg_smoother;
		int mg_smooth_iter = aspdata->mg_smooth_iter;
		int maxit = aspdata->max_iter;
		int levelNum = aspdata->max_levels;
		dCSRmat *As = aspdata->As;
		dCSRmat *Rs = aspdata->Rs;
		dCSRmat *Ps = aspdata->Ps;
		init_array(z[0].row, z[0].val, 0);

		zs[0].row = As[0].row;
		zs[1].row = As[0].row;
		rs[0].row = As[0].row;
		rs[1].row = As[0].row;
		zs[0].val = z[0].val;
		zs[1].val = z[0].val + zs[0].row;
		rs[0].val = r[0].val;
		rs[1].val = r[0].val + rs[0].row;

		int *var_smooth_iter;
		var_smooth_iter = (int*)calloc(levelNum - 1, sizeof(int));
		var_smooth_iter[0] = smooth_iter;
		for (i = 1; i < levelNum - 1; i++)
			var_smooth_iter[i] = mg_smooth_iter;

		for (i = 0; i < maxit; i++)
		{
			//	multigrid(As, &rs[0], &zs[0], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
			multigridvar(As, &rs[0], &zs[0], Rs, Ps, 0, levelNum, mg_smoother, var_smooth_iter, var_smooth_iter, 1);
		}
		for (i = 0; i < maxit; i++)
		{
			//	multigrid(As, &rs[1], &zs[1], Rs, Ps, 0, levelNum, mg_smoother, mg_smooth_iter, mg_smooth_iter, 1);
			multigridvar(As, &rs[1], &zs[1], Rs, Ps, 0, levelNum, mg_smoother, var_smooth_iter, var_smooth_iter, 1);
		}

		free(var_smooth_iter);

		/****  asP1ElasDG_PCG  ****
		precond *prec = (precond *)malloc(sizeof(precond));
		prec->data = data;
		prec->fct = precond_aspLaplaceVec2;

		// solver part
		int iter = pcg(Adg, &r[1], &z[1], 2, 1e-8, prec, 0);
		//		printf("iter=%d\n",iter);
		/****  asP1ElasDG_PCG  ****/
	}
	//	init_dvector(&z[1], 0);

	//	printf("%lf, %lf\n", scale[0],scale[1]);
//	axy_dvector(scale[0], &z[0], &z[0]);
//	axy_dvector(scale[1], &z[1], &z[1]);

	/*printf("r[1].:\n");
	for (i = 0; i < r[1].row; i++)
	printf("%lf, ", r[1].val[i]);
	printf("\n");
	printf("z[1].:\n");
	for (i = 0; i < z[1].row; i++)
	printf("%lf, ", z[1].val[i]);
	printf("\n");
	*/
}
