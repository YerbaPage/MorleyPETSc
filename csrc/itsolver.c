/*
 *  itsolver.c
 *  
 *
 *  Created by Chensong on 04/01/2009.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file itsolver.c
 *  \brief Iterative solvers (main file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "precond.h"
#include "matvec.h"

static int itsolverInCG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int itsolver_type, int itsolver_maxit, double itsolver_tol, int print_level);

/**
 * \fn int standard_CG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level)
 * \brief Solve Ax=b by standard conjugate gradient method (CG) 
 * \param *A	pointer to the dCSRmat matrix
 * \param *b	pointer to the dvector of right hand side
 * \param *x	pointer to the dvector of dofs
 * \param MaxIt  integer, maximal number of iterations
 * \param tol	double float, the tolerance for stopage of the iteration
 * \param print_level how much information to print out
 * \return the number of iterations
 */
int standard_CG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level)
{
	int iter;
	clock_t solver_start, solver_end;
	
	// solver part
	solver_start=clock();
	iter=pcg(A, b, x, MaxIt, tol, NULL, print_level);		
	solver_end=clock();
	
	double solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("CG costs %f seconds.\n", solver_duration);
		printf("Number of iterations = %d.\n", iter);
	}
	
	return iter;
}

/**
 * \fn int diag_PCG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level)
 * \brief Solve Ax=b by diagonal preconditioned conjugate gradient method (PCG) 
 * \param *A	pointer to the dCSRmat matrix
 * \param *b	pointer to the dvector of right hand side
 * \param *x	pointer to the dvector of dofs
 * \param MaxIt  integer, maximal number of iterations
 * \param tol	double float, the tolerance for stopage of the iteration
 * \param print_level how much information to print out 
 * \return the number of iterations
 */
int diag_PCG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level)
{
	int iter;
	clock_t solver_start, solver_end;
	
	dvector diag;	
	getdiag(0,A,&diag);
	
	// setup preconditioner
	precond *prec = (precond *) malloc(sizeof(precond));	
	prec->data = &diag;
	prec->fct  = precond_diag;

	// solver part
	solver_start=clock();
	iter=pcg(A, b, x, MaxIt, tol, prec, print_level);		
	solver_end=clock();
	
	double solver_duration = (double)(solver_end - solver_start)/(double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Diagonal preconditioned CG costs %f seconds.\n", solver_duration);
		printf("Number of iterations = %d.\n", iter);
	}
	
	free_dvector(&diag);
	return iter;
}

/**
 * \fn int classicAMG_PCG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level)
 * \brief Solve Ax=b by preconditioned conjugate gradient method (PCG), 
 * with Ruge and Stuben's classic AMG as precondition
 * \param *A	pointer to the dCSRmat matrix
 * \param *b	pointer to the dvector of right hand side
 * \param *x	pointer to the dvector of dofs
 * \param *param pointer to AMG parameters
 * \param print_level how much information to print out
 * \return the number of iterations
 */
int classicAMG_PCG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level)
{
	int max_levels=param->max_levels;
	int levelNum=1, iter, i, l;
	
	dCSRmat AA[max_levels];
	dCSRmat P[max_levels], PT[max_levels];	
	
	clock_t classicAMG_start, classicAMG_end;
	
	// initialize
	for(i=0;i<max_levels;i++)
	{
		AA[i].row=0;
		AA[i].col=0;
		AA[i].IA=NULL;
		AA[i].JA=NULL;
		AA[i].val=NULL;
		
		P[i].row=0;
		P[i].col=0;
		P[i].IA=NULL;
		P[i].JA=NULL;
		P[i].val=NULL;
		
		PT[i].row=0;
		PT[i].col=0;
		PT[i].IA=NULL;
		PT[i].JA=NULL;
		PT[i].val=NULL;
	}
	
	AA[0].row=A->row;
	AA[0].col=A->col;
	AA[0].IA=(int*)calloc(AA[0].row+1, sizeof(int));
	
	for(i=0;i<=AA[0].row;i++) AA[0].IA[i]=A->IA[i];
	AA[0].JA=(int*)calloc(AA[0].IA[AA[0].row], sizeof(int));
	AA[0].val=(double*)calloc(AA[0].IA[AA[0].row], sizeof(double));

	for(i=0;i<AA[0].IA[AA[0].row];i++) {
		AA[0].JA[i]=A->JA[i];
		AA[0].val[i]=A->val[i];
	}
	
	// setup preconditioner
	classicAMG_setup(AA, P, PT, &levelNum, param);

	precond_data amgData;
	amgData.max_levels = levelNum;
	amgData.max_iter = param->AMG_max_iter;
	amgData.tol      = param->AMG_tol;
	amgData.smoother = param->smoother;
	amgData.presmooth_iter  = param->presmooth_iter;
	amgData.postsmooth_iter = param->postsmooth_iter;
	amgData.postsmooth_iter = param->coarsening_type;
	amgData.Aarray    = (dCSRmat **)malloc(3*sizeof(dCSRmat*));
	amgData.Aarray[0] = AA;
	amgData.Aarray[1] = PT;
	amgData.Aarray[2] = P;
	
	precond *prec = (precond *)malloc(sizeof(precond));	
	prec->data = &amgData;
	prec->fct = precond_classicAMG;
	
	// solver part
	classicAMG_start=clock();
	
	int MaxIt=param->max_iter;
	double tol=param->tol;
	iter=pcg(A, b, x, MaxIt, tol, prec, print_level);		
	
	classicAMG_end=clock();
	
	double classicAMGduration = (double)(classicAMG_end - classicAMG_start)/(double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Classic AMG preconditioned CG costs %f seconds.\n", classicAMGduration);
		printf("Number of iterations = %d.\n", iter);		
	}
	
	for(l=0;l<levelNum;l++)
	{
		free(AA[l].IA);
		free(AA[l].JA);
		free(AA[l].val);
	}
	
	for(l=0;l<levelNum-1;l++)
	{
		free(P[l].IA);
		free(P[l].JA);
		free(P[l].val);
		free(PT[l].IA);
		free(PT[l].JA);
		free(PT[l].val);
	}
	
	free(amgData.Aarray);
	free(prec);
	
	return iter;
}

/**
* \fn int DiagAsP1StokesCR_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block diagonal preconditioned MINRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int DiagAsP1StokesCR_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[0];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assemblePressMassmatrixStokesCR2d(&M, elements, &elementDOF[1], 1);
	dvector diag;
	getdiag(M.row, &M, &diag);

	dBDmat dbdM, Minv;
	assemblePressMassmatrixStokesCRdBD2d(&dbdM, elements, &elementDOF[1], 1);
	inverse_dBDmat(&dbdM, &Minv);
	free_dbd_matrix(&dbdM);

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toCR_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &A[0], edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &A[0], elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &A[0], elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &A[0];
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &M; // mass matrix
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = &Minv;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_DiagAsP1Stokes;

	// solver part
	solve_start = clock();

	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = minres2b(A, b, x, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Auxiliary space preconditioned MINRES costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dbd_matrix(&Minv);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int DiagAsP1StokesMINI_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block diagonal preconditioned MINRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int DiagAsP1StokesMINI_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFmini = &elementDOF[0];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assemblePressMassmatrixStokesMINI2d(&M, elements, &elementDOF[1], &elementdofTran[1], 1);
	dvector diag;
	getdiag(M.row, &M, &diag);

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toMINI_2d(&tempA, &elementDOFas, elementDOFmini);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFmini, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFmini);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &A[0], edges, elementDOFmini);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &A[0], elementEdge, edges, elementDOFmini);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &A[0], elements, edges, nodes->row, elementDOFmini);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &A[0], elements, nodes->row, elementDOFmini);
		else
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFmini);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &A[0];
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &M; // mass matrix
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_DiagAsP1Stokes;

	// solver part
	solve_start = clock();

	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = minres2b(A, b, x, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Auxiliary space preconditioned MINRES costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	aspData.Minv = NULL;
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int DiagAsP1StokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block diagonal preconditioned MINRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int DiagAsP1StokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[0];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assemblePressMassmatrixStokesNcP1P02d(&M, elements, &elementDOF[1]);
	dvector diag;
	create_dvector(M.row, &diag);
	for (i = 0; i < diag.row; i++)
		diag.val[i] = M.val[i];


	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toNcP1_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ ||  param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &A[0], edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &A[0], elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &A[0], elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &A[0];
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &M; // mass matrix
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_DiagAsP1StokesNcP1_P0;

	// solver part
	solve_start = clock();

	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = minres2b(A, b, x, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Auxiliary space preconditioned MINRES costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int DiagAMGStokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block diagonal preconditioned MINRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int DiagAMGStokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[0];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assemblePressMassmatrixStokesNcP1P02d(&M, elements, &elementDOF[1]);
	dvector diag;
	create_dvector(M.row, &diag);
	for (i = 0; i < diag.row; i++)
		diag.val[i] = M.val[i];


	//	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	//	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	//	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	//	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	//	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	//	free_csr_matrix(&tempA);
	create_csr_matrix(A[0].row / 2, A[0].row / 2, A[0].nnz / 2, &As[0]);
	for (i = 0; i < As[0].row; i++)
		As[0].IA[i + 1] = A[0].IA[i + 1];
	for (i = 0; i < As[0].nnz; i++)
	{
		As[0].JA[i] = A[0].JA[i];
		As[0].val[i] = A[0].val[i];
	}
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	//	interpVecP1toCR2d(&tempA, &elementDOFas, elementDOFcr, edges);
	//	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	//	free_csr_matrix(&tempA);
	//	getTransposeOfSparse(&P, &PT);


	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &A[0];
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &M; // mass matrix
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	//	aspData.R = &PT;
	//	aspData.P = &P;
	aspData.Minv = NULL;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_DiagAMGStokesNcP1_P0;

	// solver part
	solve_start = clock();

	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = minres2b(A, b, x, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Auxiliary space preconditioned MINRES costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}

	//	free_elementDOF(&elementDOFas);
	//	free_icsr_matrix(&elementdofTranas);


	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	//	free_csr_matrix(&P);
	//	free_csr_matrix(&PT);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int AbfpAsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by GMRES solver with approximate block factorization preconditioner and auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int AbfpAsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Asc, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[0];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assemblePressMassmatrixStokesNcP1P02d(&M, elements, &elementDOF[1]);
	for (i = 0; i < M.nnz; i++)
		M.val[i] *= param->precond_scale[1];
	dvector diag, diaginv;
	//	getdiag(M.row, &M, &diag);
	create_dvector(M.row, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
	{
		diag.val[i] = M.val[i];
		diaginv.val[i] = 1.0 / diag.val[i];
	}

	/************ form Shur complement Asc = A[1]A[3]^{-1}A[2] + A[0]  *********/
	dDiagVectorMultiplydCSR(&diaginv, &A[2], &tempA);
	sparseMultiplication(&A[1], &tempA, &tempB);
	free_csr_matrix(&tempA);
	sparseAddition(&tempB, &A[0], &Asc);
	free_csr_matrix(&tempB);
	/************ form Shur complement end  *********/


	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toNcP1_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &A[0], edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &A[0], elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &A[0], elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &Asc;
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &M; // mass matrix
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_AbfpAsP1Stokes;

	// solver part
	solve_start = clock();

	int restart = param->restart;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = fgmres2b(A, b, x, restart, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Block triangular preconditioned GMRES with auxiliary space method costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_dvector(&diaginv);
	free_csr_matrix(&Asc);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int Abfp2AsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by GMRES solver with approximate block factorization preconditioner and auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int Abfp2AsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[0];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assemblePressMassmatrixStokesNcP1P02d(&M, elements, &elementDOF[1]);
	for (i = 0; i < M.nnz; i++)
		M.val[i] *= param->precond_scale[1];
	dvector diag, diaginv;
	//	getdiag(M.row, &M, &diag);
	create_dvector(M.row, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
	{
		diag.val[i] = M.val[i];
		diaginv.val[i] = 1.0 / diag.val[i];
	}

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toNcP1_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &A[0], edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &A[0], elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &A[0], elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &A[0], elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &A[0];
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &M; // mass matrix
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_AbfpAsP1Stokes;

	// solver part
	solve_start = clock();

	int restart = param->restart;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = fgmres2b(A, b, x, restart, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Block triangular preconditioned GMRES with auxiliary space method costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int DiagAsP1ElasCR_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block diagonal preconditioned MINRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int DiagAsP1ElasCR_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Acr, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[1];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assembleStressMassmatrixStokesCR2d_stressP2(&M, elements, elementDOF, mu);
	//	assemblePressStressMassmatrixStokesCR2d_stressP2(&M, elements, elementDOF, mu);

	dvector diag, diaginv;
	getdiag(M.row, &M, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
		diaginv.val[i] = 1.0 / diag.val[i];

	/************ form Shur complement Acr = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dBDmat dbdM, Minv;
	assembleStressMassmatrixStokesCRdBD2d_stressP2(&dbdM, elements, elementDOF, mu);
	inverse_dBDmat(&dbdM, &Minv);
	free_dbd_matrix(&dbdM);
	dBDMultiplydCSR(&Minv, &A[1], &tempA);
	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Acr);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Acr);
	}
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toCR_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &Acr, edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &Acr, elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &Acr, elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Acr; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = &Minv;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_DiagAsP1Elas;

	// solver part
	solve_start = clock();

	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = minres2b(A, b, x, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Auxiliary space preconditioned MINRES costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Acr);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dbd_matrix(&Minv);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int DiagAsP1ElasNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block diagonal preconditioned MINRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int DiagAsP1ElasNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Acr, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[1];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assembleStressMassmatrixStokesNcP1P02d(&M, elements, &elementDOF[0], mu);

	dvector diag, diaginv;
//	getdiag(M.row, &M, &diag);
	create_dvector(M.row, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
	{
		diag.val[i] = M.val[i];
		diaginv.val[i] = 1.0 / diag.val[i];
	}

	/************ form Shur complement Acr = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dDiagVectorMultiplydCSR(&diaginv, &A[1], &tempA);
	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Acr);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Acr);
	}
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toNcP1_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &Acr, edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &Acr, elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &Acr, elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Acr; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_DiagAsP1Elas;

	// solver part
	solve_start = clock();

	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = minres2b(A, b, x, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Auxiliary space preconditioned MINRES costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Acr);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int DiagAsP1ElasDG_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block diagonal preconditioned MINRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int DiagAsP1ElasDG_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Adg, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFdg = &elementDOF[1];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	if (elementDOF[0].nfreenodes.row == 0)
		assembleStiffmatrixHuZhangA11_2d(&M, elements, nodes, elementDOF, elementdofTran, 0, mu);
	else
	{
		assembleStiffmatrixHuZhangA11_2d(&tempM, elements, nodes, elementDOF, elementdofTran, 0, mu);
		extractFreenodesMatrix11(&tempM, &M, &elementDOF[0], &elementDOF[0]);
		free_csr_matrix(&tempM);
	}

	dvector diag, diaginv;
	getdiag(M.row, &M, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
		diaginv.val[i] = 1.0 / diag.val[i];

						   /************ form Shur complement Adg = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dDiagVectorMultiplydCSR(&diaginv, &A[1], &tempA);
	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Adg);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Adg);
	}
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/
	
	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);
	
	interpVecP1toDG2d(&tempA, &elementDOFas, elementDOFdg);
	extractFreenodesMatrix1cBlock(&tempA, &P, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);


	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/



	
	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
		printf("%lf\n", diag.val[i]);
*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Adg; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_DiagAsP1Elas;

	// solver part
	solve_start = clock();

	int MaxIt = param->max_iter;
	double tol = param->tol;
//	prec = NULL;
	iter = minres2b(A, b, x, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Auxiliary space preconditioned MINRES costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Adg);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int TriAsElasMINI_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type, int print_level)
* \brief Solve Ax=b by block triangular preconditioned GMRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *aspparam pointer to ASP parameters
* \param print_level solver type
* \param print_level how much information to print out
* \return the number of iterations
*/
int TriAsElasMINI_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type, int print_level)
{
	switch (itsolver_type) {
	case 1: case 3: aspparam->mass_precond_type = 1; break;
	case 2: case 4: aspparam->mass_precond_type = 2; break;
	default:aspparam->mass_precond_type = 2;
	}

	//	if (itsolver_type == 1 || itsolver_type == 2)
	//	{
	printf("\nBlock triangular preconditioned GMRES solver with auxiliary space method\n");
	printf("Auxiliary space: P1 Lagrangian element.\n");
	if (aspparam->mass_precond_type == 1)
		printf("Mass matrix of stress is replaced by its diagonal part.\n\n");
	else
		printf("Mass matrix of stress is inversed directly.\n\n");

	TriAsP1ElasMINI_GMRES(A, b, x, aspparam, print_level);
	//	}
}

/**
* \fn int TriAsP1ElasMINI_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block triangular preconditioned GMRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int TriAsP1ElasMINI_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int variationalform_type = param->variationalform_type;
	int stress_fem_type = param->stress_fem_type;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Amini, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;
	ELEMENT_DOF *elementDOFmini = &elementDOF[2];

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;


	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	int mass_precond_type = param->mass_precond_type;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	// pressure-stress-velocity form
	if (stress_fem_type == 1) // stress discretized by piecewise P1P0 element
		assemblePressStressMassmatrixStokesMINI2d_stressP1P0(&M, elements, elementDOF, elementdofTran, mu);
	else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
		assemblePressStressMassmatrixStokesMINI2d_stressP2(&M, elements, elementDOF, elementdofTran, mu);
	else// if (stress_fem_type == 3) // stress discretized in compact form
		;//		assemblePressStressMassmatrixStokesCR2d_stressCompact(&M, elements, elementDOF, mu);

	dvector diag, diaginv;
	getdiag(M.row, &M, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
		diaginv.val[i] = 1.0 / diag.val[i];

	//	printf("%d, %d\n", elementDOF[0].dof, elementDOF[1].dof);////////////
	//	for (i = 0; i < diag.row; i++)
	//		printf("%d, %e\n", i, diag.val[i]);///////////////////////////////////
	/************ form Shur complement Amini = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dBDmat Minv;
	if (mass_precond_type == 1 || stress_fem_type == 3)
	{
		dDiagVectorMultiplydCSR(&diaginv, &A[1], &tempA);
		Minv.blk = NULL;
		Minv.row = 0;
		Minv.col = 0;
		Minv.nb = 0;
	}
	else
	{
		dBDmat dbdM;
		if (stress_fem_type == 1) // stress discretized by piecewise P1P0 element
			assemblePressStressMassmatrixStokesMINIdBD2d_stressP1P0(&dbdM, &diag, elements, elementDOF, mu);
		else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
			assemblePressStressMassmatrixStokesMINIdBD2d_stressP2(&dbdM, &diag, elements, elementDOF, mu);
		inverse_dBDmat(&dbdM, &Minv);
		free_dbd_matrix(&dbdM);
		dBDMultiplydCSR(&Minv, &A[1], &tempA);
	}

	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Amini);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Amini);
	}
	//////////////////////////////////////////////
	/*	FILE *outputFile;
	outputFile = fopen("output/matrixtest.dat", "w");
	fprintf(outputFile,"A[1]:%d, %d\n", A[1].row, A[1].col);
	for (i = 0; i < A[1].nnz; i++)
	{
	fprintf(outputFile, "(%d, %e), ", A[1].JA[i], A[1].val[i]);
	}
	fprintf(outputFile, "\n A[1] end\n");
	fprintf(outputFile, "tempA:%d, %d\n", tempA.row, tempA.col);
	for (i = 0; i < tempA.nnz; i++)
	{
	fprintf(outputFile, "(%d,%e), ", tempA.JA[i], tempA.val[i]);
	}
	fprintf(outputFile, "\n tempA end\n");
	fclose(outputFile);
	exit(1);*/
	//////////////////////////////////////////////
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toMINI_2d(&tempA, &elementDOFas, elementDOFmini);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFmini, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &Amini, elements, nodes->row, elementDOFmini);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &Amini, edges, elementDOFmini);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &Amini, elementEdge, edges, elementDOFmini);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &Amini, elements, edges, nodes->row, elementDOFmini);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &Amini, elements, nodes->row, elementDOFmini);
		else
			getSchwarzblocksVec2_vertex(&swzB, &Amini, elements, nodes->row, elementDOFmini);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	//////////////////////////////////////////////
	/*	printf("Adg:%d, %d\n", Adg.row, Adg.col);
	for (i = 0; i < Adg.nnz; i++)
	{
	printf("(%d,%e), ", Adg.JA[i], Adg.val[i]);
	}
	printf("\n Adg end\n");
	printf("As[0]:%d, %d\n", As[0].row, As[0].col);
	for (i = 0; i < As[0].nnz; i++)
	{
	//	if (fabs(A[1].val[i])>0.25)
	printf("(%d,%e), ", As[0].JA[i], As[0].val[i]);
	}
	printf("\n As[0] end\n");
	exit(1);
	printf("P:%d, %d, ", P.row, P.col);
	for (i = 0; i < P.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", P.JA[i], P.val[i]);
	}
	printf("P end\n");
	printf("PT:%d, %d, ", PT.row, PT.col);
	for (i = 0; i < PT.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", PT.JA[i], PT.val[i]);
	}
	printf("PT end\n");*/
	/////////////////////////////////////////////



	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Amini; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	if (mass_precond_type == 1 || stress_fem_type == 3)
		aspData.Minv = NULL;
	else
		aspData.Minv = &Minv;
	aspData.mass_precond_type = param->mass_precond_type;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_TriAsP1Elas;

	// solver part
	solve_start = clock();

	int restart = param->restart;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = fgmres2b(A, b, x, restart, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Block triangular preconditioned GMRES with auxiliary space method costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Amini);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dbd_matrix(&Minv);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int TriAsElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type, int print_level)
* \brief Solve Ax=b by block triangular preconditioned GMRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *aspparam pointer to ASP parameters
* \param print_level solver type
* \param print_level how much information to print out
* \return the number of iterations
*/
int TriAsElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type, int print_level)
{
	switch (itsolver_type) {
	case 1: case 3: aspparam->mass_precond_type = 1; break;
	case 2: case 4: aspparam->mass_precond_type = 2; break;
	default:aspparam->mass_precond_type = 2;
	}

	if (itsolver_type == 1 || itsolver_type == 2)
	{
		printf("\nBlock triangular preconditioned GMRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n");
		if (aspparam->mass_precond_type == 1)
			printf("Mass matrix of stress is replaced by its diagonal part.\n\n");
		else
			printf("Mass matrix of stress is inversed directly.\n\n");

		TriAsP1ElasCR_GMRES(A, b, x, aspparam, print_level);
	}
	else
	{
		printf("\nBlock triangular preconditioned GMRES solver with auxiliary space method\n");
		printf("Auxiliary space: P2 Lagrangian element.\n");
		if (aspparam->mass_precond_type == 1)
			printf("Mass matrix of stress is replaced by its diagonal part.\n\n");
		else
			printf("Mass matrix of stress is inversed directly.\n\n");

		TriAsP2ElasCR_GMRES(A, b, x, aspparam, print_level);
	}
}

/**
* \fn int TriAsP1ElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block triangular preconditioned GMRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int TriAsP1ElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int variationalform_type = param->variationalform_type;
	int stress_fem_type = param->stress_fem_type;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Acr, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;
	ELEMENT_DOF *elementDOFcr;
	if (variationalform_type == 2)
		elementDOFcr = &elementDOF[1];
	else
		elementDOFcr = &elementDOF[2];

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;


	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	int mass_precond_type = param->mass_precond_type;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	if (variationalform_type == 2) // elasticity form
	{
		if (stress_fem_type == 1) // stress discretized in compact form
			assembleStressMassmatrixStokesCR2d_stressCompact(&M, elements, elementDOF, mu);
		else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
			assembleStressMassmatrixStokesCR2d_stressP2(&M, elements, elementDOF, mu);
	}
	else // pressure-stress-velocity form
	{
		if (stress_fem_type == 1) // stress discretized in compact form
			assemblePressStressMassmatrixStokesCR2d_stressCompact(&M, elements, elementDOF, mu);
		else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
			assemblePressStressMassmatrixStokesCR2d_stressP2(&M, elements, elementDOF, mu);
	}

	dvector diag, diaginv;
	getdiag(M.row, &M, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
		diaginv.val[i] = 1.0 / diag.val[i];

	//	printf("%d, %d\n", elementDOF[0].dof, elementDOF[1].dof);////////////
	//	for (i = 0; i < diag.row; i++)
	//		printf("%d, %e\n", i, diag.val[i]);///////////////////////////////////
	/************ form Shur complement Acr = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dBDmat Minv;
	if (mass_precond_type == 1)
	{
		dDiagVectorMultiplydCSR(&diaginv, &A[1], &tempA);
		Minv.blk = NULL;
		Minv.row = 0;
		Minv.col = 0;
		Minv.nb = 0;
	}
	else
	{
		dBDmat dbdM;
		if (variationalform_type == 2) // elasticity form
		{
			if (stress_fem_type == 1) // stress discretized in compact form
				assembleStressMassmatrixStokesCRdBD2d_stressCompact(&dbdM, elements, elementDOF, mu);
			else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
				assembleStressMassmatrixStokesCRdBD2d_stressP2(&dbdM, elements, elementDOF, mu);
		}
		else // pressure-stress-velocity form
		{
			if (stress_fem_type == 1) // stress discretized in compact form
				assemblePressStressMassmatrixStokesCRdBD2d_stressCompact(&dbdM, elements, elementDOF, mu);
			else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
				assemblePressStressMassmatrixStokesCRdBD2d_stressP2(&dbdM, elements, elementDOF, mu);
		}
		inverse_dBDmat(&dbdM, &Minv);
		free_dbd_matrix(&dbdM);
		dBDMultiplydCSR(&Minv, &A[1], &tempA);
	}
	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Acr);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Acr);
	}
	//////////////////////////////////////////////
	/*	FILE *outputFile;
	outputFile = fopen("output/matrixtest.dat", "w");
	fprintf(outputFile,"A[1]:%d, %d\n", A[1].row, A[1].col);
	for (i = 0; i < A[1].nnz; i++)
	{
	fprintf(outputFile, "(%d, %e), ", A[1].JA[i], A[1].val[i]);
	}
	fprintf(outputFile, "\n A[1] end\n");
	fprintf(outputFile, "tempA:%d, %d\n", tempA.row, tempA.col);
	for (i = 0; i < tempA.nnz; i++)
	{
	fprintf(outputFile, "(%d,%e), ", tempA.JA[i], tempA.val[i]);
	}
	fprintf(outputFile, "\n tempA end\n");
	fclose(outputFile);
	exit(1);*/
	//////////////////////////////////////////////
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toCR_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &Acr, edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &Acr, elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &Acr, elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	//////////////////////////////////////////////
	/*	printf("Adg:%d, %d\n", Adg.row, Adg.col);
	for (i = 0; i < Adg.nnz; i++)
	{
	printf("(%d,%e), ", Adg.JA[i], Adg.val[i]);
	}
	printf("\n Adg end\n");
	printf("As[0]:%d, %d\n", As[0].row, As[0].col);
	for (i = 0; i < As[0].nnz; i++)
	{
	//	if (fabs(A[1].val[i])>0.25)
	printf("(%d,%e), ", As[0].JA[i], As[0].val[i]);
	}
	printf("\n As[0] end\n");
	exit(1);
	printf("P:%d, %d, ", P.row, P.col);
	for (i = 0; i < P.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", P.JA[i], P.val[i]);
	}
	printf("P end\n");
	printf("PT:%d, %d, ", PT.row, PT.col);
	for (i = 0; i < PT.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", PT.JA[i], PT.val[i]);
	}
	printf("PT end\n");*/
	/////////////////////////////////////////////



	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Acr; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	if (mass_precond_type == 1)
		aspData.Minv = NULL;
	else
		aspData.Minv = &Minv;
	aspData.mass_precond_type = param->mass_precond_type;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_TriAsP1Elas;

	// solver part
	solve_start = clock();

	int restart = param->restart;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = fgmres2b(A, b, x, restart, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Block triangular preconditioned GMRES with auxiliary space method costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Acr);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dbd_matrix(&Minv);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int TriAsP2ElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block triangular preconditioned GMRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int TriAsP2ElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int variationalform_type = param->variationalform_type;
	int stress_fem_type = param->stress_fem_type;

	int max_levels = param->max_levels + 1;
	dCSRmat As[max_levels], Acr, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF elementDOFas, elementDOFp2;
	iCSRmat elementdofTranas, elementdofTranp2;
	ELEMENT_DOF *elementDOFcr;
	if (variationalform_type == 2)
		elementDOFcr = &elementDOF[1];
	else
		elementDOFcr = &elementDOF[2];

	double lambda = param->lambda;
	double mu = param->mu;
	//	double t = param->t;


	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	int mass_precond_type = param->mass_precond_type;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 2; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 1; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	if (variationalform_type == 2) // elasticity form
	{
		if (stress_fem_type == 1) // stress discretized in compact form
			assembleStressMassmatrixStokesCR2d_stressCompact(&M, elements, elementDOF, mu);
		else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
			assembleStressMassmatrixStokesCR2d_stressP2(&M, elements, elementDOF, mu);
	}
	else // pressure-stress-velocity form
	{
		if (stress_fem_type == 1) // stress discretized in compact form
			assemblePressStressMassmatrixStokesCR2d_stressCompact(&M, elements, elementDOF, mu);
		else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
			assemblePressStressMassmatrixStokesCR2d_stressP2(&M, elements, elementDOF, mu);
	}

	dvector diag, diaginv;
	getdiag(M.row, &M, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
		diaginv.val[i] = 1.0 / diag.val[i];

	//	printf("%d, %d\n", elementDOF[0].dof, elementDOF[1].dof);////////////
	//	for (i = 0; i < diag.row; i++)
	//		printf("%d, %e\n", i, diag.val[i]);///////////////////////////////////
	/************ form Shur complement Acr = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dBDmat Minv;
	if (mass_precond_type == 1)
	{
		dDiagVectorMultiplydCSR(&diaginv, &A[1], &tempA);
		Minv.blk = NULL;
		Minv.row = 0;
		Minv.col = 0;
		Minv.nb = 0;
	}
	else
	{
		dBDmat dbdM;
		if (variationalform_type == 2) // elasticity form
		{
			if (stress_fem_type == 1) // stress discretized in compact form
				assembleStressMassmatrixStokesCRdBD2d_stressCompact(&dbdM, elements, elementDOF, mu);
			else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
				assembleStressMassmatrixStokesCRdBD2d_stressP2(&dbdM, elements, elementDOF, mu);
		}
		else // pressure-stress-velocity form
		{
			if (stress_fem_type == 1) // stress discretized in compact form
				assemblePressStressMassmatrixStokesCRdBD2d_stressCompact(&dbdM, elements, elementDOF, mu);
			else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
				assemblePressStressMassmatrixStokesCRdBD2d_stressP2(&dbdM, elements, elementDOF, mu);
		}
		inverse_dBDmat(&dbdM, &Minv);
		free_dbd_matrix(&dbdM);
		dBDMultiplydCSR(&Minv, &A[1], &tempA);
	}
	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Acr);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Acr);
	}
	//////////////////////////////////////////////
	/*	FILE *outputFile;
	outputFile = fopen("output/matrixtest.dat", "w");
	fprintf(outputFile,"A[1]:%d, %d\n", A[1].row, A[1].col);
	for (i = 0; i < A[1].nnz; i++)
	{
	fprintf(outputFile, "(%d, %e), ", A[1].JA[i], A[1].val[i]);
	}
	fprintf(outputFile, "\n A[1] end\n");
	fprintf(outputFile, "tempA:%d, %d\n", tempA.row, tempA.col);
	for (i = 0; i < tempA.nnz; i++)
	{
	fprintf(outputFile, "(%d,%e), ", tempA.JA[i], tempA.val[i]);
	}
	fprintf(outputFile, "\n tempA end\n");
	fclose(outputFile);
	exit(1);*/
	//////////////////////////////////////////////
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/

	getElementDOF_Lagrange(&elementDOFp2, elements, elementEdge, edges, nodes->row, 2);
	getTransposeOfelementDoF(&elementDOFp2, &elementdofTranp2, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFp2);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFp2, &elementdofTranp2, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFp2, &elementDOFp2);
	free_csr_matrix(&tempA);

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[1], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As + 1, Ps + 1, Rs + 1, &levelNum, &amgparam);

	interpP1toP2_2d(&tempA, &elementDOFas, &elementDOFp2, edges);
	extractFreenodesMatrix11(&tempA, &Ps[0], &elementDOFp2, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&Ps[0], &Rs[0]);

	interpVecP2toCR_2d(&tempA, &elementDOFp2, elementDOFcr);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFp2);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &Acr, edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &Acr, elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &Acr, elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	//////////////////////////////////////////////
	/*	printf("Adg:%d, %d\n", Adg.row, Adg.col);
	for (i = 0; i < Adg.nnz; i++)
	{
	printf("(%d,%e), ", Adg.JA[i], Adg.val[i]);
	}
	printf("\n Adg end\n");
	printf("As[0]:%d, %d\n", As[0].row, As[0].col);
	for (i = 0; i < As[0].nnz; i++)
	{
	//	if (fabs(A[1].val[i])>0.25)
	printf("(%d,%e), ", As[0].JA[i], As[0].val[i]);
	}
	printf("\n As[0] end\n");
	exit(1);
	printf("P:%d, %d, ", P.row, P.col);
	for (i = 0; i < P.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", P.JA[i], P.val[i]);
	}
	printf("P end\n");
	printf("PT:%d, %d, ", PT.row, PT.col);
	for (i = 0; i < PT.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", PT.JA[i], PT.val[i]);
	}
	printf("PT end\n");*/
	/////////////////////////////////////////////



	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum + 1;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Acr; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	if (mass_precond_type == 1)
		aspData.Minv = NULL;
	else
		aspData.Minv = &Minv;
	aspData.mass_precond_type = param->mass_precond_type;
	aspData.swzB = &swzB;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_TriAsP1Elas;

	// solver part
	solve_start = clock();

	int restart = param->restart;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = fgmres2b(A, b, x, restart, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Block triangular preconditioned GMRES with auxiliary space method costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum + 1; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);
	free_elementDOF(&elementDOFp2);
	free_icsr_matrix(&elementdofTranp2);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Acr);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dbd_matrix(&Minv);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int TriAsP1ElasNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block triangular preconditioned GMRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int TriAsP1ElasNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Acr, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFcr = &elementDOF[1];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;
	double nu = param->nu;
	int problem_num = param->problem_num;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	assembleStressMassmatrixStokesNcP1P02d(&M, elements, elementDOF, mu);

	dvector diag, diaginv;
	create_dvector(M.row, &diag);
	for (i = 0; i < diag.row; i++)
		diag.val[i] = M.val[i];
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
		diaginv.val[i] = 1.0 / diag.val[i];

	/************ form Shur complement Acr = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dDiagVectorMultiplydCSR(&diaginv, &A[1], &tempA);
	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Acr);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Acr);
	}
	//////////////////////////////////////////////
	/*	FILE *outputFile;
	outputFile = fopen("output/matrixtest.dat", "w");
	fprintf(outputFile,"A[1]:%d, %d\n", A[1].row, A[1].col);
	for (i = 0; i < A[1].nnz; i++)
	{
	fprintf(outputFile, "(%d, %e), ", A[1].JA[i], A[1].val[i]);
	}
	fprintf(outputFile, "\n A[1] end\n");
	fprintf(outputFile, "tempA:%d, %d\n", tempA.row, tempA.col);
	for (i = 0; i < tempA.nnz; i++)
	{
	fprintf(outputFile, "(%d,%e), ", tempA.JA[i], tempA.val[i]);
	}
	fprintf(outputFile, "\n tempA end\n");
	fclose(outputFile);
	exit(1);*/
	//////////////////////////////////////////////
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toNcP1_2d(&tempA, &elementDOFas, elementDOFcr, edges);
	extractFreenodesMatrix11cBlock(&tempA, &P, elementDOFcr, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	dOBDmat swzB;
	if (param->smoother == MSWZ || param->smoother == ASWZ || param->smoother == SMSWZ)
	{
		if (param->schwarz_type == 1)
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 2)
			getSchwarzblocksVec2_edge(&swzB, &Acr, edges, elementDOFcr);
		else if (param->schwarz_type == 3)
			getSchwarzblocksVec2_element(&swzB, &Acr, elementEdge, edges, elementDOFcr);
		else if (param->schwarz_type == 4)
			getSchwarzblocksVec2_edgevertex(&swzB, &Acr, elements, edges, nodes->row, elementDOFcr);
		else if (param->schwarz_type == 5)
			getSchwarzblocksVec2_elementvertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);
		else
			getSchwarzblocksVec2_vertex(&swzB, &Acr, elements, nodes->row, elementDOFcr);

	}
	else
	{
		swzB.row = 0;
		swzB.col = 0;
		swzB.nb = 0;
		swzB.blk = NULL;
		swzB.rindices = NULL;
		swzB.cindices = NULL;
	}

	//////////////////////////////////////////////
	/*	printf("Adg:%d, %d\n", Adg.row, Adg.col);
	for (i = 0; i < Adg.nnz; i++)
	{
	printf("(%d,%e), ", Adg.JA[i], Adg.val[i]);
	}
	printf("\n Adg end\n");
	printf("As[0]:%d, %d\n", As[0].row, As[0].col);
	for (i = 0; i < As[0].nnz; i++)
	{
	//	if (fabs(A[1].val[i])>0.25)
	printf("(%d,%e), ", As[0].JA[i], As[0].val[i]);
	}
	printf("\n As[0] end\n");
	exit(1);
	printf("P:%d, %d, ", P.row, P.col);
	for (i = 0; i < P.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", P.JA[i], P.val[i]);
	}
	printf("P end\n");
	printf("PT:%d, %d, ", PT.row, PT.col);
	for (i = 0; i < PT.nnz; i++)
	{
	//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", PT.JA[i], PT.val[i]);
	}
	printf("PT end\n");*/
	/////////////////////////////////////////////



	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Acr; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;
	aspData.swzB = &swzB;
	
	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_TriAsP1Elas;

	// solver part
	solve_start = clock();

	int restart = param->restart;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = fgmres2b(A, b, x, restart, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Block triangular preconditioned GMRES with auxiliary space method costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Acr);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);
	free_dobd_matrix(&swzB);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
* \fn int TriAsP1ElasDG_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
* \brief Solve Ax=b by block triangular preconditioned GMRES solver with auxiliary space method
* \param *A	pointer to the dCSRmat matrix
* \param *b	pointer to the dvector of right hand side
* \param *x	pointer to the dvector of dofs
* \param *param pointer to ASP parameters
* \param print_level how much information to print out
* \return the number of iterations
*/
int TriAsP1ElasDG_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level)
{
	int levelNum = 1;
	int iter, i;

	int max_levels = param->max_levels;
	dCSRmat As[max_levels], Adg, tempA, tempB, M, tempM;
	dCSRmat P, PT, Rs[max_levels], Ps[max_levels];
	ELEMENT *elements = param->elements;
	idenmat *elementEdge = param->elementEdge;
	EDGE *edges = param->edges;
	dennode *nodes = param->nodes;
	iCSRmat *edgesTran = param->edgesTran;
	ivector *nodeCEdge = param->nodeCEdge;
	ELEMENT_DOF *elementDOF = param->elementDOF;
	iCSRmat *elementdofTran = param->elementdofTran;
	ELEMENT_DOF *elementDOFdg = &elementDOF[1];
	ELEMENT_DOF elementDOFas;
	iCSRmat elementdofTranas;

	double lambda = param->lambda;
	double mu = param->mu;

	AMG_param amgparam; /* parameters for AMG */
	amgparam.max_levels = param->max_levels;
	amgparam.coarsening_type = param->AMG_coarsening_type;
	amgparam.interpolation_type = param->AMG_interpolation_type;
	amgparam.coarse_dof = param->AMG_coarse_dof;
	amgparam.strong_threshold = param->AMG_strong_threshold;
	amgparam.truncation_threshold = param->AMG_truncation_threshold;
	amgparam.max_row_sum = param->AMG_max_row_sum;
	amgparam.print_level = 1;

	clock_t solve_start, solve_end;

	// initialize
	for (i = 1; i<max_levels; i++)
	{
		As[i].row = 0;
		As[i].col = 0;
		As[i].IA = NULL;
		As[i].JA = NULL;
		As[i].val = NULL;
	}
	for (i = 0; i<max_levels; i++)
	{
		Ps[i].row = 0;
		Ps[i].col = 0;
		Ps[i].IA = NULL;
		Ps[i].JA = NULL;
		Ps[i].val = NULL;

		Rs[i].row = 0;
		Rs[i].col = 0;
		Rs[i].IA = NULL;
		Rs[i].JA = NULL;
		Rs[i].val = NULL;
	}

	// setup preconditioner
	if (elementDOF[0].nfreenodes.row == 0)
		assembleStiffmatrixHuZhangA11_2d(&M, elements, nodes, elementDOF, elementdofTran, 0, mu);
	else
	{
		assembleStiffmatrixHuZhangA11_2d(&tempM, elements, nodes, elementDOF, elementdofTran, 0, mu);
		extractFreenodesMatrix11(&tempM, &M, &elementDOF[0], &elementDOF[0]);
		free_csr_matrix(&tempM);
	}

	dvector diag, diaginv;
	getdiag(M.row, &M, &diag);
	create_dvector(diag.row, &diaginv);
	for (i = 0; i < diag.row; i++)
		diaginv.val[i] = 1.0 / diag.val[i];

	/************ form Shur complement Adg = A[2]A[0]^{-1}A[1] - A[3]  *********/
	dDiagVectorMultiplydCSR(&diaginv, &A[1], &tempA);
	if (A[3].row>0)
	{
		sparseMultiplication(&A[2], &tempA, &tempB);
		sparseSubtraction(&tempB, &A[3], &Adg);
		free_csr_matrix(&tempB);
	}
	else
	{
		sparseMultiplication(&A[2], &tempA, &Adg);
	}
	//////////////////////////////////////////////
/*	FILE *outputFile;
	outputFile = fopen("output/matrixtest.dat", "w");
	fprintf(outputFile,"A[1]:%d, %d\n", A[1].row, A[1].col);
	for (i = 0; i < A[1].nnz; i++)
	{
		fprintf(outputFile, "(%d, %e), ", A[1].JA[i], A[1].val[i]);
	}
	fprintf(outputFile, "\n A[1] end\n");
	fprintf(outputFile, "tempA:%d, %d\n", tempA.row, tempA.col);
	for (i = 0; i < tempA.nnz; i++)
	{
		fprintf(outputFile, "(%d,%e), ", tempA.JA[i], tempA.val[i]);
	}
	fprintf(outputFile, "\n tempA end\n");
	fclose(outputFile);
	exit(1);*/
	//////////////////////////////////////////////
	free_csr_matrix(&tempA);
	/************ form Shur complement end  *********/

	getElementDOF_Lagrange(&elementDOFas, elements, elementEdge, edges, nodes->row, 1);
	getTransposeOfelementDoF(&elementDOFas, &elementdofTranas, 0);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFas);
	assembleStiffmatrixLagrange(&tempA, elements, elementEdge, edges, nodes, &elementDOFas, &elementdofTranas, mu);
	extractFreenodesMatrix11(&tempA, &As[0], &elementDOFas, &elementDOFas);
	free_csr_matrix(&tempA);
	classicAMG_setup(As, Ps, Rs, &levelNum, &amgparam);

	interpVecP1toDG2d(&tempA, &elementDOFas, elementDOFdg);
	extractFreenodesMatrix1cBlock(&tempA, &P, &elementDOFas);
	free_csr_matrix(&tempA);
	getTransposeOfSparse(&P, &PT);

	//////////////////////////////////////////////
/*	printf("Adg:%d, %d\n", Adg.row, Adg.col);
	for (i = 0; i < Adg.nnz; i++)
	{
	printf("(%d,%e), ", Adg.JA[i], Adg.val[i]);
	}
	printf("\n Adg end\n");
	printf("As[0]:%d, %d\n", As[0].row, As[0].col);
	for (i = 0; i < As[0].nnz; i++)
	{
//	if (fabs(A[1].val[i])>0.25)
	printf("(%d,%e), ", As[0].JA[i], As[0].val[i]);
	}
	printf("\n As[0] end\n");
	exit(1);
	printf("P:%d, %d, ", P.row, P.col);
	for (i = 0; i < P.nnz; i++)
	{
//	if (fabs(A[2].val[i])>0.25)
	printf("(%d,%e),", P.JA[i], P.val[i]);
	}
	printf("P end\n");
	printf("PT:%d, %d, ", PT.row, PT.col);
	for (i = 0; i < PT.nnz; i++)
	{
		//	if (fabs(A[2].val[i])>0.25)
		printf("(%d,%e),", PT.JA[i], PT.val[i]);
	}
	printf("PT end\n");*/
	/////////////////////////////////////////////



	/**************************************************
	dvector uh;
	create_dvector(b[1].row, &uh);
	param->elementDOF = &(param->elementDOF[1]);
	printf("Auxiliary space preconditioned CG solver\n");
	asP1ElasDG_PCG(&Adg, &b[1], &uh, param, print_level);

	return 0;
	/**************************************************************/




	/*printf("diag:\n");
	for (i = 0; i < diag.row; i++)
	printf("%lf\n", diag.val[i]);
	*/
	precond_data aspData;
	aspData.max_levels = levelNum;
	aspData.max_iter = param->mg_max_iter;
	aspData.tol = param->mg_tol;
	aspData.precond_type = param->precond_type;
	aspData.precond_scale = param->precond_scale;
	aspData.smoother = param->smoother;
	aspData.schwarz_type = param->schwarz_type;
	aspData.smooth_iter = param->smooth_iter;
	aspData.mg_smoother = param->mg_smoother;
	aspData.mg_smooth_iter = param->mg_smooth_iter;
	aspData.diag = &diag;
	aspData.precA[0] = &M;  // mass matrix
	aspData.precA[1] = &A[1];
	aspData.precA[2] = &A[2];
	aspData.precA[3] = &Adg; // Schur complement
	aspData.As = As;
	aspData.Rs = Rs;
	aspData.Ps = Ps;
	aspData.R = &PT;
	aspData.P = &P;
	aspData.Minv = NULL;

	precond *prec = (precond *)malloc(sizeof(precond));
	prec->data = &aspData;
	prec->fct_dvec = precond_TriAsP1Elas;

	// solver part
	solve_start = clock();

	int restart = param->restart;
	int MaxIt = param->max_iter;
	double tol = param->tol;
	//	prec = NULL;
	iter = fgmres2b(A, b, x, restart, MaxIt, tol, prec, print_level);

	solve_end = clock();

	double solve_duration = (double)(solve_end - solve_start) / (double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Block triangular preconditioned GMRES with auxiliary space method costs %f seconds.\n", solve_duration);
		printf("Number of iterations = %d.\n", iter);
	}

	for (i = 0; i<levelNum; i++)
	{
		free_csr_matrix(&As[i]);
		free_csr_matrix(&Rs[i]);
		free_csr_matrix(&Ps[i]);
	}
	free_elementDOF(&elementDOFas);
	free_icsr_matrix(&elementdofTranas);

	//	free_elementDOF(&elementDOFdg);
	free_csr_matrix(&M);
	free_dvector(&diag);
	free_csr_matrix(&Adg);
	free_csr_matrix(&P);
	free_csr_matrix(&PT);

	//	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
 * \fn int classicAMG_GMRES(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level)
 * \brief Solve Ax=b by preconditioned generalized minimum residual method (with restarts) (GMRES), 
 * with Ruge and Stuben's classic AMG as precondition
 * \param *A	pointer to the dCSRmat matrix
 * \param *b	pointer to the dvector of right hand side
 * \param *x	pointer to the dvector of dofs
 * \param *param pointer to AMG parameters
 * \param print_level how much information to print out
 * \return the number of iterations
 */
int classicAMG_GMRES(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level)
{
	int max_levels=param->max_levels;
	int levelNum=1, iter, i, l;
	
	dCSRmat AA[max_levels];
	dCSRmat P[max_levels], PT[max_levels];	
	
	clock_t classicAMG_start, classicAMG_end;
	
	// initialize
	for(i=0;i<max_levels;i++)
	{
		AA[i].row=0;
		AA[i].col=0;
		AA[i].IA=NULL;
		AA[i].JA=NULL;
		AA[i].val=NULL;
		
		P[i].row=0;
		P[i].col=0;
		P[i].IA=NULL;
		P[i].JA=NULL;
		P[i].val=NULL;
		
		PT[i].row=0;
		PT[i].col=0;
		PT[i].IA=NULL;
		PT[i].JA=NULL;
		PT[i].val=NULL;
	}
	
	AA[0].row=A->row;
	AA[0].col=A->col;
	AA[0].IA=(int*)calloc(AA[0].row+1, sizeof(int));
	for(i=0;i<=AA[0].row;i++) AA[0].IA[i]=A->IA[i];
	AA[0].JA=(int*)calloc(AA[0].IA[AA[0].row], sizeof(int));
	AA[0].val=(double*)calloc(AA[0].IA[AA[0].row], sizeof(double));
	for(i=0;i<AA[0].IA[AA[0].row];i++) {
		AA[0].JA[i]=A->JA[i];
		AA[0].val[i]=A->val[i];
	}
	
	// setup preconditioner
	classicAMG_setup(AA, P, PT, &levelNum, param);
	
	precond_data amgData;
	amgData.max_levels = levelNum;
	amgData.max_iter = param->AMG_max_iter;
	amgData.tol      = param->AMG_tol;
	amgData.smoother = param->smoother;
	amgData.presmooth_iter  = param->presmooth_iter;
	amgData.postsmooth_iter = param->postsmooth_iter;
	amgData.postsmooth_iter = param->coarsening_type;
	amgData.Aarray    = (dCSRmat **)malloc(3*sizeof(dCSRmat*));
	amgData.Aarray[0] = AA;
	amgData.Aarray[1] = PT;
	amgData.Aarray[2] = P;
	
	precond *prec = (precond *)malloc(sizeof(precond));	
	prec->data = &amgData;
	prec->fct = precond_classicAMG;
	
	// solver part
	classicAMG_start=clock();
	
	int restart=param->restart;
	int MaxIt=param->max_iter;
	double tol=param->tol;
	iter=gmres(A, b, x, restart, MaxIt, tol, prec, print_level);		
	
	classicAMG_end=clock();
	
	double classicAMGduration = (double)(classicAMG_end - classicAMG_start)/(double)(CLOCKS_PER_SEC);

	if (print_level>0) {
		printf("Classic AMG preconditioned GMRES costs %f seconds.\n", classicAMGduration);
		printf("Number of iterations = %d.\n", iter);
	}
	
	for(l=0;l<levelNum;l++)
	{
		free(AA[l].IA);
		free(AA[l].JA);
		free(AA[l].val);
	}
	
	for(l=0;l<levelNum-1;l++)
	{
		free(P[l].IA);
		free(P[l].JA);
		free(P[l].val);
		free(PT[l].IA);
		free(PT[l].JA);
		free(PT[l].val);
	}

	free(amgData.Aarray);
	free(prec);

	return iter;
}

/**
 * \fn int cg_pcg(dCSRmat *A, dvector *b, dvector *sigma, dvector *u, int MaxIt, double tol, AMG_param *param, int itsolver_type, int print_level)
 *	 \brief A preconditioned conjugate gradient (CG) method for solving Au=b 
 *	 \param *A	 pointer to the coefficient matrix
 *	 \param *b	 pointer to the dvector of right hand side
 *	 \param *u	 pointer to the dvector of DOFs
 *	 \param MaxIt integer, maximal number of iterations
 *	 \param tol double float, the tolerance for stopage
 *	 \param *pre pointer to the structure of precondition (precond) 
 *   \param print_level how much information to print out
 *	 \return the number of iterations
 */
int cg_pcg(dCSRmat *A, dvector *b, dvector *sigma, dvector *u, int MaxIt, double tol, AMG_param *param, int itsolver_type, int print_level)
{
	int iter = 0, m = A[2].row;
	double alpha, beta, error, temp1, temp2, tempb;
	double *p, *z, *r, *t;
	dvector f1, f2;

	p = (double *)calloc(m, sizeof(double));
	z = (double *)calloc(m, sizeof(double));
	r = (double *)calloc(m, sizeof(double));
	t = (double *)calloc(m, sizeof(double));
	create_dvector(A[0].row, &f1);
	create_dvector(A[0].row, &f2);

	// (b,b)
	tempb = dot_array(m, b->val, b->val);

	// r = b-A*u
	copy_array(m, b->val, r);
	//	sparse_mv(-1.0,&A[3],u->val,r);
	//	sparse_mv0(1.0,&A[1],u->val,f1.val);
	//	itsolverInCG(&A[0], &f1, &f2, param, itsolver_type, 1000, 1e-8, 1);
	//	sparse_mv(1.0,&A[2],f2.val,r);

	// z = B*r
	copy_array(m, r, z); /* No preconditioner, B=I */

						 // p = z
	copy_array(m, z, p);

	// temp1=(z_{k-1},r_{k-1})
	temp1 = dot_array(m, z, r);

	while (iter<MaxIt)
	{
		iter++;

		// t=A*p
		init_array(m, t, 0.0);
		//		sparse_mv(1.0,&A[3],p,t);
		sparse_mv0(1.0, &A[1], p, f1.val);
		itsolverInCG(&A[0], &f1, &f2, param, itsolver_type, 1000, 1e-8, 0);
		sparse_mv(-1.0, &A[2], f2.val, t);

		// comupte alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
		temp2 = dot_array(m, t, p);
		alpha = temp1 / temp2;

		// compute u_k=u_{k-1} + alpha_k*p_{k-1}
		axpy_array(m, alpha, p, u->val);

		// compute r_k=r_{k-1} - alpha_k*A*p_{k-1}
		//		sparse_mv(-alpha,&A[3],p,r);
		sparse_mv0(1.0, &A[1], p, f1.val);
		itsolverInCG(&A[0], &f1, &f2, param, itsolver_type, 1000, 1e-8, 0);
		sparse_mv(alpha, &A[2], f2.val, r);

		temp2 = dot_array(m, r, r);
		if (temp2<1e-30) {
			iter = iter*(-1) - 1;
			break;
		}

		// relative residual = ||b-Au||_2/||b||_2=||r||_2/||b||_2
		error = sqrt(temp2 / tempb);
		if (print_level>1)
			printf("Iteration %3d: relative residual = %e\n", iter, error);
		if (error<tol) break;

		// z_k = B*r_k
		copy_array(m, r, z);	 /* No preconditioner, B=I */

								 // compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
		temp2 = dot_array(m, z, r);
		beta = temp2 / temp1;
		temp1 = temp2;

		// compute p_k = z_k + beta_k*p_{k-1}
		axpby_array(m, 1.0, z, beta, p);
	}

	sparse_mv0(-1.0, &A[1], u->val, f1.val);
	itsolverInCG(&A[0], &f1, sigma, param, itsolver_type, 1000, 1e-10, print_level);


	if (iter<0) {
		iter = iter*(-1) - 1;
		free(p);
		free(r);
		free(z);
		free(t);
		free_dvector(&f1);
		free_dvector(&f2);
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
	free_dvector(&f1);
	free_dvector(&f2);

	return iter;
}

int itsolverInCG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int itsolver_type, int itsolver_maxit, double itsolver_tol, int print_level)
{

/*	FILE *meshFile;
	int i,j;
	meshFile=fopen("output/A.dat", "w");
	for(i=0;i<A->row;i++)
	{
		for(j=A->IA[i];j<A->IA[i+1];j++)
			fprintf(meshFile,"%d %d %lf\n",i+1,A->JA[j]+1,A->val[j]);
	}
	fclose(meshFile);
	meshFile=fopen("output/b.dat", "w");
	for(i=0;i<b->row;i++)
	{
		fprintf(meshFile,"%lf\n",b->val[i]);
	}
	fclose(meshFile);
	exit(1);
	*/////////////////////////////////////////////////////////////////////////////////////

	/* AMG solver */	
	if (itsolver_type == 1) {	
//		printf("AMG iterative solver\n");		
		classicAMG(A, b, x, param);		
	}
	
	/* PCG+AMG */
	else if (itsolver_type == 2) {		
//		printf("AMG preconditioned CG solver\n");	
		classicAMG_PCG(A, b, x, param, print_level);
	}
	
	/* PCG+diag */
	else if (itsolver_type == 3) {
//		printf("Diagonal preconditioned CG solver\n");			
		diag_PCG(A, b, x, itsolver_maxit, itsolver_tol, print_level);
	}	
	
	/* CG */
	else if (itsolver_type == 4) {
//		printf("Classical CG solver\n");			
		standard_CG(A, b, x, itsolver_maxit, itsolver_tol, print_level);		
	}
	
	/* GMRES+AMG */
	else if (itsolver_type == 5) {		
		
//		printf("AMG preconditioned GMRES solver\n");	
		classicAMG_GMRES(A, b, x, param, print_level);
	}

	/* GMRES */
	else if (itsolver_type == 6) {
//		printf("Classical GMRES solver\n");			
		gmres(A, b, x, param->restart, itsolver_maxit, itsolver_tol, NULL, print_level);
	}
}