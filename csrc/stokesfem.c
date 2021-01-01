/*
 *  fem.c
 *
 *  Created by Xuehai Huang on 05/21/17.
 *  Copyright 2017 WZU. All rights reserved.
 *
 */

/*! \file assemble.c
 *  \brief Assembling for stiffness matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn void stokesfem(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
 * \brief finite element methods for Stokes equation
 * \param *uh pointer to solution
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                       the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
                                       the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
 * \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \param Input input data structure
 * \return void
 */
void stokesfem(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{	
	int variationalform_type = Input->variationalform_type;
	if(variationalform_type==1) // original stokes form
	{ 
		elementDOF[2].val = NULL;
		elementDOF[2].row = 0;
		elementDOF[2].nfFlag.row = 0;
		elementDOF[2].freenodes.row = 0;
		elementDOF[2].nfreenodes.row = 0;
		elementDOF[2].index.row = 0;

		printf("Discretize Stokes equation in original stokes form.\n");
		stokesNcP1P0(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else if (variationalform_type == 2) // elasticity form
	{
		elementDOF[2].val = NULL;
		elementDOF[2].row = 0;
		elementDOF[2].nfFlag.row = 0;
		elementDOF[2].freenodes.row = 0;
		elementDOF[2].nfreenodes.row = 0;
		elementDOF[2].index.row = 0;

		printf("Discretize Stokes equation in elasticity form.\n");
		stokesfem_elas(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else // pressure-stress-velocity form
	{
		printf("The solver of the pressure-stress-velocity form doesn't exist for noconforming P1-P0 element.\n");
		exit(0);
//		printf("Discretize Stokes equation in pressure-stress-velocity form.\n");
//		stokesfem_psv(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
}

/**
* \fn void stokesfem_elas(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief finite element methods for Stokes equation in elasticity form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesfem_elas(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	stokesNcP1P0_elas(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);

/*	int stress_fem_type = Input->stress_fem_type;
	if (stress_fem_type == 1) // stress discretized in compact form
	{
		printf("Stress is discretized in compact form.\n");
		stokesCR_elas_stressCompact(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
	{
		printf("Stress is discretized by piecewise P2 element.\n");
		stokesCR_elas_stressP2(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}*/
}

/**
* \fn void stokesfem_psv(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief finite element methods for Stokes equation in pressure-stress-velocity form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesfem_psv(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	int stress_fem_type = Input->stress_fem_type;
	if (stress_fem_type == 1) // stress discretized in compact form
	{
		printf("Stress is discretized in compact form.\n\n");
		stokesCR_psv_stressCompact(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else if (stress_fem_type == 2) // stress discretized by piecewise P2 element
	{
		printf("Stress is discretized by piecewise P2 element.\n\n");
		stokesCR_psv_stressP2(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else// if (stress_fem_type == 3) // stress discretized in compact form
	{
		printf("stress_fem_type is 1 or 2.\n\n");
		exit(0);
	}
}

/**
* \fn void stokesNcP1P0(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief Crouzeix每Raviart element method for Stokes equation in original form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesNcP1P0(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// velocity u
	getElementDOF_NoncfmP1(&elementDOF[0], elementEdge, edges->row);
	// press p
	getElementDOF(&elementDOF[1], elements->row, 0);

	getFreenodesInfoNcP1Vector2(edges, &elementDOF[0]);

	getTransposeOfelementDoF(&elementDOF[0], &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesNcP1P0(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	if (itsolver_type == 1)
	{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1StokesNcP1P0_MINRES(A, b, _uh, &aspparam, print_level);
	}
	else
	{
		printf("\nAMG Block diagonal preconditioned MINRES solver with auxiliary space method.\n\n");
		DiagAMGStokesNcP1P0_MINRES(A, b, _uh, &aspparam, print_level); // sometimes better
	//	DiagAsP1StokesNcP1P0_MINRES(A, b, _uh, &aspparam, print_level);
	}

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[elementDOF[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[i] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}

/**
* \fn void stokesNcP1P0_elas(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief nonconforming P1每P0 element method for Stokes equation in elasticity form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesNcP1P0_elas(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// stress sigma
	getElementDOF(&elementDOF[0], elements->row, 0);
	// velocity u
	getElementDOF_NoncfmP1(&elementDOF[1], elementEdge, edges->row);

	getFreenodesInfoNcP1Vector2(edges, &elementDOF[1]);

	getTransposeOfelementDoF(&elementDOF[1], &elementdofTran, 0);
	
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesNcP1P0_elas(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	if (itsolver_type == 1)
	{
		printf("Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		DiagAsP1ElasNcP1P0_MINRES(A, b, _uh, &aspparam, print_level);
	}
	else
	{
		printf("Block triangular preconditioned GMRES solver with auxiliary space method\n");
		TriAsP1ElasNcP1P0_GMRES(A, b, _uh, &aspparam, print_level);
	}
	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[i] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[elementDOF[1].freenodes.val[i]] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}

/**
* \fn void stokesMINI(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief MINI element method for Stokes equation in original form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesMINI(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran[2];

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// velocity u
	getElementDOF_MINI(&elementDOF[0], elements, nodes->row);
	// press p
	getElementDOF_Lagrange(&elementDOF[1], elements, elementEdge, edges, nodes->row, 1);

	getFreenodesInfoMINIVector2(nodes, &elementDOF[0]);

	getTransposeOfelementDoF(&elementDOF[0], &elementdofTran[0], 0);
	getTransposeOfelementDoF(&elementDOF[1], &elementdofTran[1], 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesMINI(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	printf("\nBlock diagonal preconditioned MINRES solver with auxiliary space method\n");
	printf("Auxiliary space: P1 Lagrangian element.\n\n");
	DiagAsP1StokesMINI_MINRES(A, b, _uh, &aspparam, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[elementDOF[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[i] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran[0].IA);
	free(elementdofTran[0].JA);
	free(elementdofTran[1].IA);
	free(elementdofTran[1].JA);
}

/**
* \fn void stokesMINI_psv_stressP1P0(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief MINI element method for Stokes equation in pressure-stress-velocity form with stress discretized by piecewise P2 element
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesMINI_psv_stressP1P0(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran[2];

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// press p
	getElementDOF_Lagrange(&elementDOF[0], elements, elementEdge, edges, nodes->row, 1);
	// stress sigma
	getElementDOF_MINIsymTensorP1P0(&elementDOF[1], elements->row);
	// velocity u
	getElementDOF_MINI(&elementDOF[2], elements, nodes->row);

	getFreenodesInfoMINIVector2(nodes, &elementDOF[2]);

	getTransposeOfelementDoF(&elementDOF[0], &elementdofTran[0], 0);
	getTransposeOfelementDoF(&elementDOF[2], &elementdofTran[1], 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesMINI_psv_stressP1P0(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	TriAsElasMINI_GMRES(A, b, _uh, &aspparam, itsolver_type, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[i] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[elementDOF[2].freenodes.val[i]] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran[0].IA);
	free(elementdofTran[0].JA);
	free(elementdofTran[1].IA);
	free(elementdofTran[1].JA);
}

/**
* \fn void stokesMINI_psv_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief MINI element method for Stokes equation in pressure-stress-velocity form with stress discretized by piecewise P2 element
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesMINI_psv_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran[2];

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// press p
	getElementDOF_Lagrange(&elementDOF[0], elements, elementEdge, edges, nodes->row, 1);
	// stress sigma
	getElementDOF_symTensor(&elementDOF[1], elements->row, 2);
	// velocity u
	getElementDOF_MINI(&elementDOF[2], elements, nodes->row);

	getFreenodesInfoMINIVector2(nodes, &elementDOF[2]);

	getTransposeOfelementDoF(&elementDOF[0], &elementdofTran[0], 0);
	getTransposeOfelementDoF(&elementDOF[2], &elementdofTran[1], 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesMINI_psv_stressP2(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	TriAsElasMINI_GMRES(A, b, _uh, &aspparam, itsolver_type, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[i] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[elementDOF[2].freenodes.val[i]] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran[0].IA);
	free(elementdofTran[0].JA);
	free(elementdofTran[1].IA);
	free(elementdofTran[1].JA);
}

/**
* \fn void stokesCR(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief Crouzeix每Raviart element method for Stokes equation in original form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesCR(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// velocity u
	getElementDOF_CR(&elementDOF[0], elements, elementEdge, edges->row, nodes->row);
	// press p
	getElementDOF(&elementDOF[1], elements->row, 1);

	getFreenodesInfoCRVector2(edges, nodes, &elementDOF[0]);

	getTransposeOfelementDoF(&elementDOF[0], &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesCR(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, mu, t);
	
	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	printf("\nBlock diagonal preconditioned MINRES solver with auxiliary space method\n");
	printf("Auxiliary space: P1 Lagrangian element.\n\n");
	DiagAsP1StokesCR_MINRES(A, b, _uh, &aspparam, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[elementDOF[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[i] = _uh[1].val[i];
	
	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}

/**
* \fn void stokesCR_elas_stressCompact(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief Crouzeix每Raviart element method for Stokes equation in elasticity with stress discretized in compact form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesCR_elas_stressCompact(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// press p
	getElementDOF(&elementDOF[0], elements->row, 1);
	// stress sigma
	getElementDOF_CRsymTensor(&elementDOF[1], elements->row);
	// velocity u
	getElementDOF_CR(&elementDOF[2], elements, elementEdge, edges->row, nodes->row);

	getFreenodesInfoCRVector2(edges, nodes, &elementDOF[2]);

	getTransposeOfelementDoF(&elementDOF[2], &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesCR_elas_stressCompact(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF + 1;//////////////////////////////////////////////////////////////////
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	TriAsElasCR_GMRES(A, b, _uh, &aspparam, itsolver_type, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[i] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[elementDOF[2].freenodes.val[i]] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}

/**
* \fn void stokesCR_elas_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief Crouzeix每Raviart element method for Stokes equation in elasticity with stress discretized by piecewise P2 element
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesCR_elas_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// press p
	getElementDOF(&elementDOF[0], elements->row, 1);
	// stress sigma
	getElementDOF_symTensor(&elementDOF[1], elements->row, 2);
	// velocity u
	getElementDOF_CR(&elementDOF[2], elements, elementEdge, edges->row, nodes->row);

	getFreenodesInfoCRVector2(edges, nodes, &elementDOF[2]);

	getTransposeOfelementDoF(&elementDOF[2], &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesCR_elas_stressP2(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF + 1;//////////////////////////////////////////////////////////////////
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	TriAsElasCR_GMRES(A, b, _uh, &aspparam, itsolver_type, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[i] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[elementDOF[2].freenodes.val[i]] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}

/**
* \fn void stokesCR_psv_stressCompact(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief Crouzeix每Raviart element method for Stokes equation in pressure-stress-velocity form with stress discretized  in compact form
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesCR_psv_stressCompact(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// press p
	getElementDOF(&elementDOF[0], elements->row, 1);
	// stress sigma
	getElementDOF_CRsymTensor(&elementDOF[1], elements->row);
	// velocity u
	getElementDOF_CR(&elementDOF[2], elements, elementEdge, edges->row, nodes->row);

	getFreenodesInfoCRVector2(edges, nodes, &elementDOF[2]);

	getTransposeOfelementDoF(&elementDOF[2], &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesCR_psv_stressCompact(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	TriAsElasCR_GMRES(A, b, _uh, &aspparam, itsolver_type, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[i] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[elementDOF[2].freenodes.val[i]] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}

/**
* \fn void stokesCR_psv_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief Crouzeix每Raviart element method for Stokes equation in pressure-stress-velocity form with stress discretized by piecewise P2 element
* \param *uh pointer to solution
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param Input input data structure
* \return void
*/
void stokesCR_psv_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A[4];
	dvector b[2], _uh[2];
	//	ELEMENT_DOF elementDOF[3];
	iCSRmat elementdofTran;

	int print_level = Input->print_level;
	int problem_num = Input->problem_num;
	int itsolver_type = Input->itsolver_type;
	int itsolver_maxit = Input->itsolver_maxit;
	double itsolver_tol = Input->itsolver_tol;

	double nu = Input->nu;
	double lambda = Input->lambda;
	double mu = Input->mu;
	double t = Input->t;
	int glevelNum = Input->glevelNum;
	int CglevelNum = Input->CglevelNum;
	int FglevelNum = Input->FglevelNum;
	int dop1 = Input->dop1;
	double E = 100000;
	//	E = 1;
	//	lambda = E*nu / ((1 + nu)*(1 - 2 * nu));
	//	mu = E / 2.0 / (1 + nu);
	/*	printf("%lf, %lf, %lf, %lf\n", E, nu, lambda, mu);*/

	int dop2 = dop1 - 1;

	// press p
	getElementDOF(&elementDOF[0], elements->row, 1);
	// stress sigma
	getElementDOF_symTensor(&elementDOF[1], elements->row, 2);
	// velocity u
	getElementDOF_CR(&elementDOF[2], elements, elementEdge, edges->row, nodes->row);

	getFreenodesInfoCRVector2(edges, nodes, &elementDOF[2]);

	getTransposeOfelementDoF(&elementDOF[2], &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_stokesCR_psv_stressP2(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, mu, t);

	/** Step 3. Check matrix properties */
	int i;
	for (i = 0; i<3; i++)
	{
		check_symm(&A[i]);
		check_diagpos(&A[i]);
		check_diagdom(&A[i]);
	}

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
	printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
	if (A[3].row>0)
		printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);

	ASP_param aspparam; /* parameters for AMG */

	aspparam.print_level = Input->print_level;
	aspparam.max_iter = Input->itsolver_maxit;
	aspparam.tol = Input->itsolver_tol;
	aspparam.restart = Input->restart;

	aspparam.precond_type = Input->precond_type;
	aspparam.precond_scale = Input->precond_scale;
	aspparam.smoother = Input->smoother;
	aspparam.schwarz_type = Input->schwarz_type;
	aspparam.smooth_iter = Input->smooth_iter;

	aspparam.levelNum = Input->glevelNum;
	aspparam.mg_max_iter = Input->MG_maxit;
	aspparam.mg_tol = Input->MG_tol;
	aspparam.mg_smoother = Input->MG_smoother;
	aspparam.mg_smooth_iter = Input->MG_smooth_iter;

	aspparam.elements = elements;
	aspparam.elementEdge = elementEdge;
	aspparam.edges = edges;
	aspparam.nodes = nodes;
	aspparam.edgesTran = edgesTran;
	aspparam.nodeCEdge = nodeCEdge;
	aspparam.elementDOF = elementDOF;
	aspparam.elementdofTran = &elementdofTran;

	aspparam.lambda = Input->lambda;
	aspparam.mu = Input->mu;
	aspparam.t = Input->t;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;


	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}

	TriAsElasCR_GMRES(A, b, _uh, &aspparam, itsolver_type, print_level);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[i] = _uh[0].val[i];
	for (i = 0; i < _uh[1].row; i++)
		uh[1].val[elementDOF[2].freenodes.val[i]] = _uh[1].val[i];

	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);

	for (i = 0; i<3; i++)
		free_csr_matrix(&A[i]);
	if (A[3].row>0)
		free_csr_matrix(&A[3]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}