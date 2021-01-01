/*
 *  fem.c
 *
 *  Created by Xuehai Huang on 05/21/17.
 *  Copyright 2017 WZU. All rights reserved.
 *
 */

/*! \file kirchhoffplatefem.c
 *  \brief Assembling for stiffness matrix and solve it
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"
#include "checkmat.h"

 /**
 * \fn void kirchhoffplatefem(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
 * \brief finite element methods for Kirchhoff plate problem
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
void kirchhoffplatefem(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	int i;
	int variationalform_type = Input->variationalform_type;
	if (variationalform_type == 1) // primal formulation
	{
		printf("Discretize Kirchhoff plate problem in original form.\n");
		kirchhoffplateMorley(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else if (variationalform_type == 2) // Poisson-StokesDiv-Poisson formulation
	{
		printf("Discretize Kirchhoff plate problem in Poisson-StokesDiv-Poisson form.\n");
		kirchhoffplatefem_psp(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else if (variationalform_type == 3) // Poisson-StokesRot-Poisson formulation
	{
		printf("Discretize Kirchhoff plate problem in Poisson-StokesRot-Poisson form.\n");
		kirchhoffplatefem_psprot(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else if (variationalform_type == 4) // Poisson-Elasticity-Poisson formulation
	{
		printf("Discretize Kirchhoff plate problem in Poisson-Elasticity-Poisson form.\n");
		kirchhoffplatefem_pep(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	}
	else
	{
		printf("Please set variationalform_type = 1, 2, or 3!\n");
		exit(0);
	}
}

/**
* \fn void kirchhoffplateMorley(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief Morley element method for Kirchhoff plate problem in original form
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
void kirchhoffplateMorley(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	dCSRmat A;
	dvector b, _uh;
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

	getElementDOF_Morley(elementDOF, elements, elementEdge, edges->row, nodes->row);

	getFreenodesInfoMorley(edges, nodes, elementDOF);

	getTransposeOfelementDoF(elementDOF, &elementdofTran, 0);

	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_kirchhoffplateMorley(&A, &b, uh, elements, elementEdge, edges, nodes, elementDOF, &elementdofTran, nu);

	/** Step 3. Check matrix properties */
	int i;
/*	check_symm(&A);
	check_diagpos(&A);
	check_diagdom(&A);*/

	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b.row, &_uh);

	printf("A.row=%d, A.col=%d, A.nnz=%d\n", A.row, A.col, A.nnz);

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
	aspparam.nu = Input->nu;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input->print_level;
	amgparam.max_levels = Input->AMG_levels;
	amgparam.max_iter = Input->itsolver_maxit;
	amgparam.tol = Input->itsolver_tol;
	amgparam.AMG_max_iter = Input->MG_maxit;
	amgparam.AMG_tol = Input->MG_tol;

	amgparam.smoother = Input->MG_smoother;
	amgparam.presmooth_iter = Input->MG_smooth_iter;
	amgparam.postsmooth_iter = Input->MG_smooth_iter;

	amgparam.coarsening_type = Input->AMG_coarsening_type;
	amgparam.interpolation_type = Input->AMG_interpolation_type;
	amgparam.coarse_dof = Input->AMG_coarse_dof;
	amgparam.strong_threshold = Input->AMG_strong_threshold;
	amgparam.truncation_threshold = Input->AMG_truncation_threshold;
	amgparam.max_row_sum = Input->AMG_max_row_sum;

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

	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A, &b, &_uh, &amgparam, print_level);

	for (i = 0; i < _uh.row; i++)
		uh->val[elementDOF[0].freenodes.val[i]] = _uh.val[i];

	free_dvector(&_uh);

	free_csr_matrix(&A);
	free_dvector(&b);

	free(elementdofTran.IA);
	free(elementdofTran.JA);

}

/**
* \fn void kirchhoffplatefem_psp(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief finite element methods for Kirchhoff plate problem in Poisson-Stokes-Poisson form
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
void kirchhoffplatefem_psp(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
//	stokesNcP1P0_elas(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	kirchhoffplateMorley_psp(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
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
* \fn void kirchhoffplatefem_psprot(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief finite element methods for Kirchhoff plate problem in Poisson-Stokes-Poisson form
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
void kirchhoffplatefem_psprot(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	kirchhoffplateMorley_psprot(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
}

/**
* \fn void kirchhoffplatefem_pep(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
* \brief finite element methods for Kirchhoff plate problem in Poisson-Elasticity-Poisson form
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
void kirchhoffplatefem_pep(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input)
{
	kirchhoffplateMorley_pep(uh, elements, elementEdge, edges, nodes, edgesTran, nodeCEdge, elementDOF, Input);
	/*int stress_fem_type = Input->stress_fem_type;
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
	}*/
}

/**
* \fn void kirchhoffplateMorley_psp(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input)
* \brief Morley element method for Kirchhoff plate problem in Poisson-Stokes-Poisson form
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
void kirchhoffplateMorley_psp(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input)
{
	dCSRmat As[4], A1pm;
	dvector b[2], _uh[2], A2pm, uhpm[2], uhs[2];
	ELEMENT_DOF elementDOFpm[2], elementDOFs[2];
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
	int i;

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

	aspparam.lambda = Input->lambda;
	aspparam.mu = 0.5;// Input->mu;
	aspparam.t = Input->t;
	aspparam.nu = Input->nu;
	aspparam.problem_num = Input->problem_num;
	if (aspparam.problem_num == 1)
		aspparam.mu = 0.5;
	else
		aspparam.mu = (1 - aspparam.nu) / 2;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input->print_level;
	amgparam.max_levels = Input->AMG_levels;
	amgparam.max_iter = Input->itsolver_maxit;
	amgparam.tol = Input->itsolver_tol;
	amgparam.AMG_max_iter = Input->MG_maxit;
	amgparam.AMG_tol = Input->MG_tol;

	amgparam.smoother = Input->MG_smoother;
	amgparam.presmooth_iter = Input->MG_smooth_iter;
	amgparam.postsmooth_iter = Input->MG_smooth_iter;

	amgparam.coarsening_type = Input->AMG_coarsening_type;
	amgparam.interpolation_type = Input->AMG_interpolation_type;
	amgparam.coarse_dof = Input->AMG_coarse_dof;
	amgparam.strong_threshold = Input->AMG_strong_threshold;
	amgparam.truncation_threshold = Input->AMG_truncation_threshold;
	amgparam.max_row_sum = Input->AMG_max_row_sum;

	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	/************ Eqn1: Poisson equation disctetized by Morley element ************/
	/** Step 1. generate degrees of freedom */
	getElementDOF_Lagrange(&elementDOFpm[0], elements, elementEdge, edges, nodes->row, 1);
	getElementDOF_NoncfmP1(&elementDOFpm[1], elementEdge, edges->row);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFpm[0]);
	getFreenodesInfoNcP1(edges, &elementDOFpm[1]);
	getTransposeOfelementDoF(&elementDOFpm[0], &elementdofTran, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_poissonMorley(&A1pm, &A2pm, b, uhpm, elements, elementEdge, edges, nodes, elementDOFpm, &elementdofTran);
	free_icsr_matrix(&elementdofTran);
	/** Step 3. Check matrix properties */
	printf("A1pm.row=%d, A1pm.col=%d, A1pm.nnz=%d\n", A1pm.row, A1pm.col, A1pm.nnz);
	printf("A2pm.row=%d, A2pm.col=%d, A2pm.nnz=%d\n", A2pm.row, A2pm.row, A2pm.row);
	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A1pm, &b[0], &_uh[0], &amgparam, print_level);
//	classicAMG(&A1pm, &b[0], &_uh[0], &amgparam);
	for (i = 0; i < _uh[0].row; i++)
		uhpm[0].val[elementDOFpm[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < A2pm.row; i++)
		uhpm[1].val[elementDOFpm[1].freenodes.val[i]] = b[1].val[i] / A2pm.val[i];
	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn2: Stokes equation disctetized by nonconforming P1-P0 element ************/
	/** Step 1. generate degrees of freedom */
	// velocity u
	getElementDOF_NoncfmP1(&elementDOFs[0], elementEdge, edges->row);
	// press p
	getElementDOF(&elementDOFs[1], elements->row, 0);
	getFreenodesInfoNcP1Vector2(edges, &elementDOFs[0]);
	getTransposeOfelementDoF(&elementDOFs[0], &elementdofTran, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_kirchhoffplate_StokesNcP1P0(As, b, uhs, elements, elementEdge, edges, nodes, elementDOFs, &elementdofTran, uhpm, elementDOFpm, problem_num, nu);
	free_dvector(&uhpm[0]);
	free_dvector(&uhpm[1]);
	free_icsr_matrix(&elementdofTran);
	/** Step 3. Check matrix properties */
	printf("As.row=%d, As.col=%d, As.nnz=%d\n", As[0].row, As[0].col, As[0].nnz);
	printf("Bs.row=%d, Bs.col=%d, Bs.nnz=%d\n", As[2].row, As[2].col, As[2].nnz);
	/** Step 4. Solve the system */
	aspparam.elementDOF = elementDOFs;
//	aspparam.elementdofTran = &elementdofTran;
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);
	init_dvector(&_uh[0], 1);///////////////////////////////////

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	if (itsolver_type == 1)
	{
		printf("\nASP Approximate Block Factorization preconditioned GMRES solver with auxiliary space method.\n\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		AbfpAsP1StokesNcP1P0_GMRES(As, b, _uh, &aspparam, print_level);
	//	Abfp2AsP1StokesNcP1P0_GMRES(As, b, _uh, &aspparam, print_level);
	}
	else
	{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1StokesNcP1P0_MINRES(As, b, _uh, &aspparam, print_level);
	//	DiagAMGStokesNcP1P0_MINRES(As, b, _uh, &aspparam, print_level); // sometimes better
	}
	for (i = 0; i < _uh[0].row; i++)
		uhs[0].val[elementDOFs[0].freenodes.val[i]] = _uh[0].val[i];
	/*	for (i = 0; i < _uh[1].row; i++)
	{
	if (problem_num == 1)
	uhs[1].val[i] = _uh[1].val[i];
	else
	uhs[1].val[i] = _uh[1].val[i] * (1 - nu);
	}*/
//	for (i = 0; i < _uh[1].row; i++)
//		uhs[1].val[i] = _uh[1].val[i];
	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);
	for (i = 0; i<3; i++)
		free_csr_matrix(&As[i]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn3: Poisson equation disctetized by Morley element ************/
	/** Step 1. assemble the right hand side term */
	assembleRHScurl_biharmonic_PoissonMorley2d(b, &uhs[0], elements, elementEdge, edges, elementDOFpm, &elementDOFs[0]);
	free_dvector(&uhs[0]);
	free_dvector(&uhs[1]);
	free_elementDOF(&elementDOFs[0]);
	free_elementDOF(&elementDOFs[1]);
	/** Step 2. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A1pm, &b[0], &_uh[0], &amgparam, print_level);

	getElementDOF_Morley(elementDOFm, elements, elementEdge, edges->row, nodes->row);
	getFreenodesInfoMorley(edges, nodes, elementDOFm);
	create_dvector(elementDOFm->dof, uh);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[elementDOFpm[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < A2pm.row; i++)
		uh[0].val[elementDOFpm[0].dof + elementDOFpm[1].freenodes.val[i]] = b[1].val[i] / A2pm.val[i];
//	printf("%d, %d, %d\n", elementDOFpm[0].dof, elementDOFpm[1].dof, uh[0].row);
	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);
	free_csr_matrix(&A1pm);
	free_dvector(&A2pm);
	free_elementDOF(&elementDOFpm[0]);
	free_elementDOF(&elementDOFpm[1]);
}

/**
* \fn void kirchhoffplateMorley_psprot(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input)
* \brief Morley element method for Kirchhoff plate problem in Poisson-Stokes-Poisson form
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
void kirchhoffplateMorley_psprot(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input)
{
	dCSRmat As[4], A1pm;
	dvector b[2], _uh[2], A2pm, uhpm[2], uhs[2];
	ELEMENT_DOF elementDOFpm[2], elementDOFs[2];
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
	int i;

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

	aspparam.lambda = Input->lambda;
	aspparam.mu = 0.5;// Input->mu;
	aspparam.t = Input->t;
	aspparam.nu = Input->nu;
	aspparam.problem_num = Input->problem_num;
	if (aspparam.problem_num == 1)
		aspparam.mu = 0.5;
	else
		aspparam.mu = (1 - aspparam.nu) / 2;

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input->print_level;
	amgparam.max_levels = Input->AMG_levels;
	amgparam.max_iter = Input->itsolver_maxit;
	amgparam.tol = Input->itsolver_tol;
	amgparam.AMG_max_iter = Input->MG_maxit;
	amgparam.AMG_tol = Input->MG_tol;

	amgparam.smoother = Input->MG_smoother;
	amgparam.presmooth_iter = Input->MG_smooth_iter;
	amgparam.postsmooth_iter = Input->MG_smooth_iter;

	amgparam.coarsening_type = Input->AMG_coarsening_type;
	amgparam.interpolation_type = Input->AMG_interpolation_type;
	amgparam.coarse_dof = Input->AMG_coarse_dof;
	amgparam.strong_threshold = Input->AMG_strong_threshold;
	amgparam.truncation_threshold = Input->AMG_truncation_threshold;
	amgparam.max_row_sum = Input->AMG_max_row_sum;

	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	/************ Eqn1: Poisson equation disctetized by Morley element ************/
	/** Step 1. generate degrees of freedom */
	getElementDOF_Lagrange(&elementDOFpm[0], elements, elementEdge, edges, nodes->row, 1);
	getElementDOF_NoncfmP1(&elementDOFpm[1], elementEdge, edges->row);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFpm[0]);
	getFreenodesInfoNcP1(edges, &elementDOFpm[1]);
	getTransposeOfelementDoF(&elementDOFpm[0], &elementdofTran, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_poissonMorley(&A1pm, &A2pm, b, uhpm, elements, elementEdge, edges, nodes, elementDOFpm, &elementdofTran);
	free_icsr_matrix(&elementdofTran);
	/** Step 3. Check matrix properties */
	printf("A1pm.row=%d, A1pm.col=%d, A1pm.nnz=%d\n", A1pm.row, A1pm.col, A1pm.nnz);
	printf("A2pm.row=%d, A2pm.col=%d, A2pm.nnz=%d\n", A2pm.row, A2pm.row, A2pm.row);
	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A1pm, &b[0], &_uh[0], &amgparam, print_level);
	//	classicAMG(&A1pm, &b[0], &_uh[0], &amgparam);
	for (i = 0; i < _uh[0].row; i++)
		uhpm[0].val[elementDOFpm[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < A2pm.row; i++)
		uhpm[1].val[elementDOFpm[1].freenodes.val[i]] = b[1].val[i] / A2pm.val[i];
	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn2: Stokes equation disctetized by nonconforming P1-P0 element ************/
	/** Step 1. generate degrees of freedom */
	// velocity u
	getElementDOF_NoncfmP1(&elementDOFs[0], elementEdge, edges->row);
	// press p
	getElementDOF(&elementDOFs[1], elements->row, 0);
	getFreenodesInfoNcP1Vector2(edges, &elementDOFs[0]);
	getTransposeOfelementDoF(&elementDOFs[0], &elementdofTran, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_kirchhoffplate_StokesrotNcP1P0(As, b, uhs, elements, elementEdge, edges, nodes, elementDOFs, &elementdofTran, uhpm, elementDOFpm, problem_num, nu);
	free_dvector(&uhpm[0]);
	free_dvector(&uhpm[1]);
	free_icsr_matrix(&elementdofTran);
	/** Step 3. Check matrix properties */
	printf("As.row=%d, As.col=%d, As.nnz=%d\n", As[0].row, As[0].col, As[0].nnz);
	printf("Bs.row=%d, Bs.col=%d, Bs.nnz=%d\n", As[2].row, As[2].col, As[2].nnz);
	/** Step 4. Solve the system */
	aspparam.elementDOF = elementDOFs;
	//	aspparam.elementdofTran = &elementdofTran;
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);
	init_dvector(&_uh[0], 1);///////////////////////////////////

	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	if (itsolver_type == 1)
	{
		printf("\nASP Approximate Block Factorization preconditioned GMRES solver with auxiliary space method.\n\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		AbfpAsP1StokesNcP1P0_GMRES(As, b, _uh, &aspparam, print_level);
	//	Abfp2AsP1StokesNcP1P0_GMRES(As, b, _uh, &aspparam, print_level);
	}
	else
	{
		printf("\nASP Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		printf("Auxiliary space: P1 Lagrangian element.\n\n");
		DiagAsP1StokesNcP1P0_MINRES(As, b, _uh, &aspparam, print_level);
		//	DiagAMGStokesNcP1P0_MINRES(As, b, _uh, &aspparam, print_level); // sometimes better
	}
	for (i = 0; i < _uh[0].row; i++)
		uhs[0].val[elementDOFs[0].freenodes.val[i]] = _uh[0].val[i];
	/*	for (i = 0; i < _uh[1].row; i++)
	{
	if (problem_num == 1)
	uhs[1].val[i] = _uh[1].val[i];
	else
	uhs[1].val[i] = _uh[1].val[i] * (1 - nu);
	}*/
	//	for (i = 0; i < _uh[1].row; i++)
	//		uhs[1].val[i] = _uh[1].val[i];
	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);
	for (i = 0; i<3; i++)
		free_csr_matrix(&As[i]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn3: Poisson equation disctetized by Morley element ************/
	/** Step 1. assemble the right hand side term */
	assembleRHSgrad_biharmonic_PoissonMorley2d(b, &uhs[0], elements, elementEdge, edges, elementDOFpm, &elementDOFs[0]);
	free_dvector(&uhs[0]);
	free_dvector(&uhs[1]);
	free_elementDOF(&elementDOFs[0]);
	free_elementDOF(&elementDOFs[1]);
	/** Step 2. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A1pm, &b[0], &_uh[0], &amgparam, print_level);

	getElementDOF_Morley(elementDOFm, elements, elementEdge, edges->row, nodes->row);
	getFreenodesInfoMorley(edges, nodes, elementDOFm);
	create_dvector(elementDOFm->dof, uh);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[elementDOFpm[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < A2pm.row; i++)
		uh[0].val[elementDOFpm[0].dof + elementDOFpm[1].freenodes.val[i]] = b[1].val[i] / A2pm.val[i];
	//	printf("%d, %d, %d\n", elementDOFpm[0].dof, elementDOFpm[1].dof, uh[0].row);
	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);
	free_csr_matrix(&A1pm);
	free_dvector(&A2pm);
	free_elementDOF(&elementDOFpm[0]);
	free_elementDOF(&elementDOFpm[1]);
}

/**
* \fn void kirchhoffplateMorley_pep(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input)
* \brief Morley element method for Kirchhoff plate problem in Poisson-Elasticity-Poisson form
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
void kirchhoffplateMorley_pep(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input)
{
	dCSRmat Ae[4], A1pm;
	dvector b[2], _uh[2], A2pm, uhpm[2], uhe[2];
	ELEMENT_DOF elementDOFpm[2], elementDOFe[2];
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
	int i;

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

	aspparam.lambda = Input->lambda;
	aspparam.mu = 0.5; // Input->mu;
	aspparam.t = Input->t;
	aspparam.nu = Input->nu;
	aspparam.problem_num = Input->problem_num;
/*	if (aspparam.problem_num == 1)
		aspparam.mu = 0.5;
	else
		aspparam.mu = (1 - aspparam.nu) / 2;*/

	aspparam.variationalform_type = Input->variationalform_type;
	aspparam.stress_fem_type = Input->stress_fem_type;

	aspparam.max_levels = Input->AMG_levels;
	aspparam.AMG_coarsening_type = Input->AMG_coarsening_type;
	aspparam.AMG_interpolation_type = Input->AMG_interpolation_type;
	aspparam.AMG_coarse_dof = Input->AMG_coarse_dof;
	aspparam.AMG_strong_threshold = Input->AMG_strong_threshold;
	aspparam.AMG_truncation_threshold = Input->AMG_truncation_threshold;
	aspparam.AMG_max_row_sum = Input->AMG_max_row_sum;

	AMG_param amgparam; /* parameters for AMG */

	amgparam.print_level = Input->print_level;
	amgparam.max_levels = Input->AMG_levels;
	amgparam.max_iter = Input->itsolver_maxit;
	amgparam.tol = Input->itsolver_tol;
	amgparam.AMG_max_iter = Input->MG_maxit;
	amgparam.AMG_tol = Input->MG_tol;

	amgparam.smoother = Input->MG_smoother;
	amgparam.presmooth_iter = Input->MG_smooth_iter;
	amgparam.postsmooth_iter = Input->MG_smooth_iter;

	amgparam.coarsening_type = Input->AMG_coarsening_type;
	amgparam.interpolation_type = Input->AMG_interpolation_type;
	amgparam.coarse_dof = Input->AMG_coarse_dof;
	amgparam.strong_threshold = Input->AMG_strong_threshold;
	amgparam.truncation_threshold = Input->AMG_truncation_threshold;
	amgparam.max_row_sum = Input->AMG_max_row_sum;

	AFEM_param afemparam; /* parameters for AFEM */
	afemparam.aspparam = &aspparam;
	afemparam.tol = Input->AFEM_tol;
	afemparam.max_iter = Input->AFEM_maxit;
	afemparam.solver_type = Input->itsolver_type;
	afemparam.mark_threshold = Input->AFEM_mark_threshold;

	/************ Eqn1: Poisson equation disctetized by Morley element ************/
	/** Step 1. generate degrees of freedom */
	getElementDOF_Lagrange(&elementDOFpm[0], elements, elementEdge, edges, nodes->row, 1);
	getElementDOF_NoncfmP1(&elementDOFpm[1], elementEdge, edges->row);
	getFreenodesInfoLagrange(edges, nodes, &elementDOFpm[0]);
	getFreenodesInfoNcP1(edges, &elementDOFpm[1]);
	getTransposeOfelementDoF(&elementDOFpm[0], &elementdofTran, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_poissonMorley(&A1pm, &A2pm, b, uhpm, elements, elementEdge, edges, nodes, elementDOFpm, &elementdofTran);
	free_icsr_matrix(&elementdofTran);
	/** Step 3. Check matrix properties */
	printf("A1pm.row=%d, A1pm.col=%d, A1pm.nnz=%d\n", A1pm.row, A1pm.col, A1pm.nnz);
	printf("A2pm.row=%d, A2pm.col=%d, A2pm.nnz=%d\n", A2pm.row, A2pm.row, A2pm.row);
	/** Step 4. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A1pm, &b[0], &_uh[0], &amgparam, print_level);
//	classicAMG(&A1pm, &b[0], &_uh[0], &amgparam);
	for (i = 0; i < _uh[0].row; i++)
		uhpm[0].val[elementDOFpm[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < A2pm.row; i++)
		uhpm[1].val[elementDOFpm[1].freenodes.val[i]] = b[1].val[i] / A2pm.val[i];
	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn2: Elasticity problem disctetized by nonconforming P1-P0 element ************/
	/** Step 1. generate degrees of freedom */
	// stress sigma
	getElementDOF(&elementDOFe[0], elements->row, 0);
	// velocity u
	getElementDOF_NoncfmP1(&elementDOFe[1], elementEdge, edges->row);
	getFreenodesInfoNcP1Vector2(edges, &elementDOFe[1]);
	getTransposeOfelementDoF(&elementDOFe[1], &elementdofTran, 0);
	/** Step 2. assemble stiffmatrix and right hand side term */
	assemble_kirchhoffplate_ElasNcP1P0(Ae, b, uhe, elements, elementEdge, edges, nodes, elementDOFe, &elementdofTran, uhpm, elementDOFpm, problem_num, nu);
	free_dvector(&uhpm[0]);
	free_dvector(&uhpm[1]);
	free_icsr_matrix(&elementdofTran);
	/** Step 3. Check matrix properties */
	printf("Ae.row=%d, Ae.col=%d, Ae.nnz=%d\n", Ae[0].row, Ae[0].col, Ae[0].nnz);
	printf("Be.row=%d, Be.col=%d, Be.nnz=%d\n", Ae[2].row, Ae[2].col, Ae[2].nnz);
	/** Step 4. Solve the system */
	aspparam.elementDOF = elementDOFe;
	//	aspparam.elementdofTran = &elementdofTran;
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	create_dvector(b[1].row, &_uh[1]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	if (itsolver_type == 1)
	{
		printf("Block triangular preconditioned GMRES solver with auxiliary space method\n");
		TriAsP1ElasNcP1P0_GMRES(Ae, b, _uh, &aspparam, print_level);
	}
	else
	{
		printf("Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		DiagAsP1ElasNcP1P0_MINRES(Ae, b, _uh, &aspparam, print_level);
	}
/*	for (i = 0; i < _uh[0].row; i++)
	{
		if (problem_num == 1)
			uhe[0].val[i] = _uh[0].val[i];
		else
			uhe[0].val[i] = _uh[0].val[i] * (1 - nu);
	}*/
	for (i = 0; i < _uh[1].row; i++)
		uhe[1].val[elementDOFe[1].freenodes.val[i]] = _uh[1].val[i];
	free_dvector(&_uh[0]);
	free_dvector(&_uh[1]);
	for (i = 0; i<3; i++)
		free_csr_matrix(&Ae[i]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);

	/************ Eqn3: Poisson equation disctetized by Morley element ************/
	/** Step 1. assemble the right hand side term */
	assembleRHScurl_biharmonic_PoissonMorley2d(b, &uhe[1], elements, elementEdge, edges, elementDOFpm, &elementDOFe[1]);
	free_dvector(&uhe[0]);
	free_dvector(&uhe[1]);
	free_elementDOF(&elementDOFe[0]);
	free_elementDOF(&elementDOFe[1]);
	/** Step 2. Solve the system */
	printf("Solve...\n");
	create_dvector(b[0].row, &_uh[0]);
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	printf("Solving Morley by AMG preconditioned CG solver\n");
	classicAMG_PCG(&A1pm, &b[0], &_uh[0], &amgparam, print_level);

	getElementDOF_Morley(elementDOFm, elements, elementEdge, edges->row, nodes->row);
	getFreenodesInfoMorley(edges, nodes, elementDOFm);
	create_dvector(elementDOFm->dof, uh);

	for (i = 0; i < _uh[0].row; i++)
		uh[0].val[elementDOFpm[0].freenodes.val[i]] = _uh[0].val[i];
	for (i = 0; i < A2pm.row; i++)
		uh[0].val[elementDOFpm[0].dof + elementDOFpm[1].freenodes.val[i]] = b[1].val[i] / A2pm.val[i];
	//	printf("%d, %d, %d\n", elementDOFpm[0].dof, elementDOFpm[1].dof, uh[0].row);
	free_dvector(&_uh[0]);
	free_dvector(&b[0]);
	free_dvector(&b[1]);
	free_csr_matrix(&A1pm);
	free_dvector(&A2pm);
	free_elementDOF(&elementDOFpm[0]);
	free_elementDOF(&elementDOFpm[1]);
}
