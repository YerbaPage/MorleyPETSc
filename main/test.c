/**
 *		Test function for solving a sparse SPD linear system using PCG and AMG. 
 *
 *		Matrix-Solver Community Project
 *
 *------------------------------------------------------
 *
 *		Created by Xuehai Huang on 07/11/2013.
 *		Copyright 2013 WZU. All rights reserved. 
 *
 *------------------------------------------------------
 *
 */

/*! \file test.c
 *  \brief Test Function for Solvers
 */
// static char help[] = "Solves a linear system with KSP.\n\n";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <petscksp.h>
#include <petscvec.h>

#include "header.h"
#include "matvec.h"
#include "checkmat.h"

/**
 * \fn int main (int argc, const char * argv[])
 *
 * This is the main function for test purpose. It contains five steps:
 */
int main(int argc, const char * argv[], char **args)
{
//	dCSRmat A[4];
//	dvector b[2], uh[2], _uh[2];
	dvector uh[2];
	ELEMENT *elements;
	idenmat *elementEdge;
	EDGE *edges;
	dennode *nodes;
	iCSRmat *edgesTran;
	ivector nodeCEdge;
	ELEMENT_DOF elementDOF;

	/** Step 0. Read input parameters */
	char *inputfile = "ini/input.dat";
	Input_data Input;
	read_Input_data(inputfile, &Input);

	int glevelNum = Input.glevelNum;
	int domain_num = Input.domain_num;

	elements = (ELEMENT*)calloc(glevelNum, sizeof(ELEMENT));
	elementEdge = (idenmat*)calloc(glevelNum, sizeof(idenmat));
	edges = (EDGE*)calloc(glevelNum, sizeof(EDGE));
	nodes = (dennode*)calloc(glevelNum, sizeof(dennode));
	edgesTran = (iCSRmat*)calloc(glevelNum, sizeof(iCSRmat));


	/** Step 1. generate mesh */
	if (getmesh(domain_num, elements, elementEdge, edges, nodes, edgesTran, &nodeCEdge, glevelNum) == 0)
	{
		printf("It's fail to generate mesh.\n");
		return 0;
	}

	int i;

	getElementEdgeGeoInfo(&elements[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1]);
//	getTensorBasis4HuZhang2d(&edges[glevelNum - 1], &edgesTran[glevelNum - 1], &nodes[glevelNum - 1]);

//	for (i = 0; i < glevelNum; i++)
//		getElementEdgeGeoInfo(&elements[i], &edges[i], &nodes[i]);

	printf("h=%f\n", edges[glevelNum - 1].length[1]);

//	Input.variationalform_type = 1;
	kirchhoffplatefem(uh, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &edgesTran[glevelNum - 1], &nodeCEdge, &elementDOF, &Input);


//	Input.variationalform_type = 2;
//	kirchhoffplatefem(uh+1, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &edgesTran[glevelNum - 1], &nodeCEdge, &elementDOF, &Input);

/*	dvector uh1[2], uh2[2];
	Input.variationalform_type = 2;
	Input.stress_fem_type = 1;
	stokesfem(uh1, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &edgesTran[glevelNum - 1], &nodeCEdge, elementDOF, &Input);
	Input.variationalform_type = 2;
	Input.stress_fem_type = 2;
	stokesfem(uh2, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &edgesTran[glevelNum - 1], &nodeCEdge, elementDOF, &Input);
	double error = 0;
	for (i = 0; i < uh1[1].row; i++)
		error += (uh1[1].val[i]- uh2[1].val[i])*(uh1[1].val[i] - uh2[1].val[i]);
	printf("error %le\n", sqrt(error));
	error = 0;
	for (i = 0; i < uh1[1].row; i++)
		error += (uh1[1].val[i])*(uh1[1].val[i]);
	printf("energy %le\n", sqrt(error));
	error = 0;
	for (i = 0; i < uh2[1].row; i++)
		error += (uh2[1].val[i])*(uh2[1].val[i]);
	printf("energy %le\n", sqrt(error));*/

	/** Step 5. Compute the error between numerical solution and exact solution */
	double errors[3];

	geterrors_kirchhoffpalte_morley(errors, uh, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &elementDOF);
	printf("Error between Morley solution and true solution:\n");
	printf("L2 norm of error = %e\n", errors[0]);
	printf("H1 norm of error = %e\n", errors[1]);
	printf("H2 norm of error = %e\n", errors[2]);

/*	geterrors_kirchhoffpalte_morley(errors, uh+1, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &elementDOF);
	printf("Error between PSP  solution and true solution:\n");
	printf("L2 norm of error = %e\n", errors[0]);
	printf("H1 norm of error = %e\n", errors[1]);
	printf("H2 norm of error = %e\n", errors[2]);

	geterrors2discr_kirchhoffpalte_morley(errors, uh, uh+1, &elements[glevelNum - 1], &elementEdge[glevelNum - 1], &edges[glevelNum - 1], &nodes[glevelNum - 1], &elementDOF);
	printf("Error between Morley solution and PSP solution:\n");
	printf("L2 norm of error = %e\n", errors[0]);
	printf("H1 norm of error = %e\n", errors[1]);
	printf("H2 norm of error = %e\n", errors[2]);*/

	/*********************************************************************************************/
	FILE *outputFile;
	outputFile = fopen("output/error.dat", "w");
	fprintf(outputFile, "%e\t%e\t%e\n", errors[0], errors[1], errors[2]);
	fclose(outputFile);
	/********************************************************************************************/

	getElementsNodes4Matlab(&elements[glevelNum - 1], &nodes[glevelNum - 1], &uh[0]);

	
	for (i = 0; i < glevelNum; i++)
	{
		free_ELEMENT(&elements[i]);
		free_iden_matrix(&elementEdge[i]);
		free_EDGE(&edges[i]);
		free_dennode(&nodes[i]);
		free_icsr_matrix(&edgesTran[i]);
	}
	free_ivector(&nodeCEdge);
	free_dvector(&uh[0]);
//	free_dvector(&uh[1]);
	free(elements);
	free(elementEdge);
	free(edges);
	free(nodes);
	free(edgesTran);
	free_elementDOF(&elementDOF);

	return 1;

    // Vec x, b, u;    /* approx solution, RHS, exact solution */
    // Mat A;          /* linear system matrix */
    // KSP ksp;        /* linear solver context */
    // PC pc;          /* preconditioner context */
    // PetscReal norm; /* norm of solution error */
    // PetscErrorCode ierr;
    // PetscInt k, n = 100, col[3], its;
    // PetscMPIInt size;
    // PetscScalar value[3];
 
}
