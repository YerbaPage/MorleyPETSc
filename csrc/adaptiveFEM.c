/*
 *  adaptiveFEM.c
 *  
 *
 *  Created by Xuehai Huang on 08/24/2010.
 *  Copyright 2010 WZU. All rights reserved.
 *
 */

/*! \file adaptiveFEM.c
 *  \brief adaptive finite element method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "precond.h"
#include "matvec.h"

/**
 * \fn double adaptiveFEM(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, iCSRmat *edgesTran, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uh, dvector *uhstar, double lambda, double mu, AFEM_param *afemparam)
 * \brief compute the error estimator for adaptive finite element method
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *edgesTran the relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *uh pointer to numerical solution
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \param *estimator pointer to the error estimator
 * \return total posterior error
 */
double adaptiveFEM(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, iCSRmat *edgesTran, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uh, dvector *uhstar, double lambda, double mu, AFEM_param *afemparam)
{
	int i;
	dCSRmat A[4];
	dvector b[2], _uh[2];
	double totalError, markError;
	dvector estimator;
	ivector estimatorIndex;
	ivector marker;
	int iterNum = 0;

	int dop1 = elementDOF[0].dop;
	int dop2 = elementDOF[1].dop;

	ASP_param *aspparam = afemparam->aspparam;
	double afemTol = afemparam->tol;
	int afemMaxiter = afemparam->max_iter;
	int itsolver_type = afemparam->solver_type;
	double mark_threshold = afemparam->mark_threshold;

	/*********************************************************************************************/
	FILE *ofile1, *ofile2;
	ofile1 = fopen("output/priorierror.dat", "w");
	ofile2 = fopen("output/posteriorierror.dat", "w");
	double errors[10], posterrors[8];
	dvector Qhu;
	int tdofs;
	/*********************************************************************************************/


	while (1)
	{
		/* Assemble       */
		printf("Assemble...\n");
		assemble(A, b, uh, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, lambda, mu);

		printf("A.row=%d, A.col=%d, A.nnz=%d\n", A[0].row, A[0].col, A[0].nnz);
		printf("B.row=%d, B.col=%d, B.nnz=%d\n", A[1].row, A[1].col, A[1].nnz);
		if (A[3].row > 0)
			printf("C.row=%d, C.col=%d, C.nnz=%d\n", A[3].row, A[3].col, A[3].nnz);
		tdofs = A[1].row + A[1].col;
		//////////////////////////////////////////////
		/*printf("A:%d, %d\n", A[0].row, A[0].col);
		for (i = 0; i < A->nnz; i++)
		{
			if (fabs(A->val[i]) > 0.25)
				printf("(%d,%e), ", A->JA[i], A->val[i]);
		}
		printf("A end\n");
		printf("BT:%d, %d\n", A[1].row, A[1].col);
		for (i = 0; i < A[1].nnz; i++)
		{
			if (fabs(A[1].val[i])>0.25)
			printf("(%d,%e), ", A[1].JA[i], A[1].val[i]);
		}
		printf("BT end\n");
		printf("B:%d, %d, ", A[2].row, A[2].col);
		for (i = 0; i < A[2].nnz; i++)
		{
			if (fabs(A[2].val[i])>0.25)
				printf("(%d,%e),", A[2].JA[i], A[2].val[i]);
		}
		printf("B end\n");
		printf("b[0]:%d\n", b[0].row);
		for (i = 0; i < b[0].row; i++)
			printf("(%e),", b[0].val[i]);
		printf("b[0] end\n");
		printf("b[1]:%d\n", b[1].row);
		for (i = 0; i < b[1].row; i++)
			printf("(%e),", b[1].val[i]);
		printf("b[1] end\n");
				exit(1);*/
		/////////////////////////////////////////////

		/* Solve       */
		printf("Solve...\n");
		create_dvector(b[0].row, &_uh[0]);
		create_dvector(b[1].row, &_uh[1]);
		solve(A, b, _uh, aspparam, itsolver_type);

		for (i = 0; i < _uh[0].row; i++)
			uh[0].val[elementDOF[0].freenodes.val[i]] = _uh[0].val[i];
		for (i = 0; i < _uh[1].row; i++)
			uh[1].val[i] = _uh[1].val[i];

		free_dvector(&_uh[0]);
		free_dvector(&_uh[1]);

		/***************************write matrix to file**********************************************
		write_IJ_matrix(&A[0], "output/A11.dat");
		write_IJ_matrix(&A[1], "output/A12.dat");
		*********************************************************************************************/
		for (i = 0; i<3; i++)
				free_csr_matrix(&A[i]);
		if (A[3].row>0)
			free_csr_matrix(&A[3]);
		free_dvector(&b[0]);
		free_dvector(&b[1]);

		free(elementdofTran->IA);
		free(elementdofTran->JA);

		/* Postprocessing */
		printf("Postprocessing...\n");
		create_dvector(elementDOF[2].dof * 2, uhstar);////////////////////////////////////
	//	postprocess2newDisplacement(uhstar, &uh[0], &uh[1], elements, nodes, elementDOF, lambda, mu);
	
		printf("\nIterative number=%d\n", iterNum);
		iterNum++;

		/* estimate      */
		printf("Estimate...\n");
		totalError = estimate(elements, elementEdge, edges, nodes, elementDOF, &uh[0], &uh[1], uhstar, lambda, mu, &estimator);

		/*********************************output errors to files***********************************************/
		create_dvector(uh[1].row, &Qhu);
		projPiecewiseLagrangeDisplacement(&Qhu, elements, nodes, &elementDOF[1], lambda, mu);
		geterrors(errors, &uh[0], &uh[1], &Qhu, uhstar, elements, elementEdge, edges, nodes, elementDOF, lambda, mu);
		free_dvector(&Qhu);
		getposteriorierrors(posterrors, &uh[0], &uh[1], uhstar, elements, elementEdge, edges, nodes, elementDOF, lambda, mu);
		fprintf(ofile1, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", tdofs, errors[9], errors[0], errors[1], errors[2], errors[3], errors[4], errors[5], errors[6], errors[7], errors[8]);
		fprintf(ofile2, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", tdofs, posterrors[0], posterrors[1], posterrors[2], posterrors[5], posterrors[3], posterrors[4], posterrors[6], posterrors[7]);
		/*****************************************************************************************************/


		if(iterNum>afemMaxiter || totalError<afemTol)
		{
			fclose(ofile1);
			fclose(ofile2);
			free_dvector(&estimator);
			break;
		}

		free_elementDOF(elementDOF);
		free_elementDOF(elementDOF + 1);
		free_elementDOF(elementDOF + 2);
		free_dvector(&uh[0]);
		free_dvector(&uh[1]);
		free_dvector(uhstar);
			
		/* mark          */
		printf("mark...\n");
		create_ivector(elements->row, &estimatorIndex);
		for(i=0;i<estimatorIndex.row;i++)
			estimatorIndex.val[i]=i;

		quicksortIndex(estimator.val,estimatorIndex.val,0,estimator.row-1); // quicksort

		create_ivector(edges->row, &marker);
		markError = mark(elementEdge, edges, nodes, &estimator, estimatorIndex.val, mark_threshold, totalError, marker.val);

		printf("theta=%e, totalError=%e, markError=%e\n", mark_threshold, totalError, markError);///////////////////////
		
		free_EDGE(edges);
		free_icsr_matrix(edgesTran);
		free_dvector(&estimator);
		free_ivector(&estimatorIndex); 
		
		/* refine        */
		printf("refine...\n");
		refine(elements, elementEdge, edges, edgesTran, nodes, marker.val);

		free_ivector(&marker);

		getElementEdgeGeoInfo(elements, edges, nodes);
		getTensorBasis4HuZhang2d(edges, edgesTran, nodes);

		getElementDOF_HuZhang(&elementDOF[0], elements, elementEdge, edges, nodes->row, dop1);
		getElementDOF(&elementDOF[1], elements->row, dop2);
		getElementDOF(&elementDOF[2], elements->row, dop2 + 2);

		getFreenodesInfoHuZhang(edges, nodes, &elementDOF[0]);
		
		getTransposeOfelementDoF(elementDOF, elementdofTran, 0);
	
		aspparam->elements = elements;
		aspparam->elementEdge = elementEdge;
		aspparam->edges = edges;
		aspparam->nodes = nodes;
//		aspparam->edgesTran = &edgesTran[glevelNum - 1];
//		aspparam->nodeCEdge = &nodeCEdge;
		aspparam->elementDOF = elementDOF;
		aspparam->elementdofTran = elementdofTran;

		
	}

	return sqrt(totalError);
}


/**
 * \fn int solve(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type)
 * \brief Solve Ax=b
 * \param *A	pointer to the dCSRmat matrix
 * \param *b	pointer to the dvector of right hand side
 * \param *x	pointer to the dvector of dofs
 * \param *aspparam pointer to ASP parameters
 * \param itsolver_type type of iterative solver
 * \return the number of iterations
 */
int solve(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type)
{
	int iter=0;
	int print_level = aspparam->print_level;
	int itsolver_maxit = aspparam->max_iter;
	double itsolver_tol = aspparam->tol;
	
	if (print_level>0) {
		printf("Maximal iteration number = %d\n", itsolver_maxit);
		printf("Tolerance for rel. res.  = %e\n", itsolver_tol);
	}
	
	if (itsolver_type == 1)
	{
		printf("Block diagonal preconditioned MINRES solver with auxiliary space method\n");
		iter = DiagAsP1ElasDG_MINRES(A, b, x, aspparam, print_level);
	}
	else
	{
		printf("Block triangular preconditioned GMRES solver with auxiliary space method\n");
		iter = TriAsP1ElasDG_GMRES(A, b, x, aspparam, print_level);
	}
	
	return iter;
}

//ESTIMATE
// see file estimate.c

//MARK

//REFINE