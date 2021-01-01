/*
 *  multigrid.c
 *  classicAMG
 *
 *  Created by Xuehai Huang on 2/2/09.
 *  Modified by Chensong Zhang on 03/30/2009
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file multigrid.c
 *  \brief Abstract multigrid cycle
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "precond.h"
#include "matvec.h"

/**
 * \fn void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P,
 *                    int level, int levelNum, int smoother, int m1, int m2, int m0)
 * \brief Solve Ax=b by abstract multigrid cycle
 *
 * \param *A pointer to stiffness matrix of levelNum levels
 * \param *b pointer to the dvector of right hand side term
 * \param *x pointer to the dvector of dofs
 * \param *R pointer to resriction operator array of levelNum levels
 * \param *P pointer to interpolation operator array of levelNum levels
 * \param level current level
 * \param levelNum total level num of grid
 * \param smoother smoother type
 * \param m1 pre-smoothing times
 * \param m2 post-smoothing times
 * \param m0 times of correction on coarse grid
 */
void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P, 
							 int level, int levelNum, int smoother, int m1, int m2, int m0)
{
	if(level<levelNum-1)
	{
		int i, j;
		dvector newb, e;
		int Rsize=R[level].row, Asize=A[level].row;
		double *r = (double*)calloc(Asize, sizeof(double));
		
		create_dvector(Rsize,&newb);
		create_dvector(Rsize,&e);

		/** pre smoothing */
		if (smoother == GS) {
			gs(x, 0, A[level].row-1, 1, &A[level], b, m1);
		}
		else if (smoother == JACOBI) {
			jacobi(x, 0, A[level].row-1, 1, &A[level], b, m1);
		}
		else if (smoother == SGS) {
			gs(x, 0, A[level].row-1, 1,  &A[level], b, m1);
			gs(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		
		/** form residual r = b - A x */
		copy_array(Asize,b->val,r);
		sparse_mv(-1.0,&A[level],x->val,r);
		
		/** restriction */
		for(i=0;i<R[level].row;i++) {
			for(j=R[level].IA[i];j<R[level].IA[i+1];j++) {
				newb.val[i]+=R[level].val[j]*r[R[level].JA[j]];
			}
		}
				
		/** call MG recursively: m0 = 1 for V cycle, m0 = 2 for W cycle */
		for(i=0; i<m0; i++)
			multigrid(A, &newb, &e, R, P, level+1, levelNum, smoother, m1, m2, m0);
		
		/** prolongation */
		for(i=0;i<P[level].row;i++) {
			for(j=P[level].IA[i];j<P[level].IA[i+1];j++) {
				x->val[i]+=P[level].val[j]*e.val[P[level].JA[j]];
			}
		}
		
		/** post smoothing */
		if (smoother == GS) {
			gs(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		else if (smoother == JACOBI) {
			jacobi(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		else if (smoother == SGS) {
			gs(x, 0, A[level].row-1, 1,  &A[level], b, m1);
			gs(x, A[level].row-1, 0, -1, &A[level], b, m2);
		}
		
		free_dvector(&newb);
		free_dvector(&e);
		free(r);
	}
	else // coarsest level solver
	{
		int MaxNumIt = 100; 
		double CoarseTol = 1e-10;
		pcg(&A[level], b, x, MaxNumIt, CoarseTol, NULL, 0);
	}

}

/**
* \fn void multigridvar(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P,
*                    int level, int levelNum, int smoother, int *m1, int *m2, int m0)
* \brief Solve Ax=b by abstract multigrid cycle with variable smoothing
*
* \param *A pointer to stiffness matrix of levelNum levels
* \param *b pointer to the dvector of right hand side term
* \param *x pointer to the dvector of dofs
* \param *R pointer to resriction operator array of levelNum levels
* \param *P pointer to interpolation operator array of levelNum levels
* \param level current level
* \param levelNum total level num of grid
* \param smoother smoother type
* \param *m1 pointer to pre-smoothing times
* \param *m2 pointer to post-smoothing times
* \param m0 times of correction on coarse grid
*/
void multigridvar(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P,
	int level, int levelNum, int smoother, int *m1, int *m2, int m0)
{
	if (level<levelNum - 1)
	{
		int i, j;
		dvector newb, e;
		int Rsize = R[level].row, Asize = A[level].row;
		double *r = (double*)calloc(Asize, sizeof(double));

		create_dvector(Rsize, &newb);
		create_dvector(Rsize, &e);

		/** pre smoothing */
		if (smoother == GS) {
			gs(x, 0, A[level].row - 1, 1, &A[level], b, m1[level]);
		}
		else if (smoother == JACOBI) {
			jacobi(x, 0, A[level].row - 1, 1, &A[level], b, m1[level]);
		}
		else if (smoother == SGS) {
			gs(x, 0, A[level].row - 1, 1, &A[level], b, m1[level]);
			gs(x, A[level].row - 1, 0, -1, &A[level], b, m2[level]);
		}

		/** form residual r = b - A x */
		copy_array(Asize, b->val, r);
		sparse_mv(-1.0, &A[level], x->val, r);

		/** restriction */
		for (i = 0; i<R[level].row; i++) {
			for (j = R[level].IA[i]; j<R[level].IA[i + 1]; j++) {
				newb.val[i] += R[level].val[j] * r[R[level].JA[j]];
			}
		}

		/** call MG recursively: m0 = 1 for V cycle, m0 = 2 for W cycle */
		for (i = 0; i<m0; i++)
			multigridvar(A, &newb, &e, R, P, level + 1, levelNum, smoother, m1, m2, m0);

		/** prolongation */
		for (i = 0; i<P[level].row; i++) {
			for (j = P[level].IA[i]; j<P[level].IA[i + 1]; j++) {
				x->val[i] += P[level].val[j] * e.val[P[level].JA[j]];
			}
		}

		/** post smoothing */
		if (smoother == GS) {
			gs(x, A[level].row - 1, 0, -1, &A[level], b, m2[level]);
		}
		else if (smoother == JACOBI) {
			jacobi(x, A[level].row - 1, 0, -1, &A[level], b, m2[level]);
		}
		else if (smoother == SGS) {
			gs(x, 0, A[level].row - 1, 1, &A[level], b, m1[level]);
			gs(x, A[level].row - 1, 0, -1, &A[level], b, m2[level]);
		}

		free_dvector(&newb);
		free_dvector(&e);
		free(r);
	}
	else // coarsest level solver
	{
		int MaxNumIt = 100;
		double CoarseTol = 1e-10;
		pcg(&A[level], b, x, MaxNumIt, CoarseTol, NULL, 0);
	}

}

/**
* \fn double mgvVectorP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m)
* \brief V-cycle geometric multigrid method for linear elasticity discretzed by conforming linear finite element method
*        There is no Lame constants for this linear elasticity, i.e. lambda=0, mu=0.5
*
* \param *A pointer to stiffness matrix of levelNum levels
* \param *b pointer to the dvector of right hand side term
* \param *x pointer to the dvector of dofs
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *nondirichlet pointer to the index of nondirichlet nodes
* \param *index pointer to the transpose of dirichlet and nondirichlet nodes
* \param levelNum total level num of grid
* \param smoother smoother type
* \param m smoothing times
*/
double mgvVectorP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m)
{
	int i, i1, l;
	dvector *r, *e;
	r = (dvector*)calloc(levelNum, sizeof(dvector));
	e = (dvector*)calloc(levelNum, sizeof(dvector));
	for (i = 0; i < levelNum; i++)
	{
		create_dvector(A[i].row, &r[i]);
		create_dvector(A[i].row, &e[i]);
	}

	/** form residual r = b - A x */
	copy_dvector(b, &r[levelNum - 1]);
	sparse_mv(-1.0, &A[levelNum - 1], x->val, r[levelNum - 1].val);

	// Loop down through the levels: "\"
	for (l = levelNum - 1; l>0; l--)
	{
		// Pre-smoothing
		if (smoother == GS) {
			gs(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
		}
		else if (smoother == JACOBI) {
			jacobi(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
		}
		else if (smoother == SGS) {
			gs(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
			gs(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}

		// Restrict the residual and keep it
		dvector temp;
		create_dvector(r[l].row, &temp);
		/** form residual temp = r - A e */
		copy_dvector(&r[l], &temp);
		sparse_mv(-1.0, &A[l], e[l].val, temp.val);

		restrictionPTvector2d(&edgesTran[l], &index[l], &nondirichlet[l - 1], &temp, &r[l - 1]);
		free_dvector(&temp);
	}

	// Solve on the coarse grid, using Conjugate Gradient method
	int MaxNumIt = 100;
	double CoarseTol = 1e-10;
	pcg(&A[0], &r[0], &e[0], MaxNumIt, CoarseTol, NULL, 0);

	// Loop up: "/"
	for (l = 1; l<levelNum; l++)
	{
		// Correct
		dvector temp;
		create_dvector(e[l].row, &temp);
		interpolationPvector2d(&edges[l - 1], &isInNode[l - 1], &index[l - 1], &nondirichlet[l], nodeCEdge, index[l - 1].row / 2, &e[l - 1], &temp);
		axpy_dvector(1.0, &temp, &e[l]);
		free_dvector(&temp);

		// Post smoothing
		if (smoother == GS) {
			gs(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}
		else if (smoother == JACOBI) {
			jacobi(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}
		else if (smoother == SGS) {
			gs(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
			gs(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}
	}

	// Update
	axpy_dvector(1.0, &e[levelNum - 1], x);

	// computer error
	copy_dvector(b, &r[levelNum - 1]);
	sparse_mv(-1.0, &A[levelNum - 1], x->val, r[levelNum - 1].val);

	double aa = dot_dvector(&r[levelNum - 1], &r[levelNum - 1]);
	double bb = dot_dvector(b, b);
	double relres = sqrt(aa / bb);

	for (i = 0; i < levelNum; i++)
	{
		free_dvector(&r[i]);
		free_dvector(&e[i]);
	}
	free(r); free(e);

	return relres;
}

/**
* \fn double mgvP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m)
* \brief V-cycle geometric multigrid method for linear elasticity discretzed by conforming linear finite element method
*        There is no Lame constants for this linear elasticity, i.e. lambda=0, mu=0.5
*
* \param *A pointer to stiffness matrix of levelNum levels
* \param *b pointer to the dvector of right hand side term
* \param *x pointer to the dvector of dofs
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
* \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *nondirichlet pointer to the index of nondirichlet nodes
* \param *index pointer to the transpose of dirichlet and nondirichlet nodes
* \param levelNum total level num of grid
* \param smoother smoother type
* \param m smoothing times
*/
double mgvP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m)
{
	int i, i1, l;
	dvector *r, *e;
	r = (dvector*)calloc(levelNum, sizeof(dvector));
	e = (dvector*)calloc(levelNum, sizeof(dvector));
	for (i = 0; i < levelNum; i++)
	{
		create_dvector(A[i].row, &r[i]);
		create_dvector(A[i].row, &e[i]);
	}

	/** form residual r = b - A x */
	copy_dvector(b, &r[levelNum - 1]);
	sparse_mv(-1.0, &A[levelNum - 1], x->val, r[levelNum - 1].val);

	// Loop down through the levels: "\"
	for (l = levelNum - 1; l>0; l--)
	{
		// Pre-smoothing
		if (smoother == GS) {
			gs(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
		}
		else if (smoother == JACOBI) {
			jacobi(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
		}
		else if (smoother == SGS) {
			gs(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
			gs(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}

		// Restrict the residual and keep it
		dvector temp;
		create_dvector(r[l].row, &temp);
		/** form residual temp = r - A e */
		copy_dvector(&r[l], &temp);
		sparse_mv(-1.0, &A[l], e[l].val, temp.val);

		restrictionPT2d(&edgesTran[l], &index[l], &nondirichlet[l - 1], &temp, &r[l - 1]);
		free_dvector(&temp);
	}

	// Solve on the coarse grid, using Conjugate Gradient method
	int MaxNumIt = 100;
	double CoarseTol = 1e-10;
	pcg(&A[0], &r[0], &e[0], MaxNumIt, CoarseTol, NULL, 0);
	
	// Loop up: "/"
	for (l = 1; l<levelNum; l++)
	{
		// Correct
		dvector temp;
		create_dvector(e[l].row, &temp);
		interpolationP2d(&edges[l - 1], &isInNode[l - 1], &index[l - 1], &nondirichlet[l], nodeCEdge, index[l - 1].row , &e[l - 1], &temp);
		axpy_dvector(1.0, &temp, &e[l]);
		free_dvector(&temp);

		// Post smoothing
		if (smoother == GS) {
			gs(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}
		else if (smoother == JACOBI) {
			jacobi(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}
		else if (smoother == SGS) {
			gs(&e[l], 0, e[l].row - 1, 1, &A[l], &r[l], m);
			gs(&e[l], e[l].row - 1, 0, -1, &A[l], &r[l], m);
		}
	}
	
	// Update
	axpy_dvector(1.0, &e[levelNum - 1], x);

	// computer error
	copy_dvector(b, &r[levelNum - 1]);
	sparse_mv(-1.0, &A[levelNum - 1], x->val, r[levelNum - 1].val);

	double aa = dot_dvector(&r[levelNum - 1], &r[levelNum - 1]);
	double bb = dot_dvector(b, b);
	double relres = sqrt(aa / bb);

	for (i = 0; i < levelNum; i++)
	{
		free_dvector(&r[i]);
		free_dvector(&e[i]);
	}
	free(r); free(e);

	return relres;
}

/**
* \fn void interpolationPvector2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe)
* \brief the interpolation operator from coarse grid to fine grid 
* \param *Cedges pointer to edges on coarse grid: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
* \param *CisInNode pointer to boundary information of nodes on coarse grid: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *Cindex pointer to the transpose of dirichlet and nondirichlet nodes on coarse grid
* \param *Fnondirichlet pointer to the index of nondirichlet nodes on fine grid
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param cnn the number of vertices on coarse grid
* \param *e pointer to the dvector on coarse grid
* \param *Pe pointer to the dvector  on fine grid
*/
void interpolationPvector2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe)
{
	int i, j, i1, j1, k;
	int cdof = e->row / 2;
	int fdof = Pe->row / 2;

	for (i1 = 0; i1<fdof; i1++)
	{
		i = i1;
		if (Fnondirichlet != NULL)
			i = Fnondirichlet->val[i];

		if (i<cnn) // case i1 is on the level-1(coarse) grid 
		{
			j1 = i;
			if (Cindex != NULL)
				j1 = Cindex->val[j1];

			Pe->val[i1] = e->val[j1];
			Pe->val[i1 + fdof] = e->val[j1 + cdof];
		}
		else // case i1 is on the level(fine) grid 
		{
			Pe->val[i1] = 0;
			Pe->val[i1 + fdof] = 0;
			for (k = 0; k < 2; k++)
			{
				j1 = Cedges->val[nodeCEdge->val[i]][k];
				if (CisInNode != NULL)
				{
					if (CisInNode->val[j1] != -1)
					{
						if (Cindex != NULL)
							j1 = Cindex->val[j1];

						Pe->val[i1] += e->val[j1] / 2;
						Pe->val[i1 + fdof] += e->val[j1 + cdof] / 2;
					}
				}
				else
				{
					if (Cindex != NULL)
						j1 = Cindex->val[j1];
					
					Pe->val[i1] += e->val[j1] / 2;
					Pe->val[i1 + fdof] += e->val[j1 + cdof] / 2;
				}
				
			}
		}
	}
}

/**
* \fn void interpolationP2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe)
* \brief the interpolation operator from coarse grid to fine grid
* \param *Cedges pointer to edges on coarse grid: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
* \param *CisInNode pointer to boundary information of nodes on coarse grid: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *Cindex pointer to the transpose of dirichlet and nondirichlet nodes on coarse grid
* \param *Fnondirichlet pointer to the index of nondirichlet nodes on fine grid
* \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
* \param cnn the number of vertices on coarse grid
* \param *e pointer to the dvector on coarse grid
* \param *Pe pointer to the dvector  on fine grid
*/
void interpolationP2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe)
{
	int i, j, i1, j1, k;
	int cdof = e->row;
	int fdof = Pe->row;

	for (i1 = 0; i1<fdof; i1++)
	{
		i = i1;
		if (Fnondirichlet != NULL)
			i = Fnondirichlet->val[i];

		if (i<cnn) // case i1 is on the level-1(coarse) grid 
		{
			j1 = i;
			if (Cindex != NULL)
				j1 = Cindex->val[j1];

			Pe->val[i1] = e->val[j1];
		}
		else // case i1 is on the level(fine) grid 
		{
			Pe->val[i1] = 0;
			for (k = 0; k < 2; k++)
			{
				j1 = Cedges->val[nodeCEdge->val[i]][k];
				if (CisInNode != NULL)
				{
					if (CisInNode->val[j1] != -1)
					{
						if (Cindex != NULL)
							j1 = Cindex->val[j1];

						Pe->val[i1] += e->val[j1] / 2;
					}
				}
				else
				{
					if (Cindex != NULL)
						j1 = Cindex->val[j1];

					Pe->val[i1] += e->val[j1] / 2;
				}

			}
		}
	}
}

/**
* \fn void restrictionPTvector2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr)
* \brief the restriction operator from fine grid to coarse grid
* \param *FedgesTran pointer to the tranpose of edges on fine grid, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *Findex pointer to the transpose of dirichlet and nondirichlet nodes on fine grid
* \param *Cnondirichlet pointer to the index of nondirichlet nodes on coarse grid
* \param *r pointer to the dvector on fine grid
* \param *PTr pointer to the dvector  on coarse grid
*/
void restrictionPTvector2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr)
{
	int i, j, i1, j1;
	int cdof = PTr->row / 2;
	int fdof = r->row / 2;
	for (i1 = 0; i1<cdof; i1++)
	{
		if (Cnondirichlet != NULL && Findex != NULL)
		{ 
			i = Cnondirichlet->val[i1];
			j1 = Findex->val[i];
		}
		else
		{
			i = i1;
			j1 = i;
		}
		PTr->val[i1] = r->val[j1];
		PTr->val[i1 + cdof] = r->val[j1 + fdof];
		for (j = FedgesTran->IA[i]; j<FedgesTran->IA[i + 1]; j++)
		{
			j1 = FedgesTran->val[j];
			if(Findex != NULL)
				j1 = Findex->val[j1];

			PTr->val[i1] += r->val[j1] / 2;
			PTr->val[i1 + cdof] += r->val[j1 + fdof] / 2;
		}
	}
}

/**
* \fn void restrictionPT2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr)
* \brief the restriction operator from fine grid to coarse grid
* \param *FedgesTran pointer to the tranpose of edges on fine grid, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *Findex pointer to the transpose of dirichlet and nondirichlet nodes on fine grid
* \param *Cnondirichlet pointer to the index of nondirichlet nodes on coarse grid
* \param *r pointer to the dvector on fine grid
* \param *PTr pointer to the dvector  on coarse grid
*/
void restrictionPT2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr)
{
	int i, j, i1, j1;
	int cdof = PTr->row;
	int fdof = r->row;
	for (i1 = 0; i1<cdof; i1++)
	{
		if (Cnondirichlet != NULL && Findex != NULL)
		{
			i = Cnondirichlet->val[i1];
			j1 = Findex->val[i];
		}
		else
		{
			i = i1;
			j1 = i;
		}
		PTr->val[i1] = r->val[j1];
		for (j = FedgesTran->IA[i]; j<FedgesTran->IA[i + 1]; j++)
		{
			j1 = FedgesTran->val[j];
			if (Findex != NULL)
				j1 = Findex->val[j1];

			PTr->val[i1] += r->val[j1] / 2;
		}
	}
}

/**
* \fn void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P,
*                    int level, int levelNum, int smoother, int m1, int m2, int m0)
* \brief V-cycle geometric multigrid method for linear elasticity discretzed by conforming linear finite element method
*        There is no Lame constants for this linear elasticity, i.e. lambda=0, mu=0.5
*
* \param *A pointer to stiffness matrix of levelNum levels
* \param *b pointer to the dvector of right hand side term
* \param *x pointer to the dvector of dofs
* \param *R pointer to resriction operator array of levelNum levels
* \param *P pointer to interpolation operator array of levelNum levels
* \param level current level
* \param levelNum total level num of grid
* \param smoother smoother type
* \param m smoothing times
*/
double mgvVectorP1_solveOld(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *nondirichlet, ivector *index, int levelNum, int m)
{
	int i, i1, l;
	dvector *r, *r1, *e, *e1;
	r = (dvector*)calloc(levelNum, sizeof(dvector));
	r1 = (dvector*)calloc(levelNum, sizeof(dvector));
	e = (dvector*)calloc(levelNum, sizeof(dvector));
	e1 = (dvector*)calloc(levelNum, sizeof(dvector));
	for (i = 0; i < levelNum; i++)
	{
		create_dvector(b[i].row, &r[i]);
		create_dvector(b1[i].row, &r1[i]);
		create_dvector(b[i].row, &e[i]);
		create_dvector(b1[i].row, &e1[i]);
	}

	/** form residual r = b - A x */
	copy_dvector(&b[levelNum - 1], &r[levelNum - 1]);
	sparse_mv(-1.0, &A[levelNum - 1], x->val, r[levelNum - 1].val);

	// achiveve r1[levelNum-1] due to dirichlet boundary condition
	for (i1 = 0; i1<r1[levelNum - 1].row; i1++)
	{
		i = nondirichlet[levelNum - 1].val[i1];
		r1[levelNum - 1].val[i1] = r[levelNum - 1].val[i];
	}

	// Loop down through the levels: "\"
	for (l = levelNum - 1; l>0; l--)
	{
		// Pre-smoothing
		gs(&e1[l], 0, e1[l].row - 1, 1, &A11[l], &r1[l], m);
		for (i = 0; i<e1[l].row; i++)
		{
			e[l].val[nondirichlet[l].val[i]] = e1[l].val[i];
		}
		// Restrict the residual and keep it
		dvector temp;
		create_dvector(r[l].row, &temp);
		/** form residual temp = r - A e */
		copy_dvector(&r[l], &temp);
		sparse_mv(-1.0, &A[l], e[l].val, temp.val);

		restrictionPTvector2dOld(&edgesTran[l], &temp, &r[l - 1]);
		free_dvector(&temp);
		// achieve r1[l-1] due to dirichlet boundary condition
		for (i1 = 0; i1<r1[l - 1].row; i1++)
		{
			i = nondirichlet[l - 1].val[i1];
			r1[l - 1].val[i1] = r[l - 1].val[i];
		}
	}

	// Solve on the coarse grid, using Conjugate Gradient method
	int MaxNumIt = 100;
	double CoarseTol = 1e-10;
	pcg(&A11[0], &r1[0], &e1[0], MaxNumIt, CoarseTol, NULL, 0);
	for (i = 0; i<e1[0].row; i++)
	{
		e[0].val[nondirichlet[0].val[i]] = e1[0].val[i];
	}

	// Loop up: "/"
	for (l = 1; l<levelNum; l++)
	{
		// Correct
		dvector temp;
		create_dvector(e[l].row, &temp);
		interpolationPvector2dOld(&edges[l - 1], nodeCEdge, &e[l - 1], &temp);
		axpy_dvector(1.0, &temp, &e[l]);
		free_dvector(&temp);

		// achiveve e1[l] due to dirichlet boundary condition
		for (i1 = 0; i1<e1[l].row; i1++)
		{
			i = nondirichlet[l].val[i1];
			e1[l].val[i1] = e[l].val[i];
		}
		// Post smoothing
		gs(&e1[l], e1[l].row - 1, 0, -1, &A11[l], &r1[l], m);
		for (i = 0; i<e1[l].row; i++)
		{
			e[l].val[nondirichlet[l].val[i]] = e1[l].val[i];
		}
	}

	// Update
	axpy_dvector(1.0, &e[levelNum - 1], x);

	// computer error
	copy_dvector(&b[levelNum - 1], &r[levelNum - 1]);
	sparse_mv(-1.0, &A[levelNum - 1], x->val, r[levelNum - 1].val);
	for (i1 = 0; i1<r1[levelNum - 1].row; i1++)
	{
		i = nondirichlet[levelNum - 1].val[i1];
		r1[levelNum - 1].val[i1] = r[levelNum - 1].val[i];
	}

	double aa = dot_dvector(&r1[levelNum - 1], &r1[levelNum - 1]);
	double bb = dot_dvector(&b1[levelNum - 1], &b1[levelNum - 1]);
	double relres = sqrt(aa / bb);

	for (i = 0; i < levelNum; i++)
	{
		free_dvector(&r[i]);
		free_dvector(&r1[i]);
		free_dvector(&e[i]);
		free_dvector(&e1[i]);
	}
	free(r); free(r1); free(e); free(e1);

	return relres;
}

/**
* \fn void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P,
*                    int level, int levelNum, int smoother, int m1, int m2, int m0)
* \brief the interpolation operator from coarse grid level-1 to fine grid level
*  Input: Cedges, nodeCEdge, e
*  Output: Pe
* \param *A pointer to stiffness matrix of levelNum levels
* \param *b pointer to the dvector of right hand side term
* \param *x pointer to the dvector of dofs
* \param *R pointer to resriction operator array of levelNum levels
* \param *P pointer to interpolation operator array of levelNum levels
* \param level current level
* \param levelNum total level num of grid
* \param smoother smoother type
* \param m smoothing times
*/
void interpolationPvector2dOld(EDGE *Cedges, ivector *nodeCEdge, dvector *e, dvector *Pe)
{
	int i, ni[2];
	int cdof = e->row / 2;
	int fdof = Pe->row / 2;

	for (i = 0; i<fdof; i++)
	{
		if (i<cdof) // case i is on the level-1(coarse) grid 
		{
			Pe->val[i] = e->val[i];
			Pe->val[i + fdof] = e->val[i + cdof];
		}
		else // case i is on the level(fine) grid 
		{
			ni[0] = Cedges->val[nodeCEdge->val[i]][0];
			ni[1] = Cedges->val[nodeCEdge->val[i]][1];
			Pe->val[i] = (e->val[ni[0]] + e->val[ni[1]]) / 2;
			Pe->val[i + fdof] = (e->val[ni[0] + cdof] + e->val[ni[1] + cdof]) / 2;
		}
	}
}

/**
* \fn void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P,
*                    int level, int levelNum, int smoother, int m1, int m2, int m0)
* \brief the restriction operator from fine grid level to coarse grid level-1
*  Input: edgesTran, level, r
*  Output: PTr
* \param *A pointer to stiffness matrix of levelNum levels
* \param *b pointer to the dvector of right hand side term
* \param *x pointer to the dvector of dofs
* \param *R pointer to resriction operator array of levelNum levels
* \param *P pointer to interpolation operator array of levelNum levels
* \param level current level
* \param levelNum total level num of grid
* \param smoother smoother type
* \param m smoothing times
*/
void restrictionPTvector2dOld(iCSRmat *edgesTran, dvector *r, dvector *PTr)
{
	int i, j;
	int cdof = PTr->row / 2;
	int fdof = r->row / 2;
	for (i = 0; i<cdof; i++)
	{
		PTr->val[i] = r->val[i];
		PTr->val[i + cdof] = r->val[i + fdof];
		for (j = edgesTran->IA[i]; j<edgesTran->IA[i + 1]; j++)
		{
			PTr->val[i] += r->val[edgesTran->val[j]] / 2;
			PTr->val[i + cdof] += r->val[edgesTran->val[j] + fdof] / 2;
		}
	}
}