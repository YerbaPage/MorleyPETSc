/*
 *  smoother.c
 *
 *  Created by Xuehai Huang on 10/30/2008
 *  Modified by Chensong Zhang on 03/27/2009
 *  Copyright 2008 PSU. All rights reserved.
 *
 */
 
/*! \file smoother.c
 *  \brief Smoothers for multigrid method
 */
  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn void gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief Gauss-Seidel method  as the smoother in solving Au=b with multigrid method
 * \param u initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1 the index to begin with
 * \param i_n the index to end
 * \param s the step (s=1: forward, s=-1: backward)
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param L number of iterations
 * \return void
 *
 * Gauss-Seidel smoother (Smoother_Type = 1)
 */
void gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
{
	int i,j,k;
	double t,d;
	int *ia=A->IA,*ja=A->JA;
	double *aj=A->val,*bval=b->val,*uval=u->val;
	
	if (s > 0) {
		while (L--) {
			for (i=i_1;i<=i_n;i+=s) {
				t=bval[i];
				int begin_row=ia[i],end_row=ia[i+1]-1;
				for (k=begin_row;k<=end_row;k++) {
					j=ja[k];
					if (i!=j) t-=aj[k]*uval[j];
					else d=aj[k];
				} // end for k
				
				uval[i]=t/d;
			} // end for i
		} // end while		
	} // if s
	else {		
		while (L--) {
			for (i=i_1;i>=i_n;i+=s) {
				t=bval[i];
				int begin_row=ia[i],end_row=ia[i+1]-1;
				for (k=begin_row;k<=end_row;k++) {
					j=ja[k];
					if (i!=j) t-=aj[k]*uval[j];
					else d=aj[k];
				} // end for k
				
				uval[i]=t/d;
			} // end for i
		} // end while		
  } // end if
	return;

/*	int i,j,k,l=0;
	double t, d;
	
	for (l=0;l<L;l++)
	{
		for (i=i_1;i<=i_n;i+=s)
		{
			t=b->val[i];
			for (k=A->IA[i];k<A->IA[i+1];k++)
			{
				j=A->JA[k];
				if (i!=j)
					t-=A->val[k]*u->val[j];
				else
					d=A->val[k];
			}
			u->val[i]=t/d;
		}
	}*/
}

/**
 * \fn void jacobi(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
 * \brief Jacobi method as the smoother in solving Au=b with multigrid method
 * \param u initial guess and the new approximation to the solution obtained after L Gauss-Seidel steps
 * \param i_1 the index to begin with
 * \param i_n the index to end
 * \param s the step
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param L number of iterations
 * \return void
 *
 * Jacobi smoother (Smoother_Type = 2)
 */
void jacobi(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L)
{
	int i,j,k,l=0;
	dvector t, d;
	create_dvector(u->row, &t);
	create_dvector(u->row, &d);

	for (l=0;l<L;l++)
	{
		for (i=i_1;i<=i_n;i+=s)
		{
			t.val[i]=b->val[i];
			for (k=A->IA[i];k<A->IA[i+1];k++)
			{
				j=A->JA[k];
				if(i!=j)
 					t.val[i]-=A->val[k]*u->val[j];
				else
					d.val[i]=A->val[k];
			}
		}
		for (i=i_1;i<=i_n;i+=s) u->val[i]=t.val[i]/d.val[i];
	}

	free_dvector(&t);
	free_dvector(&d);

}

/**
* \fn void mulschwarz(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, dOBDmat *B, int L)
* \brief Multiplicative Schwarz smoother in solving Au=b with multigrid method
* \param u initial guess and the new approximation to the solution obtained after L smoothing steps
* \param i_1 the index to begin with
* \param i_n the index to end
* \param s the step (s=1: forward, s=-1: backward)
* \param *A pointer to stiffness matrix
* \param *B pointer to smoother of block diagnoal type
* \param *b pointer to right hand side
* \param L number of iterations
* \return void
*/
void mulschwarz(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, dOBDmat *B, int L)
{
	int i, j, k, i1, j1, ii, jj;
	double t, d;
	int *ia = A->IA, *ja = A->JA;
	double *aj = A->val, *bval = b->val, *uval = u->val;
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	double *r;

	if (s > 0) {
		while (L--) {
			for (i = i_1; i <= i_n; i += s) {
				blki = blk + i;
				ri = rindices + i;
				ci = cindices + i;
				if (blki->row < 1 || blki->col < 1)
					continue;

				// r = b - Au for block i
				r = (double*)malloc(ci->row *sizeof(double));
				for (j1 = 0; j1 < ci->row; j1++)
				{
					jj = ci->val[j1];
					r[j1] = bval[jj];
					for (k = ia[jj]; k < ia[jj + 1]; k++) {
						j = ja[k];
						r[j1] -= aj[k] * uval[j];
					} // end for k
				}
				// u = u + Br = u + B(b-Au) for block i
				for (i1 = 0; i1 < ri->row; i1++)
				{
					ii = ri->val[i1];
					for (j1 = 0; j1 < ci->row; j1++)
						uval[ii] += blki->val[i1][j1] * r[j1];
				}
				free(r);
			} // end for i
		} // end while		
	} // if s
	else {
		while (L--) {
			for (i = i_1; i >= i_n; i += s) {
				blki = blk + i;
				ri = rindices + i;
				ci = cindices + i;
				if (blki->row < 1 || blki->col < 1)
					continue;
				
				// r = b - Au for block i
				r = (double*)malloc(ci->row * sizeof(double));
				for (j1 = 0; j1 < ci->row; j1++)
				{
					jj = ci->val[j1];
					r[j1] = bval[jj];
					for (k = ia[jj]; k < ia[jj + 1]; k++) {
						j = ja[k];
						r[j1] -= aj[k] * uval[j];
					} // end for k
				}
				// u = u + Br = u + B(b-Au) for block i
				for (i1 = 0; i1 < ri->row; i1++)
				{
					ii = ri->val[i1];
					for (j1 = 0; j1 < ci->row; j1++)
						uval[ii] += blki->val[i1][j1] * r[j1];
				}
				free(r);
			} // end for i
		} // end while		
	} // end if

	return;
}

/**
* \fn void getSchwarzblocks_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks based on the vertex patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elements pointer to the structure of the triangulation
* \param nvertices number of vertices
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocks_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
{
	int i, j, k;
	int nb = nvertices;
	iCSRmat nodeElem;
	int element, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	getTransposeOfELEMENT(elements, &nodeElem, elements->col, nvertices);

	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (j = nodeElem.IA[i]; j<nodeElem.IA[i + 1]; j++)
		{
			element = nodeElem.JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element][k];
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
			}
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free_icsr_matrix(&nodeElem);

	return;
}

/**
* \fn void getSchwarzblocks_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the edge patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocks_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF)
{
	int i, j, k;
	int nb = edges->row;
	int element, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (j = 0; j<2; j++)
		{
			element = edges->val[i][2 + j];
			if (element == -1)
				continue;

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element][k];
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);

	return;
}

/**
* \fn void getSchwarzblocks_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the edge patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elements pointer to the structure of the triangulation
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param nvertices number of vertices
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocks_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF)
{
	int i, j, k, ii, jj;
	int nb = edges->row;
	iCSRmat nodeElem;
	int element, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	getTransposeOfELEMENT(elements, &nodeElem, elements->col, nvertices);

	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (jj = 0; jj<2; jj++)
		{
			ii = edges->val[i][jj];
			for (j = nodeElem.IA[ii]; j<nodeElem.IA[ii + 1]; j++)
			{
				element = nodeElem.JA[j];

				for (k = 0; k<elementDOF->col; k++)
				{
					node = elementDOF->val[element][k];
					ni = index->val[node];
					if (nfFlag->val[node] != 1 && mask[ni] == -1)
					{
						mask[ni] = istart;
						istart = ni;
						count++;
					}
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);
	free_icsr_matrix(&nodeElem);

	return;
}

/**
* \fn void getSchwarzblocks_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks based on the element patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocks_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF)
{
	int i, j, k;
	int nb = elementEdge->row;
	int element, edge, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (j = 0; j<elementEdge->col; j++)
		{
			edge = elementEdge->val[i][j];
			element = edges->val[edge][2] + edges->val[edge][3] - i;
			if (element == -1)
				continue;

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element][k];
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);

	return;
}

/**
* \fn void getSchwarzblocks_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the element-vertex patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elements pointer to the structure of the triangulation
* \param nvertices number of vertices
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocks_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
{
	int i, j, k, ii, jj;
	int nb = elements->row;
	iCSRmat nodeElem;
	int element, edge, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	getTransposeOfELEMENT(elements, &nodeElem, elements->col, nvertices);

	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (jj = 0; jj<elements->col; jj++)
		{
			ii = elements->val[i][jj];
			for (j = nodeElem.IA[ii]; j<nodeElem.IA[ii + 1]; j++)
			{
				element = nodeElem.JA[j];

				for (k = 0; k<elementDOF->col; k++)
				{
					node = elementDOF->val[element][k];
					ni = index->val[node];
					if (nfFlag->val[node] != 1 && mask[ni] == -1)
					{
						mask[ni] = istart;
						istart = ni;
						count++;
					}
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);
	free_icsr_matrix(&nodeElem);

	return;
}

/**
* \fn void getSchwarzblocksVec2_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the vertex patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elements pointer to the structure of the triangulation
* \param nvertices number of vertices
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocksVec2_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
{
	int i, j, k;
	int nb = nvertices;
	iCSRmat nodeElem;
	int element, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;
	
	getTransposeOfELEMENT(elements, &nodeElem, elements->col, nvertices);
	
	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;
	
	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (j = nodeElem.IA[i]; j<nodeElem.IA[i + 1]; j++)
		{
			element = nodeElem.JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element][k];
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
				node += elementDOF->dof;
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
			}
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
	
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);
	free_icsr_matrix(&nodeElem);
	
	return;
}

/**
* \fn void getSchwarzblocksVec2_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the edge patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocksVec2_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF)
{
	int i, j, k;
	int nb = edges->row;
	int element, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (j = 0; j<2; j++)
		{
			element = edges->val[i][2 + j];
			if (element == -1)
				continue;

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element][k];
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
				node += elementDOF->dof;
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);

	return;
}

/**
* \fn void getSchwarzblocksVec2_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the edge-vertex patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elements pointer to the structure of the triangulation
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param nvertices number of vertices
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocksVec2_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF)
{
	int i, j, k, ii, jj;
	int nb = edges->row;
	iCSRmat nodeElem;
	int element, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	getTransposeOfELEMENT(elements, &nodeElem, elements->col, nvertices);
	
	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (jj = 0; jj<2; jj++)
		{
			ii = edges->val[i][jj];
			for (j = nodeElem.IA[ii]; j<nodeElem.IA[ii + 1]; j++)
			{
				element = nodeElem.JA[j];

				for (k = 0; k<elementDOF->col; k++)
				{
					node = elementDOF->val[element][k];
					ni = index->val[node];
					if (nfFlag->val[node] != 1 && mask[ni] == -1)
					{
						mask[ni] = istart;
						istart = ni;
						count++;
					}
					node += elementDOF->dof;
					ni = index->val[node];
					if (nfFlag->val[node] != 1 && mask[ni] == -1)
					{
						mask[ni] = istart;
						istart = ni;
						count++;
					}
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);
	free_icsr_matrix(&nodeElem);

	return;
}

/**
* \fn void getSchwarzblocksVec2_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the element patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocksVec2_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF)
{
	int i, j, k;
	int nb = elementEdge->row;
	int element, edge, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (j = 0; j<elementEdge->col; j++)
		{
			edge = elementEdge->val[i][j];
			element = edges->val[edge][2] + edges->val[edge][3] - i;
			if (element == -1)
				continue;

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element][k];
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
				node += elementDOF->dof;
				ni = index->val[node];
				if (nfFlag->val[node] != 1 && mask[ni] == -1)
				{
					mask[ni] = istart;
					istart = ni;
					count++;
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);

	return;
}

/**
* \fn void getSchwarzblocksVec2_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
* \brief get Schwarz blocks in 2d vector based on the element-vertex patch
* \param *B pointer to Schwarz diagnoal blocks
* \param *A pointer to stiffness matrix
* \param *elements pointer to the structure of the triangulation
* \param nvertices number of vertices
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getSchwarzblocksVec2_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF)
{
	int i, j, k, ii, jj;
	int nb = elements->row;
	iCSRmat nodeElem;
	int element, edge, node, ni;
	ddenmat Aloc;
	int count;

	ivector *freenodes = &elementDOF->freenodes;
	//	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *index = &elementDOF->index;

	create_dobd_matrix(A->row, A->col, nb, B);
	ddenmat *blk = B->blk, *blki;
	ivector *rindices = B->rindices, *ri;
	ivector *cindices = B->cindices, *ci;

	getTransposeOfELEMENT(elements, &nodeElem, elements->col, nvertices);
	
	int *mask;
	int istart;
	mask = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		mask[i] = -1;

	for (i = 0; i<nb; i++)
	{
		blki = blk + i;
		ri = rindices + i;
		ci = cindices + i;

		count = 0;
		istart = -2;
		for (jj = 0; jj<elements->col; jj++)
		{
			ii = elements->val[i][jj];
			for (j = nodeElem.IA[ii]; j<nodeElem.IA[ii + 1]; j++)
			{
				element = nodeElem.JA[j];

				for (k = 0; k<elementDOF->col; k++)
				{
					node = elementDOF->val[element][k];
					ni = index->val[node];
					if (nfFlag->val[node] != 1 && mask[ni] == -1)
					{
						mask[ni] = istart;
						istart = ni;
						count++;
					}
					node += elementDOF->dof;
					ni = index->val[node];
					if (nfFlag->val[node] != 1 && mask[ni] == -1)
					{
						mask[ni] = istart;
						istart = ni;
						count++;
					}
				}
			}
		}

		if (count == 0)
		{
			blki->row = 0;
			blki->col = 0;
			blki->val = NULL;
			continue;
		}

		create_ivector(count, ri);
		create_ivector(count, ci);

		for (j = 0; j<count; j++)
		{
			ri->val[j] = istart;
			ci->val[j] = istart;
			istart = mask[istart];
			mask[ri->val[j]] = -1;
		}

		get_block_dden(A, ri->row, ci->row, ri->val, ci->val, &Aloc, mask);
		create_dden_matrix(ri->row, ci->row, blki);
		for (j = 0; j < blki->row; j++)
			blki->val[j][j] = 1;
		AxBrref(&Aloc, blki);
		free_dden_matrix(&Aloc);
	} // i

	free(mask);
	free_icsr_matrix(&nodeElem);

	return;
}
