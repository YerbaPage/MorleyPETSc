/*
 *  constrInterpolation.c
 *  classicAMG
 *
 *  Created by Xuehai Huang on 1/31/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */
 
/*! \file interpolation.c
 *  \brief Interpolation operators
 */
  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "header.h"
#include "matvec.h"
#include "checkmat.h"

 /**
 * \fn void interpP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg)
 * \brief the interpolation matrix from the 1st order Lagrange element to piecewise kth order polynomial in 2d
 * \param *P pointer to the vector-version interpolation matrix
 * \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
 * \param *elementDOFdg pointer to the relation between elements and degrees of freedom of the piecewise kth order polynomial
 */
void interpP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg)
{
	int i, j, ie, ii, k;
	int curnode;
	int dop = elementDOFdg->dop;

	if (elementDOFp1->dop != 1)
	{
		P = NULL;
		return;
	}

	P->row = elementDOFdg->dof;
	P->col = elementDOFp1->dof;
	P->IA = (int*)malloc((P->row + 1) * sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	P->IA[0] = 0;
	for (k = 0; k < elementDOFdg->row; k++)
	{
		if (dop == 0)
		{
			for (i = 0; i < elementDOFdg->col; i++) //  for each node
			{
				curnode = elementDOFdg->val[k][i];
				P->IA[curnode + 1] = elementDOFp1->col;
			}
			continue;
		}

		for (i = 0; i < 3; i++) //  for each vertex
		{
			curnode = elementDOFdg->val[k][i];
			P->IA[curnode + 1] = 1;
		}
		for (i = 3; i < 3 * dop; i++) //  for each node in edge
		{
			curnode = elementDOFdg->val[k][i];
			P->IA[curnode + 1] = 2;
		}
		for (i = 3 * dop; i < elementDOFdg->col; i++) //  for each node in element
		{
			curnode = elementDOFdg->val[k][i];
			P->IA[curnode + 1] = elementDOFp1->col;
		}
	}
	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];
	P->nnz = P->IA[P->row];

	// step 2P: Find the structure JA of the interpolation matrix P
	P->JA = (int*)malloc(P->nnz * sizeof(int));
	for (k = 0; k < elementDOFdg->row; k++)
	{
		if (dop == 0)
		{
			for (i = 0; i < elementDOFdg->col; i++) //  for each node
			{
				curnode = elementDOFdg->val[k][i];
				for (j = 0; j < elementDOFp1->col; j++)
				{
					P->JA[P->IA[curnode] + j] = elementDOFp1->val[k][j];
				}
			}
			continue;
		}

		for (i = 0; i < 3; i++) //  for each vertex
		{
			curnode = elementDOFdg->val[k][i];
			P->JA[P->IA[curnode]] = elementDOFp1->val[k][i];
		}
		for (i = 3; i < 3 * dop; i++) //  for each node in edge
		{
			curnode = elementDOFdg->val[k][i];
			ie = (i - 3) / (dop - 1);
			for (j = 0; j < 2; j++)
			{
				P->JA[P->IA[curnode] + j] = elementDOFp1->val[k][(ie + j) % 3];
			}
		}
		for (i = 3 * dop; i < elementDOFdg->col; i++) //  for each node in element
		{
			curnode = elementDOFdg->val[k][i];
			for (j = 0; j < elementDOFp1->col; j++)
			{
				P->JA[P->IA[curnode] + j] = elementDOFp1->val[k][j];
			}
		}
	}

	// step 3P: Loop element by element and compute the actual entries storing them in P
	P->val = (double*)malloc(P->nnz * sizeof(double));
	for (k = 0; k < elementDOFdg->row; k++)
	{
		if (elementDOFdg->dop == 0)
		{
			curnode = elementDOFdg->val[k][0];
			for (j = 0; j < 3; j++)
				P->val[P->IA[curnode] + j] = 1.0 / 3.0;
			continue;
		}

		for (i = 0; i < 3; i++) //  for each vertex
		{
			curnode = elementDOFdg->val[k][i];
			P->val[P->IA[curnode]] = 1.0;
		}

		for (ie = 0; ie < 3; ie++) //  for each dof in edge
		{
			for (ii = 0; ii < elementDOFdg->dop - 1; ii++)
			{
				curnode = elementDOFdg->val[k][3 + ie*(elementDOFdg->dop - 1) + ii];
				P->val[P->IA[curnode]] = ((double)elementDOFdg->dop - 1 - ii) / (double)elementDOFdg->dop;
				P->val[P->IA[curnode] + 1] = (1.0 + ii) / (double)elementDOFdg->dop;
				/*				P->val[P->IA[curnode] + ie] = 0;
				P->val[P->IA[curnode] + (ie + 1) % 3] = ((double)elementDOFdg->dop - 1 - ii) / (double)elementDOFdg->dop;
				P->val[P->IA[curnode] + (ie + 2) % 3] = (1.0 + ii) / (double)elementDOFdg->dop;*/
			}
		}

		if (elementDOFdg->dop > 2) //  for each dof in element
		{
			if (elementDOFdg->dop == 3)
			{
				curnode = elementDOFdg->val[k][3 * elementDOFdg->dop];
				for (j = 0; j < 3; j++)
					P->val[P->IA[curnode] + j] = 1.0 / 3.0;
			}
			else if (elementDOFdg->dop == 4)
			{
				for (i = 0; i < 3; i++)
				{
					curnode = elementDOFdg->val[k][3 * elementDOFdg->dop + i];
					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode] + j] = 1.0 / 4.0;
					P->val[P->IA[curnode] + i] = 2.0 / 4.0;
				}
			}
			else if (elementDOFdg->dop == 5)
			{
				for (i = 0; i < 3; i++)
				{
					curnode = elementDOFdg->val[k][3 * elementDOFdg->dop + i];
					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode] + j] = 1.0 / 5.0;
					P->val[P->IA[curnode] + i] = 3.0 / 5.0;
				}

				for (ie = 0; ie < 3; ie++)
				{
					curnode = elementDOFdg->val[k][3 * elementDOFdg->dop + 3 + ie];
					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode] + j] = 2.0 / 5.0;
					P->val[P->IA[curnode] + ie] = 1.0 / 5.0;
				}
			}
		} // elementDOFdg->dop > 2
	} // k
}

 /**
 * \fn void interpP1toP2_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFp2, EDGE *edges)
 * \brief the interpolation matrix from the 1st order Lagrange element to 2nd order Lagrange element in 2d
 * \param *P pointer to the vector-version interpolation matrix
 * \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
 * \param *elementDOFcr pointer to the relation between elements and degrees of freedom of the Crouzeix–Raviart element
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
 the fourth column stores -1 if the edge is on boundary
 */
void interpP1toP2_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFp2, EDGE *edges)
{
	int i, j, ie, ii, k;
	//	int curnode[2];

	if (elementDOFp1->dop != 1 || elementDOFp2->dop != 2)
	{
		P = NULL;
		return;
	}

	int nn = elementDOFp1->dof;
	int ne = edges->row;

	P->row = elementDOFp2->dof;
	P->col = elementDOFp1->dof;
	P->IA = (int*)malloc((P->row + 1) * sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	P->IA[0] = 0;
	for (i = 0; i < nn; i++)
		P->IA[i + 1] = 1;
	for (i = nn; i < nn + ne; i++)
		P->IA[i + 1] = 2;
	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];
	P->nnz = P->IA[P->row];

	// step 2P&3P: Find the structure JA and the actual entries val of the interpolation matrix P
	P->JA = (int*)malloc(P->nnz * sizeof(int));
	P->val = (double*)malloc(P->nnz * sizeof(double));
	for (i = 0; i < nn; i++)
	{
		P->JA[P->IA[i]] = i;
		P->val[P->IA[i]] = 1;
	}
	for (i = nn; i < nn + ne; i++)
	{
		for (j = 0; j < 2; j++)
		{
			P->JA[P->IA[i] + j] = edges->val[i - nn][j];
			P->val[P->IA[i] + j] = 0.5;
		}
	}
}

/**
* \fn void interpVecP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg)
* \brief the vector-version interpolation matrix from the 1st order Lagrange element to piecewise kth order polynomial in 2d
* \param *P pointer to the vector-version interpolation matrix
* \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
* \param *elementDOFdg pointer to the relation between elements and degrees of freedom of the piecewise kth order polynomial
*/
void interpVecP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg)
{
	int i, j, ie, ii, k;
	int curnode[2];
	int dop = elementDOFdg->dop;

	if (elementDOFp1->dop != 1)
	{
		P = NULL;
		return;
	}

	P->row = elementDOFdg->dof * 2;
	P->col = elementDOFp1->dof * 2;
	P->IA = (int*)malloc((P->row + 1) * sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	P->IA[0] = 0;
	for (k = 0; k < elementDOFdg->row; k++)
	{
		if (dop == 0)
		{
			for (i = 0; i < elementDOFdg->col; i++) //  for each node
			{
				curnode[0] = elementDOFdg->val[k][i];
				curnode[1] = curnode[0] + elementDOFdg->dof;
				P->IA[curnode[0] + 1] = elementDOFp1->col;
				P->IA[curnode[1] + 1] = elementDOFp1->col;
			}
			continue;
		}

		for (i = 0; i < 3; i++) //  for each vertex
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			P->IA[curnode[0] + 1] = 1;
			P->IA[curnode[1] + 1] = 1;
		}
		for (i = 3; i < 3 * dop; i++) //  for each node in edge
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			P->IA[curnode[0] + 1] = 2;
			P->IA[curnode[1] + 1] = 2;
		}
		for (i = 3 * dop; i < elementDOFdg->col; i++) //  for each node in element
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			P->IA[curnode[0] + 1] = elementDOFp1->col;
			P->IA[curnode[1] + 1] = elementDOFp1->col;
		}
	}
	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];
	P->nnz = P->IA[P->row];

	// step 2P: Find the structure JA of the interpolation matrix P
	P->JA = (int*)malloc(P->nnz * sizeof(int));
	for (k = 0; k < elementDOFdg->row; k++)
	{
		if (dop == 0)
		{
			for (i = 0; i < elementDOFdg->col; i++) //  for each node
			{
				curnode[0] = elementDOFdg->val[k][i];
				curnode[1] = curnode[0] + elementDOFdg->dof;
				for (j = 0; j < elementDOFp1->col; j++)
				{
					P->JA[P->IA[curnode[0]] + j] = elementDOFp1->val[k][j];
					P->JA[P->IA[curnode[1]] + j] = elementDOFp1->val[k][j] + elementDOFp1->dof;
				}
			}
			continue;
		}

		for (i = 0; i < 3; i++) //  for each vertex
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			P->JA[P->IA[curnode[0]]] = elementDOFp1->val[k][i];
			P->JA[P->IA[curnode[1]]] = elementDOFp1->val[k][i] + elementDOFp1->dof;
		}
		for (i = 3; i < 3 * dop; i++) //  for each node in edge
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			ie = (i - 3) / (dop - 1);
			for (j = 0; j < 2; j++)
			{
				P->JA[P->IA[curnode[0]] + j] = elementDOFp1->val[k][(ie + j) % 3];
				P->JA[P->IA[curnode[1]] + j] = elementDOFp1->val[k][(ie + j) % 3] + elementDOFp1->dof;
			}
		}
		for (i = 3 * dop; i < elementDOFdg->col; i++) //  for each node in element
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			for (j = 0; j < elementDOFp1->col; j++)
			{
				P->JA[P->IA[curnode[0]] + j] = elementDOFp1->val[k][j];
				P->JA[P->IA[curnode[1]] + j] = elementDOFp1->val[k][j] + elementDOFp1->dof;
			}
		}
	}

	// step 3P: Loop element by element and compute the actual entries storing them in P
	P->val = (double*)malloc(P->nnz * sizeof(double));
	for (k = 0; k < elementDOFdg->row; k++)
	{
		if (elementDOFdg->dop == 0)
		{
			curnode[0] = elementDOFdg->val[k][0];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			for (j = 0; j < 3; j++)
			{
				P->val[P->IA[curnode[0]] + j] = 1.0 / 3.0;
				P->val[P->IA[curnode[1]] + j] = 1.0 / 3.0;
			}
			continue;
		}

		for (i = 0; i < 3; i++) //  for each vertex
		{
			curnode[0] = elementDOFdg->val[k][i];
			curnode[1] = curnode[0] + elementDOFdg->dof;
			P->val[P->IA[curnode[0]]] = 1.0;
			P->val[P->IA[curnode[1]]] = 1.0;
		}

		for (ie = 0; ie < 3; ie++) //  for each dof in edge
		{
			for (ii = 0; ii < elementDOFdg->dop - 1; ii++)
			{
				curnode[0] = elementDOFdg->val[k][3 + ie*(elementDOFdg->dop - 1) + ii];
				curnode[1] = curnode[0] + elementDOFdg->dof;

				P->val[P->IA[curnode[0]]] = ((double)elementDOFdg->dop - 1 - ii) / (double)elementDOFdg->dop;
				P->val[P->IA[curnode[0]] + 1] = (1.0 + ii) / (double)elementDOFdg->dop;
				/*				P->val[P->IA[curnode[0]] + ie] = 0;
				P->val[P->IA[curnode[0]] + (ie + 1) % 3] = ((double)elementDOFdg->dop - 1 - ii) / (double)elementDOFdg->dop;
				P->val[P->IA[curnode[0]] + (ie + 2) % 3] = (1.0 + ii) / (double)elementDOFdg->dop;*/

				for (j = 0; j < 3; j++)
					P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
			}
		}

		if (elementDOFdg->dop > 2) //  for each dof in element
		{
			if (elementDOFdg->dop == 3)
			{
				curnode[0] = elementDOFdg->val[k][3 * elementDOFdg->dop];
				curnode[1] = curnode[0] + elementDOFdg->dof;

				for (j = 0; j < 3; j++)
				{
					P->val[P->IA[curnode[0]] + j] = 1.0 / 3.0;
					P->val[P->IA[curnode[1]] + j] = 1.0 / 3.0;
				}
			}
			else if (elementDOFdg->dop == 4)
			{
				for (i = 0; i < 3; i++)
				{
					curnode[0] = elementDOFdg->val[k][3 * elementDOFdg->dop + i];
					curnode[1] = curnode[0] + elementDOFdg->dof;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[0]] + j] = 1.0 / 4.0;

					P->val[P->IA[curnode[0]] + i] = 2.0 / 4.0;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
				}
			}
			else if (elementDOFdg->dop == 5)
			{
				for (i = 0; i < 3; i++)
				{
					curnode[0] = elementDOFdg->val[k][3 * elementDOFdg->dop + i];
					curnode[1] = curnode[0] + elementDOFdg->dof;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[0]] + j] = 1.0 / 5.0;

					P->val[P->IA[curnode[0]] + i] = 3.0 / 5.0;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
				}

				for (ie = 0; ie < 3; ie++)
				{
					curnode[0] = elementDOFdg->val[k][3 * elementDOFdg->dop + 3 + ie];
					curnode[1] = curnode[0] + elementDOFdg->dof;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[0]] + j] = 2.0 / 5.0;

					P->val[P->IA[curnode[0]] + ie] = 1.0 / 5.0;

					for (j = 0; j < 3; j++)
						P->val[P->IA[curnode[1]] + j] = P->val[P->IA[curnode[0]] + j];
				}
			}
		} // elementDOFdg->dop > 2

	} // k
}

/**
* \fn void interpVecP1toNcP1_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges)
* \brief the vector-version interpolation matrix from the 1st order Lagrange element to nonconforming P1 element in 2d
* \param *P pointer to the vector-version interpolation matrix
* \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
* \param *elementDOFcr pointer to the relation between elements and degrees of freedom of the nonconforming P1 element
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
*/
void interpVecP1toNcP1_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges)
{
	int i, j, ie, ii, k;
	int curnode[2];

	if (elementDOFp1->dop != 1)
	{
		P = NULL;
		return;
	}

	P->row = elementDOFcr->dof * 2;
	P->col = elementDOFp1->dof * 2;
	P->IA = (int*)malloc((P->row + 1) * sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	P->IA[0] = 0;
	for (i = 0; i<P->row; i++)
		P->IA[i + 1] = 2;

	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];

	P->nnz = P->IA[P->row];

	// step 2P: Find the structure JA of the interpolation matrix P
	P->JA = (int*)malloc(P->nnz * sizeof(int));
	for (i = 0; i < elementDOFcr->dof; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFcr->dof;
		for (j = 0; j < 2; j++)
		{
			P->JA[P->IA[curnode[0]] + j] = edges->val[i][j];
			P->JA[P->IA[curnode[1]] + j] = edges->val[i][j] + elementDOFp1->dof;
		}
	}

	// step 3P: Loop element by element and compute the actual entries storing them in P
	P->val = (double*)malloc(P->nnz * sizeof(double));
	for (i = 0; i<P->nnz; i++)
		P->val[i] = 0.5;
}

/**
* \fn void interpVecP1toMINI_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFmini)
* \brief the vector-version interpolation matrix from the 1st order Lagrange element to MINI element for Stokes equation in 2d
* \param *P pointer to the vector-version interpolation matrix
* \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
* \param *elementDOFmini pointer to the relation between elements and degrees of freedom of the MINI element
*/
void interpVecP1toMINI_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFmini)
{
	int i, j, ie, ii, k;
	int curnode[2];

	if (elementDOFp1->dop != 1)
	{
		P = NULL;
		return;
	}

	int nn = elementDOFp1->dof;

	P->row = elementDOFmini->dof * 2;
	P->col = elementDOFp1->dof * 2;
	P->IA = (int*)malloc((P->row + 1) * sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	P->IA[0] = 0;
	for (i = 0; i < nn; i++)
	{
		P->IA[i + 1] = 1;
		P->IA[i + 1 + elementDOFmini->dof] = 1;
	}
	for (i = nn; i < elementDOFmini->dof; i++)
	{
		P->IA[i + 1] = 3;
		P->IA[i + 1 + elementDOFmini->dof] = 3;
	}

	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];

	P->nnz = P->IA[P->row];

	// step 2P&3P: Find the structure JA and the actual entries val of the interpolation matrix P
	P->JA = (int*)malloc(P->nnz * sizeof(int));
	P->val = (double*)malloc(P->nnz * sizeof(double));
	for (i = 0; i < nn; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFmini->dof;

		P->JA[P->IA[curnode[0]]] = i;
		P->JA[P->IA[curnode[1]]] = i + nn;
		P->val[P->IA[curnode[0]]] = 1;
		P->val[P->IA[curnode[1]]] = 1;
	}
	for (i = nn; i < elementDOFmini->dof; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFmini->dof;
		for (j = 0; j < 3; j++)
		{
			P->JA[P->IA[curnode[0]] + j] = elementDOFp1->val[i - nn][j];
			P->JA[P->IA[curnode[1]] + j] = elementDOFp1->val[i - nn][j] + nn;
			P->val[P->IA[curnode[0]] + j] = 1.0 / 3.0;
			P->val[P->IA[curnode[1]] + j] = 1.0 / 3.0;
		}
	}
}

/**
* \fn void interpVecP1toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges)
* \brief the vector-version interpolation matrix from the 1st order Lagrange element to Crouzeix–Raviart element for Stokes equation in 2d
* \param *P pointer to the vector-version interpolation matrix
* \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
* \param *elementDOFcr pointer to the relation between elements and degrees of freedom of the Crouzeix–Raviart element
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
*/
void interpVecP1toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges)
{
	int i, j, ie, ii, k;
	int curnode[2];

	if (elementDOFp1->dop != 1)
	{
		P = NULL;
		return;
	}

	int nn = elementDOFp1->dof;
	int ne = edges->row;

	P->row = elementDOFcr->dof * 2;
	P->col = elementDOFp1->dof * 2;
	P->IA = (int*)malloc((P->row + 1) * sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	P->IA[0] = 0;
	for (i = 0; i < nn; i++)
	{
		P->IA[i + 1] = 1;
		P->IA[i + 1 + elementDOFcr->dof] = 1;
	}
	for (i = nn; i < nn + ne; i++)
	{
		P->IA[i + 1] = 2;
		P->IA[i + 1 + elementDOFcr->dof] = 2;
	}
	for (i = nn + ne; i < elementDOFcr->dof; i++)
	{
		P->IA[i + 1] = 3;
		P->IA[i + 1 + elementDOFcr->dof] = 3;
	}

	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];

	P->nnz = P->IA[P->row];

	// step 2P&3P: Find the structure JA and the actual entries val of the interpolation matrix P
	P->JA = (int*)malloc(P->nnz * sizeof(int));
	P->val = (double*)malloc(P->nnz * sizeof(double));
	for (i = 0; i < nn; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFcr->dof;

		P->JA[P->IA[curnode[0]]] = i;
		P->JA[P->IA[curnode[1]]] = i + nn;
		P->val[P->IA[curnode[0]]] = 1;
		P->val[P->IA[curnode[1]]] = 1;
	}
	for (i = nn; i < nn + ne; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFcr->dof;
		for (j = 0; j < 2; j++)
		{
			P->JA[P->IA[curnode[0]] + j] = edges->val[i - nn][j];
			P->JA[P->IA[curnode[1]] + j] = edges->val[i - nn][j] + nn;
			P->val[P->IA[curnode[0]] + j] = 0.5;
			P->val[P->IA[curnode[1]] + j] = 0.5;
		}
	}
	for (i = nn + ne; i < elementDOFcr->dof; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFcr->dof;
		for (j = 0; j < 3; j++)
		{
			P->JA[P->IA[curnode[0]] + j] = elementDOFp1->val[i - nn - ne][j];
			P->JA[P->IA[curnode[1]] + j] = elementDOFp1->val[i - nn - ne][j] + nn;
			P->val[P->IA[curnode[0]] + j] = 1.0 / 3.0;
			P->val[P->IA[curnode[1]] + j] = 1.0 / 3.0;

		}
	}
}

/**
* \fn void interpVecP2toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp2, ELEMENT_DOF *elementDOFcr)
* \brief the vector-version interpolation matrix from the 2st order Lagrange element to Crouzeix–Raviart element for Stokes equation in 2d
* \param *P pointer to the vector-version interpolation matrix
* \param *elementDOFp1 pointer to the relation between elements and degrees of freedom of the 1st order Lagrange element
* \param *elementDOFcr pointer to the relation between elements and degrees of freedom of the Crouzeix–Raviart element
*/
void interpVecP2toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp2, ELEMENT_DOF *elementDOFcr)
{
	int i, j, ie, ii, k;
	int curnode[2];

	if (elementDOFp2->dop != 2)
	{
		P = NULL;
		return;
	}

	//	int nn = elementDOFp2->dof;
	//	int ne = edges->row;

	P->row = elementDOFcr->dof * 2;
	P->col = elementDOFp2->dof * 2;
	P->IA = (int*)malloc((P->row + 1) * sizeof(int));
	P->JA = NULL;
	P->val = NULL;

	// step 1P: Find first the structure IA of the interpolation matrix P
	P->IA[0] = 0;
	for (i = 0; i < elementDOFp2->dof; i++)
	{
		P->IA[i + 1] = 1;
		P->IA[i + 1 + elementDOFcr->dof] = 1;
	}
	for (i = elementDOFp2->dof; i < elementDOFcr->dof; i++)
	{
		P->IA[i + 1] = 6;
		P->IA[i + 1 + elementDOFcr->dof] = 6;
	}

	for (i = 0; i<P->row; i++)
		P->IA[i + 1] += P->IA[i];

	P->nnz = P->IA[P->row];

	// step 2P&3P: Find the structure JA and the actual entries val of the interpolation matrix P
	P->JA = (int*)malloc(P->nnz * sizeof(int));
	P->val = (double*)malloc(P->nnz * sizeof(double));
	for (i = 0; i < elementDOFp2->dof; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFcr->dof;

		P->JA[P->IA[curnode[0]]] = i;
		P->JA[P->IA[curnode[1]]] = i + elementDOFp2->dof;
		P->val[P->IA[curnode[0]]] = 1;
		P->val[P->IA[curnode[1]]] = 1;
	}
	for (i = elementDOFp2->dof; i < elementDOFcr->dof; i++)
	{
		curnode[0] = i;
		curnode[1] = curnode[0] + elementDOFcr->dof;
		for (j = 0; j < 3; j++)
		{
			P->JA[P->IA[curnode[0]] + j] = elementDOFp2->val[i - elementDOFp2->dof][j];
			P->JA[P->IA[curnode[1]] + j] = elementDOFp2->val[i - elementDOFp2->dof][j] + elementDOFp2->dof;
			P->val[P->IA[curnode[0]] + j] = -1.0 / 9.0;
			P->val[P->IA[curnode[1]] + j] = -1.0 / 9.0;
		}
		for (j = 3; j < 6; j++)
		{
			P->JA[P->IA[curnode[0]] + j] = elementDOFp2->val[i - elementDOFp2->dof][j];
			P->JA[P->IA[curnode[1]] + j] = elementDOFp2->val[i - elementDOFp2->dof][j] + elementDOFp2->dof;
			P->val[P->IA[curnode[0]] + j] = 4.0 / 9.0;
			P->val[P->IA[curnode[1]] + j] = 4.0 / 9.0;
		}
	}
}

/**
 * \fn void interpolation(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param)
 * \brief Generate interpolation P 
 * \param *A pointer to the stiffness matrix
 * \param *vertices pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param *P pointer to the dCSRmat matrix of resulted interpolation
 * \param *param pointer to AMG parameters
 * \return void
 */
void interpolation(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param)
{
	int i, j, index;
	int interp_type=param->interpolation_type;
	
	// change the fine index to coarse index
	if(interp_type == 2 || interp_type == 3)
	{
		int *CoarseIndex=(int*)calloc(vertices->row, sizeof(int));
		index=0;
		for(i=0;i<vertices->row;i++)
		{
			if(vertices->val[i]==1)
			{
				CoarseIndex[i]=index;
				index++;
			}
		}
		
		for(i=0;i<P->nnz;i++)
		{
			j=P->JA[i];
			P->JA[i]=CoarseIndex[j];
		}
		free(CoarseIndex);
	}
	
	/*-- Standard interpolation operator --*/
	if ( interp_type == 1 ) 
		interpolationRS(A, vertices, P, param);
	
	/*-- Energy minimization interpolation operator in Fortran --*/		
/*	if ( interp_type == 2 ) {
		int *IBmarker = (int*)calloc(A->row, sizeof(int));		
		int *ia_xz = (int*)calloc(A->row+1, sizeof(int));
		int *ja_xz = (int*)calloc(A->IA[A->row], sizeof(int));     
		double *a_xz = (double*)calloc(A->IA[A->row], sizeof(double));
		
		for(i=0;i<=A->row;i++) ia_xz[i]=A->IA[i];
		
		for(i=0;i<A->IA[A->row];i++) {
			ja_xz[i]=A->JA[i];
			a_xz[i]=A->val[i];
		}
		
		int nf_xz = P->row; 
		int nc_xz = P->col;
		int pnz_xz = P->IA[P->row];
		int anz_xz = A->IA[A->row];
		
		int *ip_xz = (int*)calloc(nf_xz+1, sizeof(int));
		int *jp_xz = (int*)calloc(pnz_xz, sizeof(int));      
		double *p_xz = (double*)calloc(pnz_xz, sizeof(double));    
		
		for(i=0;i<=P->row;i++) ip_xz[i]=P->IA[i];
		
		for(i=0;i<P->IA[A->row];i++) {
			jp_xz[i]=P->JA[i];
			p_xz[i]=P->val[i];
		}			  
		
		int  *ipt_xz = (int*)calloc(nc_xz+1, sizeof(int));
		int  *jpt_xz = (int*)calloc(pnz_xz, sizeof(int));
		
		for(i=0; i<=nf_xz; i++) ia_xz[i] += 1;
		for(i=0; i<anz_xz; i++) ja_xz[i] += 1;
		for(i=0; i<=nf_xz; i++) ip_xz[i] += 1;
		for(i=0; i<pnz_xz; i++) jp_xz[i] += 1;
		
		get_p_xuludmil_(ia_xz,ja_xz,a_xz,&nf_xz,&nc_xz,ip_xz,jp_xz,p_xz,ipt_xz,jpt_xz,IBmarker);
		
		for(i=0; i<pnz_xz; i++) P->val[i] = p_xz[i];
		
		free(IBmarker);
		free(ia_xz);
		free(ja_xz);
		free(a_xz);
		free(ip_xz);
		free(jp_xz);
		free(p_xz);
		free(ipt_xz);
		free(jpt_xz);
	} */
	
	/*-- Energy minimization interpolation operator in C --*/		
	if ( interp_type == 3 ) {
		getiteval(A, P);
	}
}	

/**
 * \fn void interpolationRS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param)
 * \brief Direct interpolation 
 * \param *A pointer to the stiffness matrix
 * \param *vertices pointer to the indicator of CF split node is on fine(current level) or coarse grid (fine: 0; coarse: 1)
 * \param *Ptr pointer to the dCSRmat matrix of resulted interpolation
 * \param *param pointer to AMG parameters
 * \return void
 *
 * Refer to P479, U. Trottenberg, C. W. Oosterlee, and A. Sch¨uller. Multigrid. 
 *	 Academic Press Inc., San Diego, CA, 2001. 
 * With contributions by A. Brandt, P. Oswald and K. St¨uben.
 */
void interpolationRS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param)
{
	double epsilon_tr = param->truncation_threshold;

	int i,j,k,l;
	int index, indexJ;
	
	/** Generate interpolation P */
	dCSRmat P;
	P.row=Ptr->row;
	P.col=Ptr->col;
	P.nnz=Ptr->nnz;
	P.val=NULL;
	P.JA=NULL;
	P.IA=(int*)calloc(P.row+1, sizeof(int));
	/** step 1: Find first the structure IA of P */
	for(i=0;i<=P.row;i++)
		P.IA[i]=Ptr->IA[i];
	
	/** step 2: Find the structure JA of P */
	P.JA=(int*)calloc(P.nnz,sizeof(int));
	for(i=0;i<P.nnz;i++)
		P.JA[i]=Ptr->JA[i];
		
	/** step 3: Fill the data A of P */
	double amN, amP, apN, apP;
	double alpha, beta, aii;
	int countPplus;
	int diagindex;
	index=0;
	P.val=(double*)calloc(P.nnz,sizeof(double));
	for(i=0;i<A->row;i++)
	{
		for(diagindex=A->IA[i];diagindex<A->IA[i+1];diagindex++)
		{
			if(A->JA[diagindex]==i)
			{
				aii=A->val[diagindex];
				break;
			}
		}
		
		if(vertices->val[i]==0)  // if node i is on fine grid 
		{
			amN=0, amP=0, apN=0, apP=0,  countPplus=0;
			for(j=A->IA[i];j<A->IA[i+1];j++)
			{
				if(j==diagindex)
					continue;
				
				for(k=Ptr->IA[i];k<Ptr->IA[i+1];k++)
				{
					if(Ptr->JA[k]==A->JA[j])
						break;
				}
				
				if(A->val[j]>0)
				{
					apN+=A->val[j];
					if(k<Ptr->IA[i+1])
					{
						apP+=A->val[j];
						countPplus++;
					}
				}
				else
				{
					amN+=A->val[j];
					if(k<Ptr->IA[i+1])
					{
						amP+=A->val[j];
					}
				}
			} // j
			
			alpha=amN/amP;
			if(countPplus>0)
			{
				beta=apN/apP;
			}
			else
			{
				beta=0;
				aii+=apN;
			}
			
			for(j=P.IA[i];j<P.IA[i+1];j++)
			{
				k=P.JA[j];
				for(l=A->IA[i];l<A->IA[i+1];l++)
				{
					if(A->JA[l]==k)
					{
						break;
					}
				}
				
				if(A->val[l]>0)
				{
					P.val[j]=-beta*A->val[l]/aii;
				}
				else
				{
					P.val[j]=-alpha*A->val[l]/aii;
				}
			}
		}
		else if(vertices->val[i]==2)  // if node i is a special fine node
		{
		}
		else // if node i is on coarse grid 
		{
			P.val[P.IA[i]]=1;
		}
	}	
	
	free(Ptr->IA);
	free(Ptr->JA);
	free(Ptr->val);	
	
	int *CoarseIndex=(int*)calloc(vertices->row, sizeof(int));
	index=0;
	for(i=0;i<vertices->row;i++)
	{
		if(vertices->val[i]==1)
		{
			CoarseIndex[i]=index;
			index++;
		}
	}
	
	for(i=0;i<P.IA[P.row];i++)
	{
		j=P.JA[i];
		P.JA[i]=CoarseIndex[j];
	}
	free(CoarseIndex);
	
	/** Truncation of interpolation */
	double mMin, pMax;
	double mSum, pSum;
	double mTruncedSum, pTruncedSum;
	int mTruncCount, pTruncCount;
	int num_lost=0;
	Ptr->val=(double*)calloc(P.IA[Ptr->row],sizeof(double));
	Ptr->JA=(int*)calloc(P.IA[Ptr->row],sizeof(int));
	Ptr->IA=(int*)calloc(Ptr->row+1, sizeof(int));
	int index1=0, index2=0;
	for(i=0;i<P.row;i++)
	{
		mMin=0;
		pMax=0;
		mSum=0;
		pSum=0;
		mTruncedSum=0;
		pTruncedSum=0;
		mTruncCount=0;
		pTruncCount=0;
		
		Ptr->IA[i]-=num_lost;
		
		for(j=P.IA[i];j<P.IA[i+1];j++)
		{
			if(P.val[j]<0)
			{
				mSum+=P.val[j];
				if(P.val[j]<mMin)
				{
					mMin=P.val[j];
				}
			}
			
			if(P.val[j]>0)
			{
				pSum+=P.val[j];
				if(P.val[j]>pMax)
				{
					pMax=P.val[j];
				}
			}
		}
		
		for(j=P.IA[i];j<P.IA[i+1];j++)
		{
			if(P.val[j]<0)
			{
				if(P.val[j]>mMin*epsilon_tr)
				{
					mTruncCount++;
				}
				else
				{
					num_lost--;
				}
			}
			
			if(P.val[j]>0)
			{
				if(P.val[j]<pMax*epsilon_tr)
				{
					pTruncCount++;
				}
				else
				{
					num_lost--;
				}
			}
		}
				
		// step 2: Find the structure JA and fill the data A of Ptr
		for(j=P.IA[i];j<P.IA[i+1];j++)
		{
			if(P.val[j]<0)
			{
				if(!(P.val[j]>mMin*epsilon_tr))
				{
					Ptr->JA[index1]=P.JA[j];
					mTruncedSum+=P.val[j];
					index1++;
				}
			}
			
			if(P.val[j]>0)
			{
				if(!(P.val[j]<pMax*epsilon_tr))
				{
					Ptr->JA[index1]=P.JA[j];
					pTruncedSum+=P.val[j];
					index1++;
				}
			}
		}
		
		// step 3: Fill the data A of Ptr
		for(j=P.IA[i];j<P.IA[i+1];j++)
		{
			if(P.val[j]<0)
			{
				if(!(P.val[j]>mMin*epsilon_tr))
				{
					Ptr->val[index2]=P.val[j]/mTruncedSum*mSum;
					index2++;
				}
			}
			
			if(P.val[j]>0)
			{
				if(!(P.val[j]<pMax*epsilon_tr))
				{
					Ptr->val[index2]=P.val[j]/pTruncedSum*pSum;
					index2++;
				}
			}
		}
	}
	
	Ptr->IA[P.row]-=num_lost;
	Ptr->nnz=Ptr->IA[Ptr->row];
	
	Ptr->JA=(int*)realloc(Ptr->JA, Ptr->IA[Ptr->row]*sizeof(int));
	Ptr->val=(double*)realloc(Ptr->val, Ptr->IA[Ptr->row]*sizeof(double));
	
	free(P.IA);
	free(P.JA);
	free(P.val);
}



/**
	* \fn int invden(int nn, double *mat, double *invmat)
	* \brief the routine is to find the inverse of a dense matrix
	* \param nn scale of the matrix
	* \param mat the double pointer to the full matrix
	* \param invmat the double pointer to the full inverse matrix
	* \return at the end of the routine return 1
	* 
	* Note: that this routine works for symmetric matrix.
*/

int invden(int nn, double *mat, double *invmat)
{
	int *pivot;
	double *rhs,*sol;
	int indic,i,j,k;
	
	pivot=malloc(nn*sizeof(int));
	rhs=malloc(nn*sizeof(double));
	sol=malloc(nn*sizeof(double));
	
	indic=LU_Decomp(mat,pivot,nn);

	for (i=0;i<nn;i++) {
		for (j=0;j<nn;j++) rhs[j]=0.;
		rhs[i]=1.;
		LU_Solve(mat,rhs,pivot,sol,nn);
		for (j=0;j<nn;j++) invmat[i*nn+j]=sol[j];
	}
	
	free(pivot);
	free(rhs);
	free(sol);
	
	return 1;
}


/**
	* \fn int gentisquare_nomass(dCSRmat *A, int mm, int *Ii, double *ima, int *mask)
	* \brief given the row indices and col indices, to find a block submatrix and get its inverse
	* \param *A pointer to the whole matrix
	* \param mm integer of the scale of the submatrix
	* \param Ii integer to an integer array, to record the indices of row (also col)
	* \param *ima pointer to the inverse of the full submatrix, the storage is row by row
	* \param *mask working interger array
	* \return 1 if succeed
*/

int gentisquare_nomass(dCSRmat *A, int mm, int *Ii, double *ima, int *mask)
{
	int indic;
	double *ms;
	int i,j,k,l;
	
	ms=malloc(mm*mm*sizeof(double));	
		
	indic=get_block(A,mm,mm,Ii,Ii,ms,mask);
	indic=invden(mm,ms,ima);

	free(ms);
	return 1;
}

/**
	* \fn int getinonefull(int **mat, double **matval, int *lengths, int mm, int *Ii, double *ima)
	* \brief to add a small submatrix to a big matrix with respect to its row and cols in the big matrix
	* \param mat  a double pointer pointing to the structure of the matrix
	* \param matval a double pointer pointing to the values according to the structure
	* \param lengths  a 2d array, the second entry is the lengths of matval
	* \param mm the number of the rows (also the columns) 
	* \param Ii an integer pointer to the array to store the relative position of the rows and cols
	* \param ima the pointer to the full submatrix, the sequence is row by row
*/

int getinonefull(int **mat, double **matval, int *lengths, int mm, int *Ii, double *ima)
{
	int nn,nnz,tniz;
	int nr;
	int i,j,k,l;
	
	tniz=lengths[1];
	for (i=0;i<mm;i++) {
		for (j=0;j<mm;j++) {		
			mat[0][tniz+i*mm+j]=Ii[i];
			mat[1][tniz+i*mm+j]=Ii[j];
			matval[0][tniz+i*mm+j]=ima[i*mm+j];
		}
	}
	lengths[1]=tniz+mm*mm;
	
	return 1;
}


/**
	* \fn int orderone(int **mat, double **matval, int *lengths)
	* \brief Order a cluster of entries in a sequence
	* \param **mat a double pointer to the relative position of the entries
	* \param **matval a double pointer to the values corresponding to the position
	* \param lengths int array, to record the number of rows, number of cols and number of nonzerow
	* \return return 1 at the end of the routine
*/

int orderone(int **mat, double **matval, int *lengths)
//	lengths[0] for the number of rows
//	lengths[1] for the number of cols
//	lengths[2] for the number of nonzeros
{
	int *rows[2],*cols[2],nns[2],tnizs[2];
	double *vals[2];
	int tniz;
	int i,j,k,l;
	
	nns[0]=lengths[0];
	nns[1]=lengths[1];
	tnizs[0]=lengths[2];
	tniz=lengths[2];
	
	
	rows[0]=malloc(tniz*sizeof(int));
	cols[0]=malloc(tniz*sizeof(int));
	vals[0]=malloc(tniz*sizeof(double));

	for (i=0;i<tniz;i++) 
	{
		rows[0][i]=mat[0][i];
		cols[0][i]=mat[1][i];
		vals[0][i]=matval[0][i];
	}
	rows[1]=malloc(tniz*sizeof(int));
	cols[1]=malloc(tniz*sizeof(int));
	vals[1]=malloc(tniz*sizeof(double));
	
	transpose(rows,cols,vals,nns,tnizs);

	//	all the nonzeros with same col are gathering together

	for (i=0;i<tniz;i++)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	transpose(rows,cols,vals,nns,tnizs);

	//	all the nozeros with same col and row are gathering togheter	
	
	for (i=0;i<tniz;i++)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	
	tniz=tnizs[0];
	for (i=0;i<tniz-1;i++) {
		if (rows[0][i]==rows[0][i+1]&&cols[0][i]==cols[0][i+1]) {
			vals[0][i+1]+=vals[0][i];
			rows[0][i]=nns[0];
			cols[0][i]=nns[1];
		}
	}
	nns[0]=nns[0]+1;
	nns[1]=nns[1]+1;
	
	transpose(rows,cols,vals,nns,tnizs);
	

	for (i=0;i<tniz;i++)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	
	transpose(rows,cols,vals,nns,tnizs);

	for (i=0;i<tniz;i++)
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	tnizs[1]=nns[0];
	nns[0]=nns[1];
	nns[1]=tnizs[1];
	tnizs[1]=tnizs[0];
	
	tniz=0;
	for (i=0;i<tnizs[0];i++)
		if (rows[0][i]<nns[0]-1) tniz++;
	
	for (i=0;i<tniz;i++)
	{
		mat[0][i]=rows[0][i];
		mat[1][i]=cols[0][i];
		matval[0][i]=vals[0][i];
	}
	nns[0]=nns[0]-1;
	nns[1]=nns[1]-1;
	lengths[0]=nns[0];
	lengths[1]=nns[1];
	lengths[2]=tniz;
	free(rows[0]);
	free(rows[1]);
	free(cols[0]);
	free(cols[1]);
	free(vals[0]);
	free(vals[1]);
}


/**
	* \fn int genintval(dCSRmat *A, int **itmat, double **itmatval, int ittniz, int nf, int nc) // nf=number fine, nc= n coarse	
	* \brief given the structure of the interpolation, to get the evaluation of the interpolation
	* \param *A  pointer to the dCSRmat matrix
	* \param **itmat  a double integer pointer pointing to the structure of the interpolation
	* \param **itmatval a double double pointer to the evaluation of the interpolation
	* \param ittniz int, the length of interpolation
	* \param nf int, the number of fine-level nodes
	* \param nc int, the number of coarse-level nodes
	* \return return 1 at the ending of the routine
*/
int genintval(dCSRmat *A, int **itmat, double **itmatval, int ittniz, int *isol, int numiso, int nf, int nc) // nf=number fine, nc= n coarse	
//	suppose that the structure of the interpolation is known
//	it is recorded in itmat.
//	it is a N*m matrix, N>m
//	we record its row index and col index, note that the same col indices gather together
//	the itma and itmatval have a special data structure
//	to be exact, the same columns gather together
//	itmat[0] record the column number, and itmat[1] record the row number
{
	int *Ii;
	double *ima;
	double **imas;
	int *mask;
	
	int **mat;
	double **matval;
	int lengths[3];
	dCSRmat T;
	int tniz;
	dvector sol, rhs;
	double *pex;
	
	int mm;
	int sum;
	int i,j,k,l;
	
	int *iz,*izs,*izt,*izts;
	FILE *fp;
	
	mask=malloc(nf*sizeof(int));
	iz=malloc(nc*sizeof(int));
	izs=malloc(nc*sizeof(int));
	izt=malloc(nf*sizeof(int));
	izts=malloc(nf*sizeof(int));
	
	for (i=0;i<nf;i++) mask[i]=-1;
	
	for (i=0;i<nc;i++) iz[i]=0;
	
	for (i=0;i<ittniz;i++) iz[itmat[0][i]]++;
	
	izs[0]=0;
	for (i=1;i<nc;i++) izs[i]=izs[i-1]+iz[i-1];
	
	for (sum=i=0;i<nc;i++) sum+=iz[i]*iz[i];
	
	imas=malloc(nc*sizeof(double *));
	for (i=0;i<nc;i++) imas[i]=malloc(iz[i]*iz[i]*sizeof(double));
	mat=malloc(2*sizeof(int *));
	mat[0]=malloc((sum+numiso)*sizeof(int));
	mat[1]=malloc((sum+numiso)*sizeof(int));
	matval=malloc(1*sizeof(double *));
	matval[0]=malloc((sum+numiso)*sizeof(double));
	
	lengths[1]=0;

	for (i=0;i<nc;i++) {	
		mm=iz[i]; 
		Ii=malloc(mm*sizeof(int));
		for (j=0;j<mm;j++) Ii[j]=itmat[1][izs[i]+j];
		ima=malloc(mm*mm*sizeof(double));
		gentisquare_nomass(A,mm,Ii,ima,mask);
		getinonefull(mat,matval,lengths,mm,Ii,ima);

		for (j=0;j<mm*mm;j++) imas[i][j]=ima[j];
		free(ima);
		free(Ii);
	}

	for (i=0;i<numiso;i++) {
		mat[0][sum+i]=isol[i];
		mat[1][sum+i]=isol[i];
		matval[0][sum+i]=1.0;
	}
	
	lengths[0]=nf;
	lengths[2]=lengths[1]+numiso;
	lengths[1]=nf;
	orderone(mat,matval,lengths);
	tniz=lengths[2];
	sol.row=nf;
	sol.val=(double*)calloc(nf,sizeof(double));

	for (i=0;i<nf;i++) izt[i]=0;
	
	for (i=0;i<tniz;i++) izt[mat[0][i]]++;

	T.row=nf;
	T.col=nf;
	T.nnz=tniz;
	T.IA=(int*)calloc((nf+1),sizeof(int));
	T.IA[0]=0;
	for (i=1;i<nf+1;i++) T.IA[i]=T.IA[i-1]+izt[i-1];

	T.JA=(int*)calloc(tniz,sizeof(int));
	for (j=0;j<tniz;j++) T.JA[j]=mat[1][j];
	T.val=(double*)calloc(tniz,sizeof(double));
	for (j=0;j<tniz;j++) T.val[j]=matval[0][j];
	
	rhs.val=(double*)calloc(nf,sizeof(double));
	for (i=0;i<nf;i++) rhs.val[i]=1.0;
	rhs.row=nf;
	
	diag_PCG(&T,&rhs,&sol,100,1e-3,0);

	for (i=0;i<nc;i++)
	{
		mm=iz[i];
		ima=malloc(mm*mm*sizeof(double));
		pex=malloc(mm*sizeof(double));
		Ii=malloc(mm*sizeof(int));
		
		for (j=0;j<mm;j++) Ii[j]=itmat[1][izs[i]+j];
		
		for (j=0;j<mm*mm;j++) ima[j]=imas[i][j];
		
		for (k=0;k<mm;k++)
		{
			for(pex[k]=j=0;j<mm;j++) pex[k]+=ima[k*mm+j]*sol.val[Ii[j]];
		}
		for (j=0;j<mm;j++) itmatval[0][izs[i]+j]=pex[j];

		free(ima);
		free(pex);
		free(Ii);		
	}
	
	free(mask);
	free(iz);
	free(izs);
	free(izt);
	free(izts);
	free(mat[0]);
	free(mat[1]);
	free(matval[0]);
}


/**
	* \fn int getiteval(dCSRmat *A, dCSRmat *it)
	* \brief given a coarsening (in the form of an interpolation operator), inherit the structure, get new evaluation
	* \param *A ponter to the dCSRmat matrix
	* \param *it pointer to the interpolation matrix
*/
int getiteval(dCSRmat *A, dCSRmat *it)
{
	int nf,nc,ittniz;
	int *itmat[2];
	double **itmatval;
	int *rows[2],*cols[2];
	double *vals[2];
	int nns[2],tnizs[2];
	int i,j,k,l;
	int indic;
	int *isol;
	int numiso;
	
	nf=A->row;
	nc=it->col;
	ittniz=it->IA[nf]; 
		
	itmat[0]=malloc(ittniz*sizeof(int));
	itmat[1]=malloc(ittniz*sizeof(int));
	itmatval=malloc(1*sizeof(double *));
	itmatval[0]=malloc(ittniz*sizeof(double));
	
	isol=malloc(nf*sizeof(int));
	numiso=0;
	for (i=0;i<nf;i++) {
		if (it->IA[i]==it->IA[i+1]) {
			isol[numiso]=i;
			numiso++;
		}
	}

	for (i=0;i<nf;i++) {
		for (j=it->IA[i];j<it->IA[i+1];j++) itmat[0][j]=i;
	}
	
	for (j=0;j<ittniz;j++) itmat[1][j]=it->JA[j];
	
	for (j=0;j<ittniz;j++) itmatval[0][j]=it->val[j];
	
	rows[0]=malloc(ittniz*sizeof(int));
	cols[0]=malloc(ittniz*sizeof(int));
	vals[0]=malloc(ittniz*sizeof(double));
	
	for (i=0;i<ittniz;i++)
	{
		rows[0][i]=itmat[0][i];
		cols[0][i]=itmat[1][i];
		vals[0][i]=itmat[0][i];
	}

	nns[0]=nf;
	nns[1]=nc;
	tnizs[0]=ittniz;
	
	rows[1]=malloc(ittniz*sizeof(int));
	cols[1]=malloc(ittniz*sizeof(int));
	vals[1]=malloc(ittniz*sizeof(double));
	transpose(rows,cols,vals,nns,tnizs);
	for (i=0;i<ittniz;i++)
	{
		itmat[0][i]=rows[1][i];
		itmat[1][i]=cols[1][i];
		itmatval[0][i]=vals[1][i];
	}
	indic=genintval(A,itmat,itmatval,ittniz,isol,numiso,nf,nc);

	for (i=0;i<ittniz;i++)
	{
		rows[0][i]=itmat[0][i];
		cols[0][i]=itmat[1][i];
		vals[0][i]=itmatval[0][i];
	}
	nns[0]=nc;
	nns[1]=nf;
	tnizs[0]=ittniz;
		
	transpose(rows,cols,vals,nns,tnizs);
	
	for (i=0;i<ittniz;i++) it->val[i]=vals[1][i];
	
	free(isol);
	free(itmat[0]); free(itmat[1]);
	free(itmatval[0]); free(itmatval);
	
	free(rows[0]);
	free(rows[1]);
	free(cols[0]);
	free(cols[1]);
	free(vals[0]);
	free(vals[1]);
	       
	
	return 1;
}
	
void contiso(dCSRmat *m)
{
	int nn;
	int numiso;
	int *isol;
	int i;
	
	numiso=0;
	for (i=0;i<m->row;i++)
		if (m->IA[i]==m->IA[i+1]) numiso++;
	
	printf("the number of i=%d\n",numiso);
}
