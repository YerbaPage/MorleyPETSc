/*
 *  rref.c
 *  FEM
 *
 *  Created by Xuehai Huang on 11/11/12.
 *  Copyright 2012 WZU. All rights reserved.
 *
 */

/*! \file rref.c
 *  \brief produces the reduced row echelon form of matrix using Gauss Jordan elimination
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn void AxBrref(ddenmat *A, ddenmat *B)
 * \brief produces the reduced row echelon form of A using Gauss Jordan elimination
 * \param *A pointer to ddenmat matrix
 * \param *B pointer to the reduced row echelon form
 */
void AxBrref(ddenmat *A, ddenmat *B)
{	
	int i,j, maxi;
	double max, temp;
	
	for(j=0;j<A->col;j++)
	{
		max=fabs(A->val[j][j]);
		maxi=j;
		for(i=j+1;i<A->row;i++)
		{
			if(fabs(A->val[i][j])>max)
			{
				max=fabs(A->val[i][j]);
				maxi=i;
			}
		}
		rref1(A->val, j, maxi, A->row, A->col);
		rref1(B->val, j, maxi, B->row, B->col);
		temp=A->val[j][j];
	//	rref2(A->val, j, 1.0/temp, A->row, A->col);
	//	rref2(B->val, j, 1.0/temp, B->row, B->col);		
		rref2divide(A->val, j, temp, A->row, A->col);
		rref2divide(B->val, j, temp, B->row, B->col);
		
		for(i=0;i<A->row;i++)
		{
			temp=A->val[i][j];
			rref3(A->val, j, -temp, i, A->row, A->col);
			rref3(B->val, j, -temp, i, B->row, B->col);
		}

	}// j
}


void rref1(double **A, int i, int j, int row, int col)
{
	if(i<0 || i>=row)
		return;
	if(j<0 || j>=row)
		return;
	if(i==j)
		return;

	int k;
	double *Ai;

	Ai=(double*)calloc(col, sizeof(double));
	for(k=0;k<col;k++)
	{
		Ai[k]=A[i][k];
		A[i][k]=A[j][k];
		A[j][k]=Ai[k];
	}

	free(Ai);
}

void rref2(double **A, int i, double alpha, int row, int col)
{
	if(i<0 || i>=row)
		return;
	
	int k;
	for(k=0;k<col;k++)
		A[i][k]*=alpha;
}

void rref2divide(double **A, int i, double alpha, int row, int col)
{
	if (i<0 || i >= row)
		return;

	int k;
	for (k = 0; k<col; k++)
		A[i][k] /= alpha;
}

void rref3(double **A, int i, double alpha, int j, int row, int col)
{
	if(i<0 || i>=row)
		return;
	if(j<0 || j>=row)
		return;
	if(i==j)
		return;

	int k;
	for(k=0;k<col;k++)
		A[j][k]+=A[i][k]*alpha;
}