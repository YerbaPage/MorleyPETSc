/*
 *		Matrix properties subroutines
 *
 *		Matrix-Solver Community Project
 *
 *------------------------------------------------------
 *
 *		Created by Shuo Zhang on 04/01/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------
 *
 */

/*! \file checkmat.c
 *  \brief Check matrix properties
 */

#include <stdio.h>
#include <stdlib.h>

#include "header.h"
#include "matvec.h"
#include "math.h"

/**
 * \fn int check_diagpos(dCSRmat *A)
 * \brief Check positivity of diagonal entries of a CSR sparse matrix.
 * \param *A pointer to the dCSRmat matrix
 * \return number of negative entries 
 */
int check_diagpos(dCSRmat *A)
{
	int m=A->row, num_neg, i;
	dvector diag;
	
	create_dvector(m,&diag);
  getdiag(m,A,&diag);
	
	for (num_neg=i=0;i<m;i++) 
		if (diag.val[i]<0) num_neg++;
		
	printf("check_diagpos: number of negative diagonal entries = %d\n", num_neg);
	
	free_dvector(&diag);
	return num_neg;
}

/* Private function
 * \fn void transpose(int *row[2], int *col[2], double *val[2], int *nn, int *tniz)
 * \brief Transpose a matrix. input the row column value and output row column value.
 * \param *row[2] to record the pointer of the rows of the matrix and its transpose
 * \param *col[2] to record the pointer of the columns of the matrix and its transpose
 * \param *val[2] to record the pointers to the values of the matrix and its transpose
 * \param *nn to record the number of rows and columns of the matrix
 * \param *tniz to record the number of the nonzeros in the matrix
 */
void transpose(int *row[2], int *col[2], double *val[2], int *nn, int *tniz)
{
  int i,k,m;
  int *izc,*izcaux;
  int nca;
  int *orow,*ocol,*trow,*tcol;
  double *oval,*tval;
  
  nca=nn[1];
  izc=malloc(nn[1]*sizeof(int));
  izcaux=malloc(nn[1]*sizeof(int));
	
  for (i=0;i<nca;i++) izc[i]=0;
  for (i=0;i<tniz[0];i++) izc[col[0][i]]++;
	
  izcaux[0]=0;
  for (i=1;i<nca;i++) izcaux[i]=izcaux[i-1]+izc[i-1];
  for (i=0;i<nca;i++) izc[i]=0;

  for (i=0;i<tniz[0];i++) 
  {
		m=col[0][i]; //ocol[i];
		k=row[0][i]; //orow[i];;
		row[1][izcaux[m]+izc[m]]=m;
		col[1][izcaux[m]+izc[m]]=k;
		val[1][izcaux[m]+izc[m]]=val[0][i]; //oval[i];
		izc[m]++;
  }

	free(izc);
	free(izcaux);
}

/* Private function
 * /fn int transpspm(dCSRmat *A, dCSRmat *B)
 * /prama *A pointer to dCSRmat matrix A
 * /prama *B pointer to the transpose of A
 * /return 1 if succeed, 0 if fail
 */
int transpspm(dCSRmat *A, dCSRmat *B)
{
	
	int i,j;
	int *rowp;
	int nnz,nn;
	int *rows[2],*cols[2];
	double *vals[2];
	int nns[2],tnizs[2];
	nnz=A->nnz;
	nn=A->row;
	rowp=malloc(nnz*sizeof(int));
	for (i=0;i<nn;i++)
	{
		for (j=A->IA[i];j<A->IA[i+1];j++) rowp[j]=i;
	}
	rows[0]=rowp;
	cols[0]=A->JA;
	vals[0]=A->val;
	nns[0]=nn;
	nns[1]=A->col;
	tnizs[0]=nnz;

	transpose(rows,cols,vals,nns,tnizs);
	
	B->nnz=tnizs[1];
	B->row=A->col;
	B->col=A->row;
	
	B->IA=(int*)calloc(B->row+1,sizeof(int));
	B->IA=rows[1];
	
	B->JA=(int*)calloc(B->nnz,sizeof(int));
	B->JA=cols[1];
	
	B->val=(double*)calloc(B->nnz, sizeof(double));
	B->val=vals[1];
	
	return 1;
}

/**
 * \fn int check_symm(dCSRmat *A)
 * \brief Check symmetry of a sparse matrix of CSR format.
 * \param *A pointer to the dCSRmat matrix
 * \return 1 and 2 if the structure of the matrix is not symmetric;
 * \return 0 if the structure of the matrix is symmetric,
 * 
 * Print the maximal relative difference between matrix and its transpose.
 */
int check_symm(dCSRmat *A)
{
	int i,j,k,mdi,mdj,nnz,nn;
	int *rowp,*rows[2],*cols[2];
	double *vals[2];
	int nns[2],tnizs[2];
	double maxdif,dif;
	
	nn=A->row;
	nnz=A->IA[nn];
	if (nnz!=A->nnz) printf("check_symm: nnz of the matrix is wrong!!!\n");

	rowp=malloc(nnz*sizeof(int));
	for (i=0;i<nn;i++) {
		for (j=A->IA[i];j<A->IA[i+1];j++) rowp[j]=i;
	}

	rows[0]=malloc(nnz*sizeof(int));
	cols[0]=malloc(nnz*sizeof(int));
	vals[0]=malloc(nnz*sizeof(double));
	for (i=0;i<nnz;i++)
	{
		rows[0][i]=rowp[i];
		cols[0][i]=A->JA[i];
		vals[0][i]=A->val[i];
	}
	nns[0]=nn;
	nns[1]=A->col;
	tnizs[0]=nnz;	
	
	rows[1]=malloc(nnz*sizeof(int));
	cols[1]=malloc(nnz*sizeof(int));
	vals[1]=malloc(nnz*sizeof(double));
	transpose(rows,cols,vals,nns,tnizs);
	
	for (i=0;i<nnz;i++) 
	{
		rows[0][i]=rows[1][i];
		cols[0][i]=cols[1][i];
		vals[0][i]=vals[1][i];
	}
	nns[0]=A->col;
	nns[1]=nn;
	transpose(rows,cols,vals,nns,tnizs);
	
	k=0;
	maxdif=0.;
	mdi=0;
	mdj=0;
	for (i=0;i<nnz;i++)
	{
		rows[0][i]=rows[1][i]-rows[0][i];
		if (rows[0][i]*rows[0][i]>0) {
			k=1;
			mdi=rows[1][i];
			break;
		}
		
		cols[0][i]=cols[1][i]-cols[0][i];
		if (cols[0][i]*cols[0][i]>0) {
			k=2;
			mdj=cols[1][i];
			break;
		}
		
		if (fabs(vals[0][i])>1e-15||fabs(vals[1][i])>1e-15) {
			dif=fabs(vals[1][i]-vals[0][i])/(fabs(vals[0][i])+fabs(vals[1][i]));
			if (dif>maxdif) {
				maxdif=dif;
				mdi=rows[0][i];
				mdj=cols[0][i];
			}
		}
//		if (i%100==0) printf("%d	%1.11le	%1.11le	%1.11le	%1.11le\n",i,fabs(vals[0][i]),fabs(vals[1][i]),dif,maxdif);
	}
	
	if (k==0) printf("check_symm: matrix is symmetric and max relative difference is %1.3le\n",maxdif);
	if (k==1) printf("check_symm: matrix is not symmetric, check the %d-th, %d-th and %d-th rows and cols\n",mdi-1,mdi,mdi+1);
	if (k==2) printf("check_symm: matrix is not symmetric, check the %d-th, %d-th and %d-th cols and rows\n",mdj-1,mdj,mdj+1);
	
	free(rowp);
	
	free(rows[1]);
	free(cols[1]);
	free(vals[1]);	

	free(rows[0]);
	free(cols[0]);
	free(vals[0]);
	
	return k;
}


/**
 * int check_diagdom(dCSRmat *A)
 * \brief Check whether a matrix is diagonal dominant.
 * \param *A pointer to the dCSRmat matrix
 * \return print the percentage of the rows which are diagonal dominant and not
 * the number of the rows which are diagonal dominant
 *
 * The routine chechs whether the sparse matrix is diagonal dominant on every row.
 *	 It will print out the percentage of the rows which are diagonal dominant and 
 * which are not; the routine will return the number of the rows which are diagonal 
 * dominant.
 */
int check_diagdom(dCSRmat *A)
{
	int i,j,k=0;
	double sum;
	int *rowp;
	int nnz,nn;
	
	nn=A->row;
	nnz=A->IA[nn]-A->IA[0];		
	rowp=(int *)malloc(nnz*sizeof(int));

	for (i=0;i<nn;i++) {
		for (j=A->IA[i];j<A->IA[i+1];j++) rowp[j]=i;
	}
	
	for (i=0;i<nn;i++) {
		sum=0.0;
		for (j=A->IA[i];j<A->IA[i+1];j++) {
			if (A->JA[j]==i) sum=sum+A->val[j];
			if (A->JA[j]!=i) sum=sum-fabs(A->val[j]);
		}
		if (sum<-1.e-16) k++;
	}
	
	printf("check_diagdom: percentage of the diagonal-dominant rows is %3.2lf%s\n", 
				 100.0*(double)(nn-k)/(double)nn,"%");
	
	free(rowp);
	
	return k;
}
