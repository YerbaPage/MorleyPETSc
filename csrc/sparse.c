/*
 *  sparse.c
 *  FEM
 *
 *  Created by Xuehai Huang on 10/23/08.
 *  Copyright 2008 PSU. All rights reserved.
 *
 */

/*! \file sparse.c
 *  \brief Functions for sparse matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask)
 * \brief Get a local block from a CSR sparse matrix
 * \param *A pointer to a sparse matrix
 * \param m number of rows of the local block matrix
 * \param n number of columns of the local block matrix
 * \param rows indices to local rows
 * \param cols indices to local columns
 * \param Aloc local dense matrix saved as an array
 * \param *mask working array, which should be a negative number initially
 * \return 1 if succeed, 0 if fail
 */
int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask)
{
	int i, j, k, ii, ij, jj;
	
	for ( i=0; i<m; i++ ) {
		for ( j=0; j<n; j++ ) {
			Aloc[i*n+j] = 0.0; // initialize Aloc
		}			
	}
	
	for ( j=0; j<n; j++ ) {
		if (cols[j]>=A->row) return 0;
		mask[cols[j]] = j; // initialize mask
	}		
	
	for ( i=0; i<m; i++ ) {
		ii=rows[i];
		for ( k=A->IA[ii]; k<A->IA[ii+1]; k++ ) {
			jj = A->JA[k];
			if (mask[jj]>=0) {
				ij = i*n+mask[jj];
				Aloc[ij]=A->val[k];
			}
		}
	}
	
	for ( j=0; j<n; j++ ) mask[cols[j]] = -1; // re-initialize mask
	
	return 1;
}

/**
* \fn int get_block_dden(dCSRmat *A, int m, int n, int *rows, int *cols, ddenmat *Aloc, int *mask)
* \brief Get a local block from a CSR sparse matrix
* \param *A pointer to a sparse matrix
* \param m number of rows of the local block matrix
* \param n number of columns of the local block matrix
* \param rows indices to local rows
* \param cols indices to local columns
* \param Aloc local dense matrix
* \param *mask working array, which should be a negative number initially
* \return 1 if succeed, 0 if fail
*/
int get_block_dden(dCSRmat *A, int m, int n, int *rows, int *cols, ddenmat *Aloc, int *mask)
{
	int i, j, k, ii, ij, jj;

/*	for (j = 0; j < n; j++) {
		if (cols[j] >= A->col) return 0;
	}*/

	for (j = 0; j<n; j++) {
		mask[cols[j]] = j; // initialize mask
	}

	create_dden_matrix(m, n, Aloc);

	for (i = 0; i<m; i++) {
		ii = rows[i];
		for (k = A->IA[ii]; k<A->IA[ii + 1]; k++) {
			jj = A->JA[k];
			if (mask[jj] >= 0) {
				Aloc->val[i][mask[jj]] = A->val[k];
			}
		}
	}

	for (j = 0; j<n; j++) mask[cols[j]] = -1; // re-initialize mask

	return 1;
}

/**
 * \fn iCSRmat getCSRfromIntT(idenmat *T, int nNode)
 * \brief Get CSR structure of trianglution T
 * \param *T pointer to the idenmat
 * \param nNode the number of nodes
 * \return the CSR structure of trianglution T
 */
iCSRmat getCSRfromIntT(idenmat *T, int nNode)
{
	iCSRmat element;
	element.val=NULL;
	element.JA=NULL;
	element.row=T->row;
	element.col=nNode;
	element.IA=(int*)calloc(element.row+1, sizeof(int));
	
	int i,j;
	for(i=0;i<element.row;i++)
	{
		element.IA[i+1]=element.IA[i]+T->col;
	}	
	element.JA=(int*)realloc(element.JA, (element.IA[element.row])*sizeof(int));
	for(i=0;i<T->row;i++)
	{
		for(j=0;j<T->col;j++)
		{
			element.JA[element.IA[i]+j]=T->val[i][j];
		}
	}	
	return element;
}


/**
 * \fn iCSRmat getTransposeOfA(iCSRmat *A)
 * \brief find A^T from given iCSRmat matrix A
 * \param *A pointer to the iCSRmat matrix
 * \return the transpose of iCSRmat matrix A
 */
iCSRmat getTransposeOfA(iCSRmat *A)
{
	int n,m;
	iCSRmat AT;
	n=A->row;
	m=A->col;
	AT.row=m;
	AT.col=n;
	AT.IA=(int*)calloc(m+1, sizeof(int));
	AT.JA=(int*)calloc(A->IA[n], sizeof(int));
	if(A->val==NULL)
	{
		AT.val=NULL;
	}
	else
	{
		AT.val=(int*)calloc(A->IA[n], sizeof(int));
	}
	int i,j,k,p;
	// find the number of nonzeros in the first m-1 columns of A 
	// and store these numbers in the array AT.IA in positions from 1 to m-1
	for(i=0;i<A->IA[n];i++)
	{
		j=A->JA[i];
		if(j<m-1)
		{
			AT.IA[j+2]=AT.IA[j+2]+1;
		}
	}
	
	for(i=2;i<=m;i++)
	{
		AT.IA[i]=AT.IA[i]+AT.IA[i-1];
	}
	
	int BEGIN_ROW, END_ROW;
	for(i=0;i<n;i++)
	{
		BEGIN_ROW=A->IA[i];
		END_ROW=A->IA[i+1];
		for(p=BEGIN_ROW;p<END_ROW;p++)
		{
			j=A->JA[p]+1;
			k=AT.IA[j];
			AT.JA[k]=i;
			if(A->val!=NULL)
			{
				AT.val[k]=A->val[p];
			}
			AT.IA[j]=k+1;
		}
	}
	return AT;
}


/**
 * \fn void getTransposeOfSparse(dCSRmat *A, dCSRmat *AT)
 * \brief find A^T from given dCSRmat matrix A
 * \param *A pointer to the dCSRmat matrix
 * \param *AT pointer to the transpose of dCSRmat matrix A
 * \return void
 */
void getTransposeOfSparse(dCSRmat *A, dCSRmat *AT)
{
	int n,m;
	n=A->row;
	m=A->col;
	AT->row=m;
	AT->col=n;
	AT->nnz=A->nnz;
	AT->IA=(int*)calloc(m+1, sizeof(int));
	AT->JA=(int*)calloc(AT->nnz, sizeof(int));
	if(A->val==NULL)
	{
		AT->val=NULL;
	}
	else
	{
		AT->val=(double*)calloc(AT->nnz, sizeof(double));
	}	
	int i,j,k,p;
	// find the number of nonzeros in the first m-1 columns of A 
	// and store these numbers in the array AT.IA in positions from 1 to m-1
	for(i=0;i<A->IA[n];i++)
	{
		j=A->JA[i];
		if(j<m-1)
		{
			AT->IA[j+2]++;
		}
	}
	
	for(i=2;i<=m;i++)
	{
		AT->IA[i]+=AT->IA[i-1];
	}
	
	int BEGIN_ROW, END_ROW;
	for(i=0;i<n;i++)
	{
		BEGIN_ROW=A->IA[i];
		END_ROW=A->IA[i+1];
		for(p=BEGIN_ROW;p<END_ROW;p++)
		{
			j=A->JA[p]+1;
			k=AT->IA[j];
			AT->JA[k]=i;
			if(A->val!=NULL)
			{
				AT->val[k]=A->val[p];
			}
			AT->IA[j]=k+1;
		}
	}
}

/**
 * \fn void getTransposeOfiden(idenmat *A, iCSRmat *AT, int cols, int nNode)
 * \brief find A^T from given idenmat matrix A
 * \param *A pointer to the idenmat matrix
 * \param *AT pointer to the transpose of idenmat matrix A
 * \param cols the number of columns which will be transposed in A
 * \param nNode the number of nodes
 * \return void
 */
void getTransposeOfiden(idenmat *A, iCSRmat *AT, int cols, int nNode)
{
	int n,m;
	n=A->row;
	m=nNode;
	AT->row=m;
	AT->col=n;
	AT->IA=(int*)calloc(m+1, sizeof(int));
	AT->nnz=n*cols;
	AT->JA=(int*)calloc(n*cols, sizeof(int));
	AT->val=NULL;
	
	int i,j,k,p;
	// find the number of nonzeros in the first m-1 columns of A 
	// and store these numbers in the array AT.IA in positions from 1 to m-1
	for(i=0;i<n;i++)
	{
		for(k=0;k<cols;k++)
		{
			j=A->val[i][k];
			if(j<m-1)
			{
				AT->IA[j+2]=AT->IA[j+2]+1;
			}
		}
	}
	
	for(i=2;i<=m;i++)
	{
		AT->IA[i]=AT->IA[i]+AT->IA[i-1];
	}
	
	int l;
	for(i=0;i<n;i++)
	{
		for(p=0;p<cols;p++)
		{
			j=A->val[i][p]+1;
			k=AT->IA[j];
			AT->JA[k]=i;
			AT->IA[j]=k+1;
		}
	}
}

/**
* \fn void getTransposeOfELEMENT(ELEMENT *A, iCSRmat *AT, int cols, int nNode)
* \brief find A^T from given idenmat matrix A
* \param *A pointer to the ELEMENT matrix
* \param *AT pointer to the transpose of idenmat matrix A
* \param cols the number of columns which will be transposed in A
* \param nNode the number of nodes
* \return void
*/
void getTransposeOfELEMENT(ELEMENT *A, iCSRmat *AT, int cols, int nNode)
{
	int n, m;
	n = A->row;
	m = nNode;
	AT->row = m;
	AT->col = n;
	AT->nnz = n*cols;
	AT->IA = (int*)calloc(m + 1, sizeof(int));
	AT->JA = (int*)calloc(AT->nnz, sizeof(int));
	AT->val = NULL;

	int i, j, k, p;
	// find the number of nonzeros in the first m-1 columns of A 
	// and store these numbers in the array AT.IA in positions from 1 to m-1
	for (i = 0; i<n; i++)
	{
		for (k = 0; k<cols; k++)
		{
			j = A->val[i][k];
			if (j<m - 1)
			{
				AT->IA[j + 2] = AT->IA[j + 2] + 1;
			}
		}
	}

	for (i = 2; i <= m; i++)
	{
		AT->IA[i] = AT->IA[i] + AT->IA[i - 1];
	}

	//	int l;
	for (i = 0; i<n; i++)
	{
		for (p = 0; p<cols; p++)
		{
			j = A->val[i][p] + 1;
			k = AT->IA[j];
			AT->JA[k] = i;
			AT->IA[j] = k + 1;
		}
	}
}

/**
 * \fn void getTransposeOfelementDoF(ELEMENT_DOF *A, iCSRmat *AT, int cols)
 * \brief find A^T from given ELEMENT_DOF matrix A
 * \param *A pointer to the ELEMENT_DOF matrix
 * \param *AT pointer to the transpose of ELEMENT_DOF matrix A
 * \param cols the number of columns which will be transposed in A
 * \return void
 */
void getTransposeOfelementDoF(ELEMENT_DOF *A, iCSRmat *AT, int cols)
{
	int n,m;
	if(cols<1 || cols>A->col)
		cols=A->col;

	n=A->row;
	m=A->dof;
	AT->row=m;
	AT->col=n;
	AT->IA=(int*)calloc(m+1, sizeof(int));
	AT->nnz=n*cols;
	AT->JA=(int*)calloc(AT->nnz, sizeof(int));
	AT->val=NULL;
	
	int i,j,k,p;
	// find the number of nonzeros in the first m-1 columns of A 
	// and store these numbers in the array AT.IA in positions from 1 to m-1
	for(i=0;i<n;i++)
	{
		for(k=0;k<cols;k++)
		{
			j=A->val[i][k];
			if(j<m-1)
			{
				AT->IA[j+2]=AT->IA[j+2]+1;
			}
		}
	}
	
	for(i=2;i<=m;i++)
	{
		AT->IA[i]=AT->IA[i]+AT->IA[i-1];
	}
	
	int l;
	for(i=0;i<n;i++)
	{
		for(p=0;p<cols;p++)
		{
			j=A->val[i][p]+1;
			k=AT->IA[j];
			AT->JA[k]=i;
			AT->IA[j]=k+1;
		}
	}
}

/**
 * \fn void sparseAddition(dCSRmat *A, dCSRmat *B, dCSRmat *C)
 * \brief Sparse matrix summation C=A+B
 * \param *A pointer to the dCSRmat matrix
 * \param *B pointer to the dCSRmat matrix
 * \param *C pointer to dCSRmat matrix equal to A+B
 * \return 1 or 0
 */
int sparseAddition(dCSRmat *A, dCSRmat *B, dCSRmat *C)
{
	if(A->row!=B->row || A->col!=B->col)
	{
		C=NULL;
		return 0;
	}

	int i,j,k,l;
	C->row=A->row;
	C->col=A->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int count, istart;
	int *index;
	index=(int*)calloc(C->col, sizeof(int));
	for(i=0;i<C->col;i++)
		index[i]=-1;

	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			count++;
		} // k
		for(k=B->IA[i];k<B->IA[i+1];k++)
		{
			if(index[B->JA[k]]==-1)
			{
				index[B->JA[k]]=istart;
				istart=B->JA[k];
				count++;
			}
		} // k

		C->IA[i+1]=count;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i	
	
	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];

	// step 2: Find the structure JA and val of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	C->val=(double*)calloc(C->nnz,sizeof(double));
	count=0;
	for(i=0;i<C->row;i++)
	{
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			C->JA[count]=istart;
			C->val[count]+=A->val[k];
			count++;
		} // k
		for(k=B->IA[i];k<B->IA[i+1];k++)
		{
			if(index[B->JA[k]]==-1)
			{
				index[B->JA[k]]=istart;
				istart=B->JA[k];
				C->JA[count]=istart;
				C->val[count]+=B->val[k];
				count++;
			}
			else
			{
				j=index[B->JA[k]];
				l=0;
				while(j!=-2)
				{
					j=index[j];
					l++;
				}
				C->val[C->IA[i]+l]+=B->val[k];
			}						
		} // k

		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}
	}	
	free(index);

	return 1;
}

/**
 * \fn void sparseSubtraction(dCSRmat *A, dCSRmat *B, dCSRmat *C)
 * \brief Sparse matrix difference C=A-B
 * \param *A pointer to the dCSRmat matrix
 * \param *B pointer to the dCSRmat matrix
 * \param *C pointer to dCSRmat matrix equal to A-B
 * \return 1 or 0
 */
int sparseSubtraction(dCSRmat *A, dCSRmat *B, dCSRmat *C)
{
	if(A->row!=B->row || A->col!=B->col)
	{
		C=NULL;
		return 0;
	}

	int i,j,k,l;
	C->row=A->row;
	C->col=A->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int count, istart;
	int *index;
	index=(int*)calloc(C->col, sizeof(int));
	for(i=0;i<C->col;i++)
		index[i]=-1;

	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			count++;
		} // k
		for(k=B->IA[i];k<B->IA[i+1];k++)
		{
			if(index[B->JA[k]]==-1)
			{
				index[B->JA[k]]=istart;
				istart=B->JA[k];
				count++;
			}
		} // k

		C->IA[i+1]=count;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i	
	
	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];

	// step 2: Find the structure JA and val of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	C->val=(double*)calloc(C->nnz,sizeof(double));
	count=0;
	for(i=0;i<C->row;i++)
	{
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			C->JA[count]=istart;
			C->val[count]+=A->val[k];
			count++;
		} // k
		for(k=B->IA[i];k<B->IA[i+1];k++)
		{
			if(index[B->JA[k]]==-1)
			{
				index[B->JA[k]]=istart;
				istart=B->JA[k];
				C->JA[count]=istart;
				C->val[count]-=B->val[k];
				count++;
			}
			else
			{
				j=index[B->JA[k]];
				l=0;
				while(j!=-2)
				{
					j=index[j];
					l++;
				}
				C->val[C->IA[i]+l]-=B->val[k];
			}						
		} // k

		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}
	}	
	free(index);

	return 1;
}

/**
 * \fn int dCSRPlusdDiagVector(dCSRmat *A, dvector *B, dCSRmat *C)
 * \brief Sparse matrix plus diagonal marix C=A+B
 * \param *A pointer to the dCSRmat matrix
 * \param *B pointer to the dvector matrix
 * \param *C pointer to dCSRmat matrix equal to A+B
 * \return 1 or 0
 */
int dCSRPlusdDiagVector(dCSRmat *A, dvector *B, dCSRmat *C)
{
	if(A->row!=B->row || A->col!=B->row)
	{
		C=NULL;
		return 0;
	}

	int i,j,k,l;
	C->row=A->row;
	C->col=A->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int count, istart;
	int *index;
	index=(int*)calloc(C->col, sizeof(int));
	for(i=0;i<C->col;i++)
		index[i]=-1;

	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			count++;
		} // k
		
		if(index[i]==-1)
		{
			index[i]=istart;
			istart=i;
			count++;
		}

		C->IA[i+1]=count;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i	
	
	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];

	// step 2: Find the structure JA and val of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	C->val=(double*)calloc(C->nnz,sizeof(double));
	count=0;
	for(i=0;i<C->row;i++)
	{
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			C->JA[count]=istart;
			C->val[count]+=A->val[k];
			count++;
		} // k
		
		if(index[i]==-1)
		{
			index[i]=istart;
			istart=i;
			C->JA[count]=istart;
			C->val[count]+=B->val[i];
			count++;
		}
		else
		{
			j=index[i];
			l=0;
			while(j!=-2)
			{
				j=index[j];
				l++;
			}
			C->val[C->IA[i]+l]+=B->val[i];
		}

		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}
	}	
	free(index);

	return 1;
}

/**
* \fn int dCSRMultiplydDiagVector(dCSRmat *A, dvector *D, dCSRmat *C)
* \brief Sparse matrix multiply diagonal marix C=AD
* \param *A pointer to the dCSRmat matrix
* \param *D pointer to the dvector matrix
* \param *C pointer to dCSRmat matrix equal to AD
* \return 1 or 0
*/
int dCSRMultiplydDiagVector(dCSRmat *A, dvector *D, dCSRmat *C)
{
	if (A->col != D->row)
	{
		C = NULL;
		return 0;
	}

	int i, j, k;
	C->row = A->row;
	C->col = A->col;
	C->nnz = A->nnz;
	C->IA = (int*)calloc(C->row + 1, sizeof(int));
	C->JA = (int*)calloc(C->nnz, sizeof(int));
	C->val = (double*)calloc(C->nnz, sizeof(double));

	copy_iarray(C->row + 1, A->IA, C->IA);
	copy_iarray(C->nnz, A->JA, C->JA);

	for (i = 0; i < C->row; i++)
	{
		for (k = C->IA[i]; k < C->IA[i + 1]; k++)
		{
			j = C->JA[k];
			C->val[k] = A->val[k] * D->val[j];
		}
	}

	return 1;
}

/**
* \fn int dDiagVectorMultiplydCSR(dvector *D, dCSRmat *A, dCSRmat *C)
* \brief Sparse matrix multiply diagonal marix C=DA
* \param *D pointer to the dvector matrix
* \param *A pointer to the dCSRmat matrix
* \param *C pointer to dCSRmat matrix equal to DA
* \return 1 or 0
*/
int dDiagVectorMultiplydCSR(dvector *D, dCSRmat *A, dCSRmat *C)
{
	if (A->row != D->row)
	{
		C = NULL;
		return 0;
	}

	int i, j, k;
	C->row = A->row;
	C->col = A->col;
	C->nnz = A->nnz;
	C->IA = (int*)calloc(C->row + 1, sizeof(int));
	C->JA = (int*)calloc(C->nnz, sizeof(int));
	C->val = (double*)calloc(C->nnz, sizeof(double));

	copy_iarray(C->row + 1, A->IA, C->IA);
	copy_iarray(C->nnz, A->JA, C->JA);

	for (i = 0; i < C->row; i++)
	{
		for (k = C->IA[i]; k < C->IA[i + 1]; k++)
			C->val[k] = D->val[i]*A->val[k];
	}

	return 1;
}


/**
 * \fn int dCSRPlusdBD(dCSRmat *A, dBDmat *B, dCSRmat *C)
 * \brief Sparse matrix plus block diagonal marix C=A+B
 * \param *A pointer to the dCSRmat matrix
 * \param *B pointer to the dBDmat matrix
 * \param *C pointer to dCSRmat matrix equal to A+B
 * \return 1 or 0
 */
int dCSRPlusdBD(dCSRmat *A, dBDmat *B, dCSRmat *C)
{
	if(A->row!=B->row || A->col!=B->col)
	{
		C=NULL;
		return 0;
	}

	int i,j,k,l,curb;
	C->row=A->row;
	C->col=A->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int lds=B->row/B->nb;
	int count, istart;
	int *index;
	index=(int*)calloc(C->col, sizeof(int));
	for(i=0;i<C->col;i++)
		index[i]=-1;

	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			count++;
		} // k
		curb=i/lds;
		for(k=0;k<lds;k++)
		{
			if(index[curb*lds+k]==-1)
			{
				index[curb*lds+k]=istart;
				istart=curb*lds+k;
				count++;
			}
		} // k

		C->IA[i+1]=count;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i	
	
	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];

	// step 2: Find the structure JA and val of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	C->val=(double*)calloc(C->nnz,sizeof(double));
	count=0;
	for(i=0;i<C->row;i++)
	{
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			index[A->JA[k]]=istart;
			istart=A->JA[k];
			C->JA[count]=istart;
			C->val[count]+=A->val[k];
			count++;
		} // k
		curb=i/lds;
		for(k=0;k<lds;k++)
		{
			if(index[curb*lds+k]==-1)
			{
				index[curb*lds+k]=istart;
				istart=curb*lds+k;
				C->JA[count]=istart;
				C->val[count]+=(B->blk+curb)->val[i%lds][k];
				count++;
			}
			else
			{
				j=index[curb*lds+k];
				l=0;
				while(j!=-2)
				{
					j=index[j];
					l++;
				}
				C->val[C->IA[i]+l]+=(B->blk+curb)->val[i%lds][k];
			}
		} // k

		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}
	}	
	free(index);

	return 1;
}

/**
 * \fn int dBDMinusdCSR(dBDmat *A, dCSRmat *B, dCSRmat *C)
 * \brief block diagonal matrix minus Sparse marix C=A-B
 * \param *A pointer to the dBDmat matrix
 * \param *B pointer to the dCSRmat matrix
 * \param *C pointer to dCSRmat matrix equal to A+B
 * \return 1 or 0
 */
int dBDMinusdCSR(dBDmat *A, dCSRmat *B, dCSRmat *C)
{
	if(A->row!=B->row || A->col!=B->col)
	{
		C=NULL;
		return 0;
	}

	int i,j,k,l,curb;
	C->row=A->row;
	C->col=A->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int lds=A->row/A->nb;
	int count, istart;
	int *index;
	index=(int*)calloc(C->col, sizeof(int));
	for(i=0;i<C->col;i++)
		index[i]=-1;

	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=B->IA[i];k<B->IA[i+1];k++)
		{
			index[B->JA[k]]=istart;
			istart=B->JA[k];
			count++;
		} // k
		curb=i/lds;
		for(k=0;k<lds;k++)
		{
			if(index[curb*lds+k]==-1)
			{
				index[curb*lds+k]=istart;
				istart=curb*lds+k;
				count++;
			}
		} // k

		C->IA[i+1]=count;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i	
	
	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];

	// step 2: Find the structure JA and val of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	C->val=(double*)calloc(C->nnz,sizeof(double));
	count=0;
	for(i=0;i<C->row;i++)
	{
		istart=-2;
		for(k=B->IA[i];k<B->IA[i+1];k++)
		{
			index[B->JA[k]]=istart;
			istart=B->JA[k];
			C->JA[count]=istart;
			C->val[count]-=B->val[k];
			count++;
		} // k
		curb=i/lds;
		for(k=0;k<lds;k++)
		{
			if(index[curb*lds+k]==-1)
			{
				index[curb*lds+k]=istart;
				istart=curb*lds+k;
				C->JA[count]=istart;
				C->val[count]+=(A->blk+curb)->val[i%lds][k];
				count++;
			}
			else
			{
				j=index[curb*lds+k];
				l=0;
				while(j!=-2)
				{
					j=index[j];
					l++;
				}
				C->val[C->IA[i]+l]+=(A->blk+curb)->val[i%lds][k];
			}
		} // k

		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}
	}	
	free(index);

	return 1;
}

/**
 * \fn void sparseMultiplication(dCSRmat *A, dCSRmat *B, dCSRmat *C)
 * \brief Sparse matrix multiplication C=A*B
 * \param *A pointer to the dCSRmat matrix
 * \param *B pointer to the dCSRmat matrix
 * \param *C pointer to dCSRmat matrix equal to A*B
 * \return void
 */
void sparseMultiplication(dCSRmat *A, dCSRmat *B, dCSRmat *C)
{
	int i,j,k,l;
	C->row=A->row;
	C->col=B->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int count, istart;
	int *index;
	index=(int*)calloc(B->col, sizeof(int));
	for(i=0;i<B->col;i++)
		index[i]=-1;

	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			for(j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];j++)
			{
				if(index[B->JA[j]]==-1)
				{
					index[B->JA[j]]=istart;
					istart=B->JA[j];
					count++;
				}
			}
		} // k

		C->IA[i+1]=count;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i	
	
	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];

	// step 2: Find the structure JA of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			for(j=B->IA[A->JA[k]];j<B->IA[A->JA[k]+1];j++)
			{
				if(index[B->JA[j]]==-1)
				{
					index[B->JA[j]]=istart;
					istart=B->JA[j];
					count++;
				}
			} // j
		} // k

		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			C->JA[j]=istart;
			istart=index[istart];
			index[C->JA[j]]=-1;
		}
	}	
	free(index);

	// step 3: Find the structure A of C
	C->val=(double*)calloc(C->nnz,sizeof(double));
	for(i=0;i<C->row;i++)
	{
		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			C->val[j]=0;
			for(k=A->IA[i];k<A->IA[i+1];k++)
			{
				for(l=B->IA[A->JA[k]];l<B->IA[A->JA[k]+1];l++)
				{
					if(B->JA[l]==C->JA[j])
					{
						C->val[j]+=A->val[k]*B->val[l];
					}
				}
			}
		}
	}	
}

 /**
 * \fn void dBDMultiplydCSR(dBDmat *A, dCSRmat *B, dCSRmat *C)
 * \brief block diagonal marix multiplies Sparse matrix C=A*B
 * \param *A pointer to the dBDmat matrix
 * \param *B pointer to the dCSRmat matrix
 * \param *C pointer to dCSRmat matrix equal to A*B
 * \return void
 */
void dBDMultiplydCSR(dBDmat *A, dCSRmat *B, dCSRmat *C)
{
	int i,j,k;
	int i1,j1,k1;
	int nblk;
	int nb=A->nb;
	ddenmat *blk;

	if(A->col!=B->row)
	{
		C->row=0;
		C->col=0;
		C->nnz=0;
		C->val=NULL;
		C->JA=NULL;
		C->IA=NULL;
		return;
	}

	C->row=A->row;
	C->col=B->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int count, istart;
	int *index;
	index=(int*)calloc(B->col, sizeof(int));
	for(i=0;i<B->col;i++)
		index[i]=-1;
	
	// step 1: Find first the structure IA of C
	int rcount=0;
	int lcount=0;
	for(nblk=0;nblk<nb;nblk++)
	{
		blk=A->blk+nblk;

		istart=-2;
		count=0;
		for(k1=0;k1<blk->col;k1++)
		{
			k = lcount+k1;
			for(j1=B->IA[k];j1<B->IA[k+1];j1++)
			{
				j=B->JA[j1];

				if(index[j]==-1)
				{
					index[j]=istart;
					istart=j;
					count++;
				}
			}
		} // k1

		for(i1=0;i1<blk->row;i1++)
			C->IA[rcount+i1+1]=count;

		for(j=0;j<count;j++)
		{
			j1=istart;
			istart=index[j1];
			index[j1]=-1;
		}

		rcount += blk->row;
		lcount += blk->col;
	}

	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];
	
	// step 2: Find the structure JA of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	rcount=0;
	lcount=0;
	for(nblk=0;nblk<nb;nblk++)
	{
		blk=A->blk+nblk;
		
		i=rcount;
		istart=-2;
		count=0;
		for(k1=0;k1<blk->col;k1++)
		{
			k = lcount+k1;
			for(j1=B->IA[k];j1<B->IA[k+1];j1++)
			{
				j=B->JA[j1];

				if(index[j]==-1)
				{
					index[j]=istart;
					istart=j;
					count++;
				}
			}
		} // k1

		for(j=C->IA[rcount];j<C->IA[rcount+1];j++)
		{
			C->JA[j]=istart;
			istart=index[istart];
			index[C->JA[j]]=-1;
		}

		for(i1=1;i1<blk->row;i1++)
		{
			i=rcount+i1;
			for(j=C->IA[i];j<C->IA[i+1];j++)
				C->JA[j]=C->JA[C->IA[rcount]+j-C->IA[i]];
		}

		rcount += blk->row;
		lcount += blk->col;
	}

	// step 3: Find the structure A of C
	C->val=(double*)calloc(C->nnz,sizeof(double));
	rcount=0;
	lcount=0;
	for(nblk=0;nblk<nb;nblk++)
	{
		blk=A->blk+nblk;
		for(i1=0;i1<blk->row;i1++)
		{
			i=rcount+i1;
			for(j=C->IA[i];j<C->IA[i+1];j++)
				index[C->JA[j]]=j;
			for(k1=0;k1<blk->col;k1++)
			{
				k=lcount+k1;
				for(j1=B->IA[k];j1<B->IA[k+1];j1++)
				{
					j=B->JA[j1];
					C->val[index[j]]+=blk->val[i1][k1]*B->val[j1];
				}
			} // k1
		} // i1

		rcount += blk->row;
		lcount += blk->col;
	}

	free(index);
}

/**
 * \fn void dBD2MultiplydCSR(dBDmat *A, int srow, int scol, dCSRmat *B, dCSRmat *C)
 * \brief block diagonal marix multiplies Sparse matrix C=A*B
 * \param *A pointer to the dBDmat matrix
 * \param srow number of rows of the left-top subblock in each block of A
 * \param scol number of columns of the left-top subblock in each block of A
 * \param *B pointer to the dCSRmat matrix
 * \param *C pointer to dCSRmat matrix equal to A*B
 * \return void
 */
void dBD2MultiplydCSR(dBDmat *A, int srow, int scol, dCSRmat *B, dCSRmat *C)
{
	int i,j,k;
	int i1,j1,k1;
	int nblk;
	int nb=A->nb;
	ddenmat *blk;

	if(A->col!=B->row)
	{
		C->row=0;
		C->col=0;
		C->nnz=0;
		C->val=NULL;
		C->JA=NULL;
		C->IA=NULL;
		return;
	}

	C->row=A->row;
	C->col=B->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int count, istart;
	int *index;
	index=(int*)calloc(B->col, sizeof(int));
	for(i=0;i<B->col;i++)
		index[i]=-1;
	
	// step 1: Find first the structure IA of C
	int rcount=0;
	int lcount=0;
	for(nblk=0;nblk<nb;nblk++)
	{
		blk=A->blk+nblk;

		istart=-2;
		count=0;
		for(k1=0;k1<scol;k1++)
		{
			k = lcount+k1;
			for(j1=B->IA[k];j1<B->IA[k+1];j1++)
			{
				j=B->JA[j1];

				if(index[j]==-1)
				{
					index[j]=istart;
					istart=j;
					count++;
				}
			}
		} // k1

		for(i1=0;i1<srow;i1++)
			C->IA[rcount+i1+1]=count;

		for(j=0;j<count;j++)
		{
			j1=istart;
			istart=index[j1];
			index[j1]=-1;
		}

		istart=-2;
		count=0;
		for(k1=scol;k1<blk->col;k1++)
		{
			k = lcount+k1;
			for(j1=B->IA[k];j1<B->IA[k+1];j1++)
			{
				j=B->JA[j1];

				if(index[j]==-1)
				{
					index[j]=istart;
					istart=j;
					count++;
				}
			}
		} // k1

		for(i1=srow;i1<blk->row;i1++)
			C->IA[rcount+i1+1]=count;

		for(j=0;j<count;j++)
		{
			j1=istart;
			istart=index[j1];
			index[j1]=-1;
		}

		rcount += blk->row;
		lcount += blk->col;
	}

	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];
	
	// step 2: Find the structure JA of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	rcount=0;
	lcount=0;
	for(nblk=0;nblk<nb;nblk++)
	{
		blk=A->blk+nblk;
	
		istart=-2;
		count=0;
		for(k1=0;k1<scol;k1++)
		{
			k = lcount+k1;
			for(j1=B->IA[k];j1<B->IA[k+1];j1++)
			{
				j=B->JA[j1];

				if(index[j]==-1)
				{
					index[j]=istart;
					istart=j;
					count++;
				}
			}
		} // k1

		for(j=C->IA[rcount];j<C->IA[rcount+1];j++)
		{
			C->JA[j]=istart;
			istart=index[istart];
			index[C->JA[j]]=-1;
		}

		for(i1=1;i1<srow;i1++)
		{
			i=rcount+i1;
			for(j=C->IA[i];j<C->IA[i+1];j++)
				C->JA[j]=C->JA[C->IA[rcount]+j-C->IA[i]];
		}

		istart=-2;
		count=0;
		for(k1=scol;k1<blk->col;k1++)
		{
			k = lcount+k1;
			for(j1=B->IA[k];j1<B->IA[k+1];j1++)
			{
				j=B->JA[j1];

				if(index[j]==-1)
				{
					index[j]=istart;
					istart=j;
					count++;
				}
			}
		} // k1

		for(j=C->IA[rcount+srow];j<C->IA[rcount+srow+1];j++)
		{
			C->JA[j]=istart;
			istart=index[istart];
			index[C->JA[j]]=-1;
		}

		for(i1=srow+1;i1<blk->row;i1++)
		{
			i=rcount+i1;
			for(j=C->IA[i];j<C->IA[i+1];j++)
				C->JA[j]=C->JA[C->IA[rcount+srow]+j-C->IA[i]];
		}		
		
		rcount += blk->row;
		lcount += blk->col;
	}

	// step 3: Find the structure A of C
	C->val=(double*)calloc(C->nnz,sizeof(double));
	rcount=0;
	lcount=0;
	for(nblk=0;nblk<nb;nblk++)
	{
		blk=A->blk+nblk;
		for(i1=0;i1<srow;i1++)
		{
			i=rcount+i1;
			for(j=C->IA[i];j<C->IA[i+1];j++)
				index[C->JA[j]]=j;
			for(k1=0;k1<scol;k1++)
			{
				k=lcount+k1;
				for(j1=B->IA[k];j1<B->IA[k+1];j1++)
				{
					j=B->JA[j1];
					C->val[index[j]]+=blk->val[i1][k1]*B->val[j1];
				}
			} // k1
		} // i1

		for(i1=srow;i1<blk->row;i1++)
		{
			i=rcount+i1;
			for(j=C->IA[i];j<C->IA[i+1];j++)
				index[C->JA[j]]=j;
			for(k1=scol;k1<blk->col;k1++)
			{
				k=lcount+k1;
				for(j1=B->IA[k];j1<B->IA[k+1];j1++)
				{
					j=B->JA[j1];
					C->val[index[j]]+=blk->val[i1][k1]*B->val[j1];
				}
			} // k1
		} // i1

		rcount += blk->row;
		lcount += blk->col;
	}

	free(index);
}

/**
 * \fn void dCSRMultiplydBD(dCSRmat *A, dBDmat *B, dCSRmat *C)
 * \brief Sparse matrix multiplies block diagonal marix C=A*B
 * \param *A pointer to the dCSRmat matrix
 * \param *B pointer to the dBDmat matrix
 * \param *C pointer to dCSRmat matrix equal to A*B
 * \return void
 */
void dCSRMultiplydBD(dCSRmat *A, dBDmat *B, dCSRmat *C)
{
	int i,j,k,l,j1;
	int nblk;
	int lrs=B->row/B->nb;
	int lcs=B->col/B->nb;
	
	C->row=A->row;
	C->col=B->col;
	C->val=NULL;
	C->JA=NULL;
	C->IA=(int*)calloc(C->row+1, sizeof(int));
	
	int count, istart;
	int *index;  // index of blocks
	index=(int*)calloc(B->nb, sizeof(int));
	for(i=0;i<B->nb;i++)
		index[i]=-1;
	
	// step 1: Find first the structure IA of C
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			nblk=A->JA[k]/lrs;
			if(index[nblk]==-1)
			{
				index[nblk]=istart;
				istart=nblk;
				count++;
			}
		} // k

		C->IA[i+1]=count*lcs;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i
	
	for(i=0;i<C->row;i++)
		C->IA[i+1]+=C->IA[i];
	C->nnz=C->IA[C->row];
	
	// step 2: Find the structure JA of C
	C->JA=(int*)calloc(C->nnz,sizeof(int));
	for(i=0;i<C->row;i++)
	{
		count=0;
		istart=-2;
		for(k=A->IA[i];k<A->IA[i+1];k++)
		{
			nblk=A->JA[k]/lrs;
			if(index[nblk]==-1)
			{
				index[nblk]=istart;
				istart=nblk;
				count++;
			}
		} // k

		for(j=0;j<count;j++)
		{
			l=istart;
			for(j1=0;j1<lcs;j1++)
				C->JA[C->IA[i]+lcs*j+j1]=lcs*l+j1;
			istart=index[l];
			index[l]=-1;
		}		
	} // i

	// step 3: Find the structure A of C
	C->val=(double*)calloc(C->nnz,sizeof(double));
	for(i=0;i<C->row;i++)
	{
		for(j=C->IA[i];j<C->IA[i+1];j++)
		{
			nblk=C->JA[j]/lcs;

			for(k=A->IA[i];k<A->IA[i+1];k++)
			{
				if(nblk == A->JA[k]/lrs)
				{
					C->val[j]+=A->val[k]*B->blk[nblk].val[A->JA[k]%lrs][C->JA[j]%lcs];
				}
			}
		}
	}

	free(index);
}


/**
 * \fn sparseTripleMultiplication(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
 * \brief Triple sparse matrix multiplication B=R*A*P
 *
 * this triple multiplication is faster than the following 3 conrrespoding ones
 * the algorithm of forming B.IA, B.JA is provided by 
 * ref:	Sparse Matrix Triple Multiply,Captain Trips Implementation CS521, Spring 2002 
 *	 SM4:Adam Zornes,Chris Rambicure,Danny Thorne,Chengdong Li,Vandana Chopra.
 * the algorithm of forming B.val is similar as that of B.JA, but it is different from the reference 
 * whose algorithm in forming B.IA in not very fast
 * \param *R pointer to the dCSRmat matrix
 * \param *A pointer to the dCSRmat matrix
 * \param *P pointer to the dCSRmat matrix
 * \param *B pointer to dCSRmat matrix equal to R*A*P
 * \return void
 */
void sparseTripleMultiplication(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
{
	int row=R->row, col=P->col;
	int j,jj,k,length,i,istart,iistart,count;
	int index[A->col];
	int *iindex=(int*)calloc(col,sizeof(int));
	//int iindex[col];
	B->row=row;
	B->col=col;
	B->IA=(int*)calloc(row+1,sizeof(int));
	B->JA=NULL;
	B->val=NULL;
	
	//initialize some variables
	
	for (i=0; i < A->col; i++)
		index[i] = -2;
	
	for (i=0; i < col; i++)
		iindex[i] = -2;
	
	B->IA[0] = 0;
	
	// here the calculation begins
	for (i=0; i < row; i++){
		// reset istart and length at the begining of each loop
		istart = -1;
		length = 0;
		
		//go across the rows in R
		for ( jj = R->IA[i]; jj < R->IA[i+1]; jj++ ){
			j = R->JA[jj];
			//for each column in A
			for (k=A->IA[j]; k < A->IA[j+1]; k++){
				if (index[A->JA[k]] == -2){
					index[A->JA[k]] = istart;
					istart = A->JA[k];
					length++;
				}
			}
		}    
		
		//book-keeping [reseting length and setting iistart]
		count = length;
		length = 0;
		iistart = -1;
		// use each column that would have resulted from R*A
		for (j=0; j < count; j++){
			jj = istart;
			istart = index[istart];
			index[jj] = -2;
			
			// go across the row of P
			for (k=P->IA[jj]; k < P->IA[jj+1]; k++){
				//pull out the appropriate columns of P
				if (iindex[P->JA[k]] == -2){
					iindex[P->JA[k]] = iistart;
					iistart = P->JA[k];
					length++;
				}
			}
		}
		
		// set B->IA
		B->IA[i+1]= B->IA[i]+length;
		
		B->JA=(int*)realloc(B->JA, (B->IA[i+1])*sizeof(int));
		// put the correct columns of p into the column list of the products
		for (j=B->IA[i]; j< B->IA[i+1]; j++){
			//put the value in B->JA
			B->JA[j] = iistart;
			// set istart to the next value
			iistart = iindex[iistart];
			// set the iindex spot to 0
			iindex[B->JA[j]] = -2;
		}
		
	}

	free(iindex);
	
	B->nnz=B->IA[row];
	B->val=(double*)calloc(B->IA[row],sizeof(double));
	
	double *temp;
	int *BTindex;
	
	temp=(double*)calloc(A->col,sizeof(double));
	BTindex=(int*)calloc(B->col,sizeof(int));
	
	// across each row in r*a*P
	for (i=0; i < row; i++){
		
		for( j = B->IA[i]; j < B->IA[i+1]; j++)
		{
			BTindex[B->JA[j]]=j;
		}
		// reset istart and length at the begining of each loop
		istart = -1;
		length = 0;
		//go across the rows in R
		for ( jj = R->IA[i]; jj < R->IA[i+1]; jj++ ){
			j = R->JA[jj];
			
			//for each column in A
			for (k=A->IA[j]; k < A->IA[j+1]; k++){
				if (index[A->JA[k]] == -2){
					index[A->JA[k]] = istart;
					istart = A->JA[k];
					length++;
				}
				temp[A->JA[k]]+=R->val[jj]*A->val[k];
			}
		} 
		
		//book-keeping [reseting length and setting iistart]
		// use each column that would have resulted from R*A
		for (j=0; j < length; j++){
			jj = istart;
			istart = index[istart];
			index[jj] = -2;
			
			// go across the row of P
			for (k=P->IA[jj]; k < P->IA[jj+1]; k++){
				//pull out the appropriate columns of P
				B->val[BTindex[P->JA[k]]]+=temp[jj]*P->val[k];
			}
			temp[jj]=0;
		}
	}
	
	free(temp);
	free(BTindex);
}

/**
 * \fn sparseTripleMultiplication1(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
 * \brief triple sparse matrix multiplication B=R*A*P
 *
 * form B.IA, B.JA parallelly
 * the algorithm of forming B.IA, B.JA, B.val is provided by 
 *    ref:	Sparse Matrix Triple Multiply,Captain Trips Implementation CS521, Spring 2002 
 *			SM4:Adam Zornes,Chris Rambicure,Danny Thorne,Chengdong Li,Vandana Chopra.
 * \param *R pointer to the dCSRmat matrix
 * \param *A pointer to the dCSRmat matrix
 * \param *P pointer to the dCSRmat matrix
 * \param *B pointer to dCSRmat matrix equal to R*A*P
 * \return void
 */
void sparseTripleMultiplication1(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
{
	int row=R->row, col=P->col;
	int j,jj,k,length,i,istart,iistart,count, fcol, Ai, m, find;
	int index[A->col];
	int iindex[col];
	double temp, sum, found;
	
	B->row=row;
	B->col=col;
	B->IA=(int*)calloc(row+1,sizeof(int));
	B->JA=NULL;
	B->val=NULL;
	
	//initialize some variables
	
	for (i=0; i < A->col; i++)
		index[i] = -2;
	
	for (i=0; i < col; i++)
		iindex[i] = -2;
	
	B->IA[0] = 0;
	
	// here the calculation begins
	for (i=0; i < row; i++){
		// reset istart and length at the begining of each loop
		istart = -1;
		length = 0;
		
		//go across the rows in R
		for ( jj = R->IA[i]; jj < R->IA[i+1]; jj++ ){
			j = R->JA[jj];
			
			//for each column in A
			for (k=A->IA[j]; k < A->IA[j+1]; k++){
				if (index[A->JA[k]] == -2){
					index[A->JA[k]] = istart;
					istart = A->JA[k];
					length++;
				}
			}
			
		}    
		
		//book-keeping [reseting length and setting iistart]
		count = length;
		length = 0;
		iistart = -1;
		// use each column that would have resulted from R*A
		for (j=0; j < count; j++){
			jj = istart;
			istart = index[istart];
			index[jj] = -2;
			
			// go across the row of P
			for (k=P->IA[jj]; k < P->IA[jj+1]; k++){
				//pull out the appropriate columns of P
				if (iindex[P->JA[k]] == -2){
					iindex[P->JA[k]] = iistart;
					iistart = P->JA[k];
					length++;
				}
			}
			
		}
		
		// set B->IA
		B->IA[i+1]= B->IA[i]+length;
		
		B->JA=(int*)realloc(B->JA, (B->IA[i+1])*sizeof(int));
		// put the correct columns of p into the column list of the products
		for (j=B->IA[i]; j< B->IA[i+1]; j++){
			//put the value in B->JA
			B->JA[j] = iistart;
			// set istart to the next value
			iistart = iindex[iistart];
			// set the iindex spot to 0
			iindex[B->JA[j]] = -2;
		}
		
	}
	
	B->val=(double*)calloc(B->IA[row],sizeof(double));
	
	//traspose P so that it is possible to go across its row first
	dCSRmat PT;
	getTransposeOfSparse(P,&PT);
	
	// across each row in r*a*P
	for (i=0; i < row; i++){
		
		//for each nz entry in the row
		for (j=B->IA[i]; j < B->IA[i+1]; j++){
			//here is where the calculation is done            
			sum = 0;
			//sum across the row of PT
			fcol = B->JA[j];
			for (k=PT.IA[fcol]; k < PT.IA[fcol+1]; k++){
				//calculate the corresponding entry in the product of R*A
				//sum across the ith row of r with the fcol column of A
				temp = 0;
				for (jj=R->IA[i]; jj < R->IA[i+1]; jj++){
					//see if the corresponding entry in a exists (it will have the value of find)
					find = PT.JA[k];
					found = 0;
					Ai = R->JA[jj];
					//the row of A to look in is the column of R->  go across it and see if the row (the row of PT) is there
					for (m=A->IA[Ai]; m < A->IA[Ai+1]; m++){
						//if it is, then find becomes its value
						if (A->JA[m] == find){
							found = A->val[m];
							break;
						} // if
					} // for m
					
					//sum for the product entry in R*A
					temp+= R->val[jj]*found;
				} // for jj
				
				//sum for the product entry in R*A*P                
				sum+=PT.val[k]*temp;
			} // for k
			
			//save the total sum
			B->val[j] = sum;
			
			//end calculation   
			
		} // for j
		
	}
	
	free(PT.IA);
	free(PT.JA);
	free(PT.val);
}


/**
 * \fn sparseTripleMultiplication2(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
 * \brief triple sparse matrix multiplication B=R*A*P
 *
 * form B.IA, B.JA sequecencely
 * the algorithm of forming B.IA, B.JA, B.val is provided by 
 *    ref:	Sparse Matrix Triple Multiply,Captain Trips Implementation CS521, Spring 2002 
 *			SM4:Adam Zornes,Chris Rambicure,Danny Thorne,Chengdong Li,Vandana Chopra.
 * \param *R pointer to the dCSRmat matrix
 * \param *A pointer to the dCSRmat matrix
 * \param *P pointer to the dCSRmat matrix
 * \param *B pointer to dCSRmat matrix equal to R*A*P
 * \return void
 */
void sparseTripleMultiplication2(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
{
	int row=R->row, col=P->col;
	int j,jj,k,length,i,istart,iistart,count, fcol, Ai, m, find;
	int index[A->col];
	int iindex[col];
	double temp, sum, found;
	
	B->row=row;
	B->col=col;
	B->IA=(int*)calloc(row+1,sizeof(int));
	B->JA=NULL;
	B->val=NULL;
	
	//initialize some variables
	
	for (i=0; i < A->col; i++)
		index[i] = -2;
	
	for (i=0; i < col; i++)
		iindex[i] = -2;
	
	B->IA[0] = 0;
	
	// here the calculation begins
	// form B->IA
	for (i=0; i < row; i++){
		// reset istart and length at the begining of each loop
		istart = -1;
		length = 0;
		
		//go across the rows in R
		for ( jj = R->IA[i]; jj < R->IA[i+1]; jj++ ){
			j = R->JA[jj];
			
			//for each column in A
			for (k=A->IA[j]; k < A->IA[j+1]; k++){
				if (index[A->JA[k]] == -2){
					index[A->JA[k]] = istart;
					istart = A->JA[k];
					length = length+1;
				}
			}
			
		}    
		
		//book-keeping [reseting length and setting iistart]
		count = length;
		length = 0;
		iistart = -1;
		// use each column that would have resulted from R*A
		for (j=0; j < count; j++){
			jj = istart;
			istart = index[istart];
			index[jj] = -2;
			
			// go across the row of P
			for (k=P->IA[jj]; k < P->IA[jj+1]; k++){
				//pull out the appropriate columns of P
				if (iindex[P->JA[k]] == -2){
					iindex[P->JA[k]] = iistart;
					iistart = P->JA[k];
					length = length+1;
				}
			}
			
		}
		
		// set B->IA
		B->IA[i+1]= B->IA[i]+length;
		
		// put the correct columns of p into the column list of the products
		for (j=B->IA[i]; j< B->IA[i+1]; j++){
			jj=iistart;
			// set istart to the next value
			iistart = iindex[iistart];
			// set the iindex spot to 0
			iindex[jj] = -2;
		}
	}
	
	// form B->JA
	B->JA=(int*)realloc(B->JA, (B->IA[row])*sizeof(int));
	for (i=0; i < row; i++){
		// reset istart and length at the begining of each loop
		istart = -1;
		length = 0;
		
		//go across the rows in R
		for ( jj = R->IA[i]; jj < R->IA[i+1]; jj++ ){
			j = R->JA[jj];
			
			//for each column in A
			for (k=A->IA[j]; k < A->IA[j+1]; k++){
				if (index[A->JA[k]] == -2){
					index[A->JA[k]] = istart;
					istart = A->JA[k];
					length = length+1;
				}
			}
			
		}    
		
		//book-keeping [reseting length and setting iistart]
		count = length;
		length = 0;
		iistart = -1;
		// use each column that would have resulted from R*A
		for (j=0; j < count; j++){
			jj = istart;
			istart = index[istart];
			index[jj] = -2;
			
			// go across the row of P
			for (k=P->IA[jj]; k < P->IA[jj+1]; k++){
				//pull out the appropriate columns of P
				if (iindex[P->JA[k]] == -2){
					iindex[P->JA[k]] = iistart;
					iistart = P->JA[k];
					length = length+1;
				}
			}
			
		}
		
		// put the correct columns of p into the column list of the products
		for (j=B->IA[i]; j< B->IA[i+1]; j++){
			//put the value in B->JA
			B->JA[j] = iistart;
			// set istart to the next value
			iistart = iindex[iistart];
			// set the iindex spot to 0
			iindex[B->JA[j]] = -2;
		}
	}
	
	B->val=(double*)calloc(B->IA[row],sizeof(double));
	
	//traspose P so that it is possible to go across its row first
	dCSRmat PT;
	getTransposeOfSparse(P,&PT);
	
	// across each row in r*a*P
	for (i=0; i < row; i++){
		
		//for each nz entry in the row
		for (j=B->IA[i]; j < B->IA[i+1]; j++){
			//here is where the calculation is done            
			sum = 0;
			//sum across the row of PT
			fcol = B->JA[j];
			for (k=PT.IA[fcol]; k < PT.IA[fcol+1]; k++){
				//calculate the corresponding entry in the product of R*A
				//sum across the ith row of r with the fcol column of A
				temp = 0;
				for (jj=R->IA[i]; jj < R->IA[i+1]; jj++){
					//see if the corresponding entry in a exists (it will have the value of find)
					find = PT.JA[k];
					found = 0;
					Ai = R->JA[jj];
					//the row of A to look in is the column of R->  go across it and see if the row (the row of PT) is there
					for (m=A->IA[Ai]; m < A->IA[Ai+1]; m++){
						//if it is, then find becomes its value
						if (A->JA[m] == find){
							found = A->val[m];
							break;
						} // if
					} // for m
					
					//sum for the product entry in R*A
					temp = temp + R->val[jj]*found;
				} // for jj
				
				//sum for the product entry in R*A*P                
				sum = sum + PT.val[k]*temp;
			} // for k
			
			//save the total sum
			B->val[j] = sum;
			
			//end calculation   
			
		} // for j
		
	}
	
	free(PT.IA);
	free(PT.JA);
	free(PT.val);
}


/**
 * \fn sparseTripleMultiplication3(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
 * \brief triple sparse matrix multiplication B=R*A*P using the fact that R=PT
 *
 * form B.IA, B.JA parallelly
 * the algorithm of forming B.IA, B.JA, B.val is provided by 
 *    ref:	Sparse Matrix Triple Multiply,Captain Trips Implementation CS521, Spring 2002 
 *			SM4:Adam Zornes,Chris Rambicure,Danny Thorne,Chengdong Li,Vandana Chopra.
 * \param *R pointer to the dCSRmat matrix
 * \param *A pointer to the dCSRmat matrix
 * \param *P pointer to the dCSRmat matrix
 * \param *B pointer to dCSRmat matrix equal to R*A*P
 * \return void
 */
void sparseTripleMultiplication3(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B)
{
	int row=R->row, col=P->col;
	int j,jj,k,length,i,istart,iistart,count, fcol, Ai, m, find;
	int index[A->col];
	int iindex[col];
	double temp, sum, found;
	
	B->row=row;
	B->col=col;
	B->IA=(int*)calloc(row+1,sizeof(int));
	B->JA=NULL;
	B->val=NULL;
	
	//initialize some variables
	
	for (i=0; i < A->col; i++)
		index[i] = -2;
	
	for (i=0; i < col; i++)
		iindex[i] = -2;
	
	B->IA[0] = 0;
	
	// here the calculation begins
	for (i=0; i < row; i++){
		// reset istart and length at the begining of each loop
		istart = -1;
		length = 0;
		
		//go across the rows in R
		for ( jj = R->IA[i]; jj < R->IA[i+1]; jj++ ){
			j = R->JA[jj];
			
			//for each column in A
			for (k=A->IA[j]; k < A->IA[j+1]; k++){
				if (index[A->JA[k]] == -2){
					index[A->JA[k]] = istart;
					istart = A->JA[k];
					length++;
				}
			}
			
		}    
		
		//book-keeping [reseting length and setting iistart]
		count = length;
		length = 0;
		iistart = -1;
		// use each column that would have resulted from R*A
		for (j=0; j < count; j++){
			jj = istart;
			istart = index[istart];
			index[jj] = -2;
			
			// go across the row of P
			for (k=P->IA[jj]; k < P->IA[jj+1]; k++){
				//pull out the appropriate columns of P
				if (iindex[P->JA[k]] == -2){
					iindex[P->JA[k]] = iistart;
					iistart = P->JA[k];
					length++;
				}
			}
			
		}
		
		// set B->IA
		B->IA[i+1]= B->IA[i]+length;
		
		B->JA=(int*)realloc(B->JA, (B->IA[i+1])*sizeof(int));
		// put the correct columns of p into the column list of the products
		for (j=B->IA[i]; j< B->IA[i+1]; j++){
			//put the value in B->JA
			B->JA[j] = iistart;
			// set istart to the next value
			iistart = iindex[iistart];
			// set the iindex spot to 0
			iindex[B->JA[j]] = -2;
		}
		
	}
	
	B->val=(double*)calloc(B->IA[row],sizeof(double));
	
	// across each row in r*a*P
	for (i=0; i < row; i++){
		
		//for each nz entry in the row
		for (j=B->IA[i]; j < B->IA[i+1]; j++){
			
			
			//here is where the calculation is done            
			sum = 0;
			//sum across the row of R
			fcol = B->JA[j];
			for (k=R->IA[fcol]; k < R->IA[fcol+1]; k++){
				//calculate the corresponding entry in the product of R*A
				//sum across the ith row of r with the fcol column of A
				temp = 0;
				for (jj=R->IA[i]; jj < R->IA[i+1]; jj++){
					//see if the corresponding entry in a exists (it will have the value of find)
					find = R->JA[k];
					found = 0;
					Ai = R->JA[jj];
					//the row of A to look in is the column of R->  go across it and see if the row (the row of PT) is there
					for (m=A->IA[Ai]; m < A->IA[Ai+1]; m++){
						//if it is, then find becomes its value
						if (A->JA[m] == find){
							found = A->val[m];
							break;
						} // if
					} // for m
					
					//sum for the product entry in R*A
					temp+= R->val[jj]*found;
				} // for jj
				
				//sum for the product entry in R*A*P                
				sum+=R->val[k]*temp;
			} // for k
			
			//save the total sum
			B->val[j] = sum;
			
			//end calculation   
			
		} // for j
		
	}
}


/**
 * \fn int dIJtoCSR(dIJmat *A, dCSRmat *B)
 * \brief Transform a double matrix from its IJ format to its CSR format.
 * \param *A pointer to IJ matrix
 * \param *B pointer to CSR matrix
 * \return 1 if succeed, 0 if fail 
 */
int dIJtoCSR(dIJmat *A, dCSRmat *B)
{
	int i, m=A->row, n=A->col, nnz=A->nnz;
	
	int *iz=malloc(m*sizeof(int));
	
	for (i=0;i<m;i++) iz[i]=0;
	
	for (i=0;i<nnz;i++) iz[A->I_dij[i]]++; // number of nonzeros in each row
	
	B->row=m; B->col=n; B->nnz=nnz;
	
	B->IA=(int *)calloc(m+1,sizeof(int));
	B->IA[0]=0;
	for (i=1;i<m+1;i++) B->IA[i]=B->IA[i-1]+iz[i-1];
	
	B->JA=(int *)calloc(nnz,sizeof(int));
	for (i=0;i<nnz;i++) B->JA[i]=A->J[i];
	
	B->val=(double *)calloc(nnz,sizeof(double));
	for (i=0;i<nnz;i++) B->val[i]=A->val[i];
	
	free(iz);
	
	return 1;
}


/**
 * \fn int dCSRtoIJ(dCSRmat *A, dIJmat *B)
 * \brief Transform a double matrix from its CSR format to its IJ format.
 * \param *A pointer to CSR matrix
 * \param *B pointer to IJ matrix
 * \return 1 if succeed, 0 if fail 
 */
int dCSRtoIJ(dCSRmat *A, dIJmat *B)
{
	int i, j, m=A->row, n=A->col, nnz=A->nnz;
	
	B->I_dij=calloc(nnz,sizeof(int));
	B->J=calloc(nnz,sizeof(int));
	B->val=calloc(nnz,sizeof(double));
	
	for (i=0;i<m;i++) {
		for (j=A->IA[i];j<A->IA[i+1];j++) B->I_dij[j]=i;
	}
	
	for (i=0;i<nnz;i++) {
		B->J[i]=A->JA[i];
		B->val[i]=A->val[i];
	}
	
	return 1;
}
