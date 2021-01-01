/*
 *  matvecio.c
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 03/29/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------
 */

/*! \file matvecio.c
 *  \brief Matrix-vector input/output subroutines
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "header.h"

/**
 * \fn int print_dvector(int n, dvector *u) 
 * \brief Print first n entries of a dvector
 * \param n an interger (if n=0, then print all entries)
 * \param *u pointer to a dvector
 * \return 1 if succeed, 0 if fail
 */
int print_dvector(int n, dvector *u) 
{
	int i,ret;	
	
	if (n==0) n=u->row;
	for(i=0;i<n;i++) printf("%d: %le\n", i,u->val[i]);
	
	return 1;
}

/**
 * \fn int print_ivector(int n, ivector *u) 
 * \brief Print first n entries of an ivector of int type
 * \param n an interger (if n=0, then print all entries)
 * \param *u pointer to an ivector
 * \return 1 if succeed, 0 if fail
 */
int print_ivector(int n, ivector *u) 
{
	int i,ret;	
	
	if (n==0) n=u->row;
	for(i=0;i<n;i++) printf("%d: %d\n", i,u->val[i]);
	
	return 1;
}

/**
 * \fn int read_matvec(char *filename, dCSRmat *A, dvector *b)
 * \brief Read A and b from a SINGLE disk file
 * \param *filename char for file name
 * \param *A pointer to the CSR matrix
 * \param *b pointer to the dvector
 * \return 1 if succeed, 0 if fail
 *
 * File format: 
 */
int read_matvec(char *filename, dCSRmat *A, dvector *b)
{
	FILE *inputFile;
	int i,m,n,nnz;
	
	int idata; 
	double ddata;
	
	// Open input disk file
	inputFile=fopen(filename, "r");
	if(inputFile==NULL)
	{
		printf("read_matvec: opening file %s fails!\n", filename);
		return 0;
	}
	
	// Read CSR matrix
	printf("read_matvec: reading file %s...\n", filename);
	fscanf(inputFile, "%d %d", &m, &n);
	A->row=m;
	A->col=n;
	
	A->IA=(int*)calloc(m+1, sizeof(int));
	for(i=0;i<=m;i++) {
		fscanf(inputFile, "%d", &idata);
		A->IA[i]=idata;
	}
	
	nnz=A->IA[m]-A->IA[0];		
	A->nnz=nnz;
	
	A->JA=(int*)calloc(nnz, sizeof(int));
	for(i=0;i<nnz;i++) {
		fscanf(inputFile, "%d", &idata);
		A->JA[i]=idata;
	}
	
	A->val=(double*)calloc(nnz, sizeof(double));
	for(i=0;i<nnz;i++) {
		fscanf(inputFile, "%lf", &ddata);
		A->val[i]=ddata;
	}
	
	// Read vector
	fscanf(inputFile, "%d %d", &m, &i);
	b->row=m;
	b->val=(double*)calloc(m, sizeof(double));
	for(i=0;i<m;i++) {
		fscanf(inputFile, "%lf", &ddata);
		b->val[i]=ddata;
	}
	
	fclose(inputFile);	
	return 1;
}


/**
 * \fn int read_ruth_file(char *filemat, char *filerhs, dCSRmat *B, dvector *rhs)
 * \brief Read A and b from matrix and rhs disk files for Ruth test problem
 * \param *filemat char for matrix file name
 * \param *filerhs char for rhs file name
 * \param *A pointer to the CSR matrix
 * \param *b pointer to the dvector
 * \return 1 if succeed, 0 if fail
 *
 * File format: 
 */
int read_ruth_file(char *filemat, char *filerhs, dCSRmat *A, dvector *rhs)
{
	int i,j,k,l;
	
	FILE *fp;
	
	fp=fopen(filemat,"r");
	if (fp == NULL) {
		printf("read_ruth_file: opening file %s fails!\n", filemat);
		return 0;
	}
	
	fscanf(fp,"%d",&l);
	A->row=l;
	A->col=l;
	
	A->IA=(int*)calloc(l+1, sizeof(int));
	for (i=0;i<l+1;i++) fscanf(fp,"%d",&A->IA[i]);
	A->nnz=A->IA[l]-A->IA[0];
	
	A->JA=(int*)calloc(A->nnz,sizeof(int));
	for (i=0;i<A->nnz;i++) fscanf(fp,"%d",&A->JA[i]);
	
	A->val=(double*)calloc(A->nnz,sizeof(double));
	for (i=0;i<A->nnz;i++) fscanf(fp,"%le", &A->val[i]);
	
	fclose(fp);
	
	fp=fopen(filerhs,"r");
	if (fp == NULL) {
		printf("read_ruth_file: opening file %s fails!\n", filerhs);
		return 0;
	}
	
	fscanf(fp,"%d",&l);
	rhs->row=l;
	rhs->val=(double *)calloc(l,sizeof(double));
	
	for (i=0;i<l;i++) fscanf(fp,"%le",&rhs->val[i]);
	fclose(fp);

	if (A->IA[0]!=0) /* if does not start from 0, make it start from 0 */
	   for (i=0;i<A->row+1;i++) A->IA[i]--;
	if (A->JA[0]!=0) /* if does not start from 0, make it start from 0 */
	   for (i=0;i<A->nnz;i++) A->JA[i]--;
	
	return 1;
}


extern void  read_matrix_(char *,char *,int **,int **,
												 double **,double **,int *,int *,double **);

/**
 * \fn int read_ruth_file_b(char *filemat, char *filerhs, dCSRmat *B, dvector *rhs, dvector *xapp)
 * \brief Read A and b from matrix and rhs disk files for Ruth test problem
 * \param *filemat char for matrix file name
 * \param *filerhs char for rhs file name
 * \param *A pointer to the matrix
 * \param *b pointer to the vector
 * \param *xapp pointer to the vector of application solution
 * \return 1 if succeed, 0 if fail
 *
 * File format: 
 */
//int read_ruth_file_b(char *filemat, char *filerhs, dCSRmat *A, dvector *rhs, 
//										 dvector *xsolve)
//{
//	int i,j,k,l;
	
//	int fp;
//	read_matrix_(filemat,filerhs,&A->IA,&A->JA,&A->val,&rhs->val,&A->row,
//							 &A->nnz,&xsolve->val);
//	A->col = A->row;
//	rhs->row = A->row;
//	xsolve->row = A->row;
	
	/*** make IA and JA start from 0 ***/
//	if (A->IA[0]!=0) /* if does not start from 0, make it start from 0 */
//		for (i=0;i<A->row+1;i++) A->IA[i]--;
//	if (A->JA[0]!=0) /* if does not start from 0, make it start from 0 */
//		for (i=0;i<A->nnz;i++) A->JA[i]--;
	
//	return 1;
//}


/**
 * \fn int read_IJ_matrix(char *filename, dCSRmat *A)
 * \brief Read A from matrix disk file in IJ format
 * \param *filename char for matrix file name
 * \param *A pointer to the CSR matrix
 * \return 1 if succeed, 0 if fail
 *
 * File format: 
 */
int read_IJ_matrix(char *filename, dCSRmat *A)
{
	int i,j,k,l,nnz;
	double value;
	
	dIJmat Atmp; 
	
	FILE *fp;
	fp=fopen(filename,"r");
	if (fp == NULL) {
		printf("read_IJ_matrix: opening file %s fails!\n", filename);
		return 0;
	}
	
	fscanf(fp,"%d %d",&l,&nnz);
	
	Atmp.row=l; Atmp.col=l; Atmp.nnz=nnz;
	Atmp.I=(int *)calloc(nnz, sizeof(int)); 
	Atmp.J=(int *)calloc(nnz, sizeof(int)); 
	Atmp.val=(double *)calloc(nnz, sizeof(double)); 
	
	k = 0;
	while ( fscanf(fp, "%d %d %le", &i, &j, &value) != EOF ) {
		Atmp.I[k]=i; Atmp.J[k]=j; Atmp.val[k]=value; k++;
	}
	
	dIJtoCSR(&Atmp,A);
	
	fclose(fp);
	
	free(Atmp.I); free(Atmp.J); free(Atmp.val);
	
	return 1;
}


/**
 * \fn int read_IJ_vector(char *filename, dvector *b)
 * \brief Read b from matrix disk file in IJ format
 * \param *filename char for vector file name
 * \param *b pointer to the dvector
 * \return 1 if succeed, 0 if fail
 *
 * File format: 
 */
int read_IJ_vector(char *filename, dvector *b)
{
	int i, n;
	double value;
	
	FILE *fp;
	fp=fopen(filename,"r");
	if (fp == NULL) {
		printf("read_IJ_vector: opening file %s fails!\n", filename);
		return 0;
	}
	
	fscanf(fp,"%d",&n);	
	b->row=n;
	
	b->val=(double *)calloc(n,sizeof(double)); 	
//	while ( fscanf(fp, "%d %le", &i, &value) != EOF ) b->val[i]=value;
	for(i=0;i<n;i++)
	{
		fscanf(fp, "%le", &value);
		b->val[i]=value;
	}
	
	fclose(fp);
	
	return 1;
}


/**
 * /fn int write_bmp16(const char *fname, int m, int n, const char map[])
 * /brief Write a sparse matrix structure to a BMP file
 * 
 * See definition in graphics.c. 
 */
int write_bmp16(const char *fname, int m, int n, const char map[]);
/**
 * /fn int plot_matrix(const dCSRmat *A, const char *fname, int size)
 * /brief Write sparse matrix pattern in BMP file format
 * /prama *A pointer to the dCSRmat matrix in CSR format
 * /prama size integer size*size is the picture size for the picture
 * /prama fname char for vector file name
 * 
 * The routine spm_show_mat writes pattern of the specified dCSRmat
 * matrix in uncompressed BMP file format (Windows bitmap) to a binary
 * file whose name is specified by the character string fname.
 *
 * Each pixel corresponds to one matrix element. The pixel colors have
 * the following meaning:
 *
 *  White    structurally zero element
 *  Blue     positive element
 *  Red      negative element
 *  Brown    nearly zero element
 *
 */
int plot_matrix(const dCSRmat *A, const char *fname, int size)
{     
  int m = A->row;
  int n = A->col;
  int i, j, k, l, ret;
  char *map;
	
	if (size>min(m,n)) size=min(m,n);
	
  printf("plot_matrix: writing matrix pattern to `%s'...\n",fname);
	
	map = malloc(size * size);
	memset(map, 0x0F, size * size);
	for (i = 0; i < size; i++) {
		for (j = A->IA[i]; j < A->IA[i+1]; j++) {
			if (A->JA[j]<size) {
				k = size*i + A->JA[j];
				if (map[k] != 0x0F)
					map[k] = 0x0F;
				else if (A->val[j] > 1e-20)
					map[k] = 0x09; /* bright blue */
				else if (A->val[j] < -1e-20)
					map[k] = 0x0C; /* bright red */
				else
					map[k] = 0x06; /* brown */
			}
		}
  }
	ret = write_bmp16(fname, size, size, map);
	
	free(map);
  return ret;
}


/**
 * \fn int write_IJ_matrix(dCSRmat *A, char *fname)
 * \brief Write a matrix to disk file in IJ format (coordinate format)
 * \param *A pointer to the dCSRmat matrix
 * \param *fname char for vector file name
 * \return 1 if succeed, 0 if fail
 *
 * File format: 
 *
 *	 The routine write_IJ_matrix writes the specified double vector in 
 *	 IJ format. The first line of the file gives the number of rows and
 * the number of columns. And then gives nonzero values in i,j,a(i,j)
 * order.      
 */
int write_IJ_matrix(dCSRmat *A, char *fname)
{     
  int m = A->row;
  int n = A->col;
  int i, j;
  char *map;
	
	FILE *outputFile;
	outputFile=fopen(fname, "w");
	if(outputFile==NULL)
	{
		printf("write_IJ_matrix: opening file %s fails!\n", fname);
		return 0;
	}
	
  printf("write_IJ_matrix: writing matrix to `%s'...\n",fname);
	fprintf(outputFile,"%d  %d\n",m,n);	
  for (i = 0; i < m; i++) {
		for (j = A->IA[i]; j < A->IA[i+1]; j++)
      fprintf(outputFile,"%d  %d  %f\n",i+1,A->JA[j]+1,A->val[j]);
  }
	
	fclose(outputFile);
	
  return 1;
}


/**
 * \fn int write_IJ_dvector(dvector *vec, char *fname)
 * \brief Write a dvector to disk file in IJ format (coordinate format)
 * \param *vec pointer to the dvector
 * \param *fname char for vector file name
 * \return 1 if succeed, 0 if fail
 *
 * File format: 
 *
 *	 The routine write_IJ_dvector writes the specified double vector in 
 *	 IJ format. The first line of the file is the length of the vector;
 * and after that, each line gives index and value of the entries.     
 */
int write_IJ_dvector(dvector *vec, char *fname)
{
	int m = vec->row, i;
	
	FILE *outputFile;
	outputFile=fopen(fname,"w");
	if (outputFile==NULL) {
		printf("write_IJ_vec: opening file %s fails!!\n", fname);
		return 0;
	}
	
	printf("write_IJ_vec: writing vector to '%s'...\n",fname);
	fprintf(outputFile,"%d\n",m);
	
	for (i=0;i<m;i++) fprintf(outputFile,"%d %le\n",i,vec->val[i]);
	
	fclose(outputFile);
	
	return 1;
}
