/*
 *  matvec.h
 *  Header file for sparse matrix and vector opeartions
 *  
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 03/29/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------
 */

/** 
 * \file matvec.h
 * \brief Header file for matrix-vector operations
 */
 
/* matvecio.c */
int read_matvec(char *filename, dCSRmat *A, dvector *u); 
int plot_matrix(const dCSRmat *A, const char *fname, int size);
int write_IJ_matrix(dCSRmat *A, char *fname);
int write_IJ_vector(dvector *v, char *fname);
int read_IJ_matrix(char *filename, dCSRmat *A);
int read_IJ_vector(char *filename, dvector *b);
int read_ruth_file(char *filemat, char *filerhs, dCSRmat *A, dvector *rhs);
//int read_ruth_file_b(char *filemat, char *filerhs, dCSRmat *A, dvector *rhs, 
//										 dvector *xsolve);

/* matvecop.c */
int create_csr_matrix(int m, int n, int nnz, dCSRmat *A);
int free_csr_matrix(dCSRmat *A);
int free_icsr_matrix(iCSRmat *A);
int read_csr_matrix(char *filename, dCSRmat *A);

int create_iden_matrix(int m, int n, idenmat *A);
int free_iden_matrix(idenmat *A);
int create_ELEMENT(int m, int n, ELEMENT *A);
int free_ELEMENT(ELEMENT *A);
int create_EDGE(int m, int n, EDGE *A);
int free_EDGE(EDGE *A);
int create_dden_matrix(int m, int n, ddenmat *A);
int free_dden_matrix(ddenmat *A);
int create_dden_matrix3(int l, int m, int n, ddenmat3 *A);
int free_dden_matrix3(ddenmat3 *A);
int create_dbd_matrix(int m, int n, int nb, dBDmat *A);
int free_dbd_matrix(dBDmat *A);
int create_dobd_matrix(int m, int n, int nb, dOBDmat *A);
int free_dobd_matrix(dOBDmat *A);
int create_dennode(int m, int n, dennode *A);
int free_dennode(dennode *A);
int create_elementDOF(int dop, int dof, int row, int col, ELEMENT_DOF *A);
int free_elementDOF(ELEMENT_DOF *A);

int create_dvector(int m, dvector *u);
int create_ivector(int m, ivector *u);
int free_dvector(dvector *u);
int free_ivector(ivector *u);
int print_dvector(int n, dvector *u); 
int print_ivector(int n, ivector *u);

int init_array(int n, double *x, double val);
int copy_iarray(int n, int *x, int *y);
int copy_array(int n, double *x, double *y);
int axy_array(int n, double a, double *x, double *y);
int axpy_array(int n, double a, double *x, double *y);
int axpby_array(int n, double a, double *x, double b, double *y);
int axpyz_array(int n, double a, double *x, double *y, double *z);
double dot_array(int n, double *x, double *y);

int init_dvector(dvector *x, double val);
int init_dvector2b(dvector *x, double val);
int copy_dvector(dvector *x, dvector *y);
int copy_dvector2b(dvector *x, dvector *y);
int axpy_dvector(double a, dvector *x, dvector *y);
int axpy_dvector2b(double a, dvector *x, dvector *y);
int axpyz_dvector(double a, dvector *x, dvector *y, dvector *z);
int axy_dvector(double a, dvector *x, dvector *y);
int axy_dvector2b(double a, dvector *x, dvector *y);
double dot_dvector(dvector *x, dvector *y);
double dot_dvector2b(dvector *x, dvector *y);
double maxnorm_dvector(dvector *x);
double maxdiff_dvector(dvector *a, dvector *b);
double onenorm_dvector(dvector *x);
double twonorm_dvector(dvector *x); 

int sparse_mv0(double alpha, dCSRmat *A, double *x, double *y);
int sparse_mv(double alpha, dCSRmat *A, double *x, double *y);
int sparse_mv2b(double alpha, dCSRmat *A, dvector *x, dvector *y);
int denmat_mv(double alpha, ddenmat *A, double *x, double *y);
int dBDmat_mv(double alpha, dBDmat *A, dvector *x, dvector *y);
int dBDmat_mv0(double alpha, dBDmat *A, dvector *x, dvector *y);
int getdiag(int n, dCSRmat *A, dvector *diag);

int inverse_dBDmat(dBDmat *A, dBDmat *Ainv);
int inverse_dBDmat2(dBDmat *A, dBDmat *Ainv, int col1, int col2);


/* rref.c */
void AxBrref(ddenmat *A, ddenmat *B);
void rref1(double **A, int i, int j, int row, int col);
void rref2(double **A, int i, double alpha, int row, int col);
void rref2divide(double **A, int i, double alpha, int row, int col);
void rref3(double **A, int i, double alpha, int j, int row, int col);