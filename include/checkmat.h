/**
 *		checkmat.h 
 *
 *		Matrix-Solver Community Project
 *
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 04/02/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------
 *
 */

/*! \file checkmat.h
 *  \brief Header file for property-check subroutines
 */

/* checkmat.c */
int check_diagpos(dCSRmat *A);
void transpose(int *row[2], int *col[2], double *val[2], int *nn, int *tniz);
int transpspm(dCSRmat *A, dCSRmat *B);
int check_symm(dCSRmat *A);
int check_diagdom(dCSRmat *A);
