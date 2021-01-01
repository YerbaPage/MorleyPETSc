/*
 *  precond.h
 *  Header file for preconditioners
 *  
 *------------------------------------------------------
 *
 *		Created by Chensong Zhang on 03/29/2009.
 *		Copyright 2009 PSU. All rights reserved. 
 *
 *------------------------------------------------------ 
 */

/** 
 * \file precond.h
 * \brief Header file for preconditioners
 */

#include "header.h"

/** 
 * \brief precond_data: data passed to the preconditioner for AMG.
 *
 * This is needed for the AMG preconditioner.
 */
typedef struct {

	int     print_level; /**< print level */	
	double  tol;	/**< solve params, tolerance for solver */
	int     max_iter; /**< solve params, max number of iterations */
	int     max_levels; /**< max number of levels */

	int     mass_precond_type; /**< type of preconditioning mass matrix: 1 diagonal; 2 full  */
	int     precond_type; /**< type of precondition in ASP: 1 additive; 2 multiplicative  */
	double  *precond_scale; /**< scale for the preconditioned variables in precondition  */
	int     smoother; /**< type of smoother */
	int     schwarz_type; /**< type of Schwarz smoother */
	int     smooth_iter; /**< number of smoothing */
	int     mg_smoother; /**< type of smoother */
	int     mg_smooth_iter; /**< number of smoothing in mg */
	int     presmooth_iter; /**< number of presmoothing */
	int     postsmooth_iter; /**< number of postsmoothing */
	int     coarsening_type; /**< coarsening type */

	dvector  *diag; /**< preconditioning diagonal data */
	
	dCSRmat  *precA[4]; /**< preconditioning data, the sparse matrix */
	dBDmat  *Minv; /**< preconditioning data, the dBDmat matrix */
	dOBDmat *swzB; /**< block diagonal matrix for Schwarz smoother, the dOBDmat matrix */

	dCSRmat  *As; /**< problem data, the sparse matrix */
	dCSRmat  *Rs; /**< auxiliary restriction matrix */
	dCSRmat  *Ps; /**< auxiliary prolongation matrix */
	dvector  *r;	/**< problem data, the right-hand side vector */
	dCSRmat  **Aarray; /**< data generated in the setup phase */
	dCSRmat  *R; /**< restriction matrix  */
	dCSRmat  *P; /**< prolongation matrix */

	EDGE *edges;
	iCSRmat *edgesTran;
	ivector *nodeCEdge;
	ivector *isInNode, *nondirichlet, *index;


//	edges, edgesTran, nodeCEdge, isInNode, nondirichlet, index

} precond_data;


/** 
 * \brief precond: preconditioner data and action.
 *
 * This is the preconditioner structure for PCG.
 */ 
typedef struct {
	/** data for preconditioner, void pointer */
	void *data; 
	
	/** action for preconditioner, void function pointer */
	void (*fct)(dCSRmat *, double *, double *, void *);

	/** action for preconditioner, void function pointer */
	void (*fct_dvec)(dvector *, dvector *, void *);
} precond;


/** 
 * \brief precond: preconditioner data and action.
 *
 * This is the preconditioner structure for PCG.
 */ 
typedef struct {
	/** data for preconditioner, void pointer */
	void *data; 
	
	/** action for preconditioner, void function pointer */
	void (*fct)(ddenmat *, double *, double *, void *);
} den_precond;


/* cg.c */
int pcg(dCSRmat *A, dvector *b, dvector *u, int MaxIt, double tol, precond *pre, int print_level);
int den_pcg(ddenmat *A, dvector *b, dvector *u, int MaxIt, double tol, den_precond *pre, int print_level);

/* precond.c */
void precond_null(dCSRmat *A, double *r, double *z, void *data);
void precond_diag(dCSRmat *A, double *r, double *z, void *data);
void precond_classicAMG(dCSRmat *A, double *r, double *z, void *data);
void precond_aspLaplaceVec2(dCSRmat *A, double *r, double *z, void *data);
void precond_TriAsP1Elas(dvector *r, dvector *z, void *data);
void precond_AbfpAsP1Stokes(dvector *r, dvector *z, void *data);
void precond_DiagAsP1Elas(dvector *r, dvector *z, void *data);
void precond_DiagAsP1Stokes(dvector *r, dvector *z, void *data);
void precond_DiagAsP1StokesNcP1_P0(dvector *r, dvector *z, void *data);
void precond_DiagAMGStokesNcP1_P0(dvector *r, dvector *z, void *data);

/* gmres.c */
int gmres(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level);
int gmres2b(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level);
int fgmres(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level);
int fgmres_den(ddenmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level);
int fgmres2b(dCSRmat *A, dvector *b, dvector *u, int restart, int MaxIt, double tol, precond *pre, int print_level);
int gmresBCD(dCSRmat *B, dCSRmat *C, dCSRmat *D, dvector *b, dvector *u, int restart, int MaxIt, double tol, int print_level);

/* minres.c */
int minres(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, precond *pre, int print_level);
int minres2b(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, precond *pre, int print_level);
