/*
 *  header.h
 *  MGV
 *
 *  Created by Xuehai Huang on 03/27/2009.
 *  Modified by Chensong Zhang on 03/29/2009.
 *
 *  Copyright 2008 PSU. All rights reserved.
 *
 */

/*! \file header.h
 *  \brief main header file for the MSC package
 */

/** 
 * \brief Definition of max, min, abs
 */
#define max(a,b) (((a)>(b))?(a):(b)) /**< bigger one in a and b */
#define min(a,b) (((a)<(b))?(a):(b)) /**< smaller one in a and b */
#define abs(a) (((a)>0.0)?(a):-(a)) /**< absolute value of a */

#define PI 3.1415926535897932384626
#define pi 3.1415926535897932384626

#define DBL_DIG         15                      /* # of decimal digits of precision */
#define DBL_EPSILON     2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
#define DBL_MANT_DIG    53                      /* # of bits in mantissa */
#define DBL_MAX         1.7976931348623158e+308 /* max value */
#define DBL_MAX_10_EXP  308                     /* max decimal exponent */
#define DBL_MAX_EXP     1024                    /* max binary exponent */
#define DBL_MIN         2.2250738585072014e-308 /* min positive value */
#define DBL_MIN_10_EXP  (-307)                  /* min decimal exponent */
#define DBL_MIN_EXP     (-1021)                 /* min binary exponent */
#define _DBL_RADIX      2                       /* exponent radix */
#define _DBL_ROUNDS     1                       /* addition rounding: near */


/** 
 * \brief Smoother type
 */
#define JACOBI 1  /**< Jacobi smoother */
#define GS     2  /**< Gauss-Seidel smoother */
#define SGS    3  /**< symm Gauss-Seidel smoother */
#define MSWZ   4  /**< multiplicative Schwarz smoother */
#define ASWZ   5  /**< additive Schwarz smoother */
#define SMSWZ   6  /**< symm multiplicative Schwarz smoother */

/**
 * \brief input parameters 
 *
 * Input parameters, reading from disk file
 */
typedef struct {
	
	// problem parameters	
	int problem_num;	/**< problem number to solve */
	
	int domain_num; /**< domain number */
					
	// parameters for iterative solvers
	int itsolver_type; /**< type of iterative solvers */
	double itsolver_tol; /**< tolerance for iterative linear solver */
	int itsolver_maxit; /**< maximal number of iterations for iterative solvers */
	int restart; /**< restart number used in GMRES */

	int precond_type; /**< type of precondition in ASP: 1 additive; 2 multiplicative  */
	double precond_scale[2]; /**< scale for the preconditioned variables in precondition */
	int smoother; /**< type of smoother in ASP */
	int schwarz_type; /**< type of Schwarz smoother */
	int smooth_iter; /**< number of smoothing in ASP */

	// parameters for MG
	double MG_tol; /**< tolerance for MG if used as preconditioner */
	int MG_maxit; /**< max number of iterations for MG if used as preconditioner */
	int MG_smoother; /**< type of smoother in MG */
	int MG_smooth_iter; /**< number of smoothing in MG */
//	int MG_presmooth_iter; /**< number of presmoothing */
//	int MG_postsmooth_iter; /**< number of postsmoothing */
	int AMG_levels; /**< maximal number of levels */
	int AMG_coarsening_type; /**< coarsening type */
	int AMG_interpolation_type; /**< interpolation type */
	int AMG_coarse_dof;	/**< minimal coarsest level dof */
	double AMG_strong_threshold; /**< strong threshold for coarsening */
	double AMG_truncation_threshold; /**< truncation factor for interpolation */
	double AMG_max_row_sum; /**< maximal row sum */

	// parameters for AFEM
	double AFEM_tol; /**< tolerance for AFEM */
	int AFEM_maxit; /**< max number of iterations for AFEM */
	double AFEM_mark_threshold; /**< mark threshold for AFEM */
	
	// output flags
	int print_level; /**< print level */	
	
	// temp
	double nu;
	double lambda;
	double mu;
	double t;
	double alpha1;
	double alpha2;
	double alpha3;
	double beta1;
	double beta2;
	double beta3;
	int glevelNum;
	int CglevelNum;
	int FglevelNum;
	int cDOP;
	int rDOP;
	int dop1; /**< degree of polynomial for stress tensor */
	int dop2; /**< degree of polynomial for displacement */
	int dop3; /**< degree of polynomial for trace of normal derivative of displacement */
	int dop4; /**< degree of polynomial for trace of displacement */

	int variationalform_type; /**< type of variational form */
	int stress_fem_type; /**< type of fem space for stress */

} Input_data;

/** 
 * \brief AMG_param: parameters for AMG solver.
 *
 * This is needed for the AMG solver.
 *
 * What coarsening type available? 
 * What are the parameters in detail?
 */
typedef struct {

	int print_level;								/**< print level for AMG */
	int max_levels;								/**< setup param, max number of levels */
	int coarse_dof;								/**< minimal coarsest level dof */
	int max_iter;									/**< solve params, max number of iterations */
	double tol;										/**< solve params, tolerance for solver */
	int AMG_max_iter;							/**< AMG params, max number of iterations */
	double AMG_tol;								/**< AMG params, tolerance for solver */
	int restart;                                /**< restart number used in GMRES */

	int smoother;									/**< smoother type */
	int presmooth_iter;						/**< number of presmoothers */
	int postsmooth_iter;						/**< number of postsmoothers */
	
	int coarsening_type;						/**< coarsening type */
	int interpolation_type;			  /**< interpolation type */	
	double strong_threshold;				/**< strong connection threshold for coarsening */
	double max_row_sum;						/**< maximal row sum parameter */
	double truncation_threshold;		/**< truncation threshold */

} AMG_param;

/** 
 * \brief Dense matrix of double type.
 *
 * A dense double matrix
 */ 
typedef struct ddenmat{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	double **val;
}ddenmat;

/** 
 * \brief Dense matrix of int type.
 *
 * A dense int matrix
 */ 
typedef struct idenmat{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	int **val;
}idenmat;

/** 
 * \brief ELEMENT struct.
 *
 * ELEMENT struct
 */ 
typedef struct ELEMENT{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** the first 3 columns store the indexes of vertices,
	    the last 3 columns store the indexes of edge's midpoints */
	int **val;
	/** parent element of current element, which will be the root if its parent = -1*/
	int *parent;
	/** volume of the element*/
	double *vol;
	double **eta; 
	double **xi;
	/** unit normal vector of three edges*/
	double ***nvector;
	/** unit tangential vector of three edges*/
	double ***tvector;
	/** length of three edges*/
	double **edgeslength;
}ELEMENT;

/** 
 * \brief EDGE struct.
 *
 * EDGE struct
 */ 
typedef struct EDGE{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** the first two columns store the two vertice, the third and fourth columns store the affiliated elements£¬
        the fourth column stores -1 if the edge is on boundary */
	int **val;
	/** sizes of the matrix entries*/
	double *h;
	/** x1-x2*/
	double *xi;
	/** y1-y2*/
	double *eta;
	/** unit normal vector of edge*/
	double **nvector;
	/** unit tangential vector of edge*/
	double **tvector;
	/** length of edge*/
	double *length;
	/** boundary type of edge
	    0 : non - boundary, i.e., an interior edge or face.
		1 : first type, i.e., a Dirichlet boundary edge or face.
		2 : second type, i.e., a Neumann boundary edge or face.
		3 : third type, i.e., a Robin boundary edge or face.
    **/
	int *bdFlag;
}EDGE;

/** 
 * \brief Sparse matrix of double type in CSR format.
 *
 * CSR Format (IA,JA,A)
 *
 * The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct dCSRmat{
	//! row number of matrix A, m
	int row;   
	//! column of matrix A, n
	int col;   
	//! number of nonzero entries
	int nnz;
	//! integer array of row pointers, the size is m+1
	int *IA;   
	//! integer array of column indexes, the size is nnz
	int *JA;    
	//! nonzero entries of A
	double *val;
}dCSRmat;

/** 
 * \brief Sparse matrix of int type in CSR format.
 *
 * CSR Format (IA,JA,A)
 *
 * The starting index of A is 0, other data stuctures also take this convention.  
 */
typedef struct iCSRmat{
	//! row number of matrix A, m
	int row;   
	//! column of matrix A, n
	int col;   
	//! number of nonzero entries
	int nnz;
	//! integer array of row pointers, the size is m+1
	int *IA;   
	//! integer array of column indexes, the size is nnz
	int *JA;    
	//! nonzero entries of A
	int *val;
}iCSRmat;

/**
* \brief Vector with n entries of double type.
*/
typedef struct dvector {
	//! number of rows
	int row;
	//! actual vector entries
	double *val;
}dvector;

/**
* \brief Vector with n entries of int type.
*/
typedef struct ivector {
	//! number of rows
	int row;
	//! actual vector entries
	int *val;
}ivector;

/** 
 * \brief Block diagonal matrix of double type.
 *
 * A dense double block diagonal matrix
 */ 
typedef struct dBDmat{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** number of blocks */
	int nb;	
	/** blocks */
	ddenmat *blk;
}dBDmat;

/**
* \brief Overlapped block diagonal matrix of double type.
*
* A dense overlapped double block diagonal matrix
*/
typedef struct dOBDmat {
	/** number of rows */
	int row;
	/** number of columns */
	int col;
	/** number of blocks */
	int nb;
	/** blocks */
	ddenmat *blk;
	/** indices to local rows */
	ivector *rindices;
	/** indices to local columns */
	ivector *cindices;
}dOBDmat;

/** 
 * \brief Sparse matrix of double type in IJ format.
 *
 * Coordinate Format (I,J,A)
 *
 */
typedef struct dIJmat{
	//! row number of matrix A, m
	int row;   
	//! column of matrix A, n
	int col;   
	//! number of nonzero entries
	int nnz;
	//! integer array of row indices, the size is nnz
	int *I;   
	//! integer array of column indices, the size is nnz
	int *J;    
	//! nonzero entries of A
	double *val;
}dIJmat;

/** 
 * \brief Dennode type.
 */ 
typedef struct dennode{
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	double **val;
	/** boundary type of vetex
	    0 : non - boundary, i.e., an interior vertex.
	    1 : first type, i.e., a Dirichlet boundary vertex.
	    2 : second type, i.e., a Neumann boundary vertex.
	    3 : third type, i.e., a Robin boundary vertex.
		12: a Dirichlet-Neumann boundary vertex.
		22: a Neumann-Neumann boundary vertex.
	**/
	int *bdFlag;
	/** tensor basis of nodes for Hu-Zhang element in 2d */
	double ***tensorBasis;
}dennode;

/** 
 * \brief ELEMENT_DOF type.
 */ 
typedef struct ELEMENT_DOF{
	/** number of degree of polynomial */
	int dop;	  
	/** number of global degree of freedom */
	int dof;	  
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	int **val;
	/** non-free flag
	1 : non-free node
	0 : free node
	**/
	ivector nfFlag;
	/** free nodes */
	ivector freenodes;
	/** non-free nodes */
	ivector nfreenodes;
	/** index of free and non-free nodes */
	ivector index;
}ELEMENT_DOF;

/**
* \brief AS_param: parameters for auxiliary space preconditioner.
*
* This is needed for the AS preconditioner.
*
* What coarsening type available?
* What are the parameters in detail?
*/
typedef struct {

	int problem_num;
	int print_level;								/**< print level for ASP */
	int max_iter;									/**< solve params, max number of iterations */
	double tol;										/**< solve params, tolerance for solver */
	int restart;                                /**< restart number used in GMRES */

	int mass_precond_type; /**< type of preconditioning mass matrix: 1 diagonal; 2 full  */
	int precond_type; /**< type of precondition in ASP: 1 additive; 2 multiplicative  */
	double *precond_scale; /**< scale for the preconditioned variables in precondition  */
	int smoother;									/**< smoother type */
	int schwarz_type;                       /**< type of Schwarz smoother */
	int smooth_iter;						/**< number of smoothers */

	int levelNum;								/**< setup param, max number of levels */
	int mg_max_iter;							/**< MG params, max number of iterations */
	double mg_tol;								/**< MG params, tolerance for solver */
	int mg_smoother;									/**< smoother type */
	int mg_smooth_iter;						/**< number of smoothers */

											/************* Geometric information *************/
	ELEMENT *elements;
	idenmat *elementEdge;
	EDGE *edges;
	dennode *nodes;
	iCSRmat *edgesTran;
	ivector *nodeCEdge;
	ELEMENT_DOF *elementDOF;
	iCSRmat *elementdofTran;

	/************* Material information *************/
	double lambda; /**< Lame constant */
	double mu; 
	double t;
	double nu;

	int variationalform_type; /**< type of variational form */
	int stress_fem_type; /**< type of fem space for stress */

	int max_levels;								/**< setup param, max number of levels */
	int AMG_levels; /**< maximal number of levels */
	int AMG_coarsening_type; /**< coarsening type */
	int AMG_interpolation_type; /**< interpolation type */
	int AMG_coarse_dof;	/**< minimal coarsest level dof */
	double AMG_strong_threshold; /**< strong threshold for coarsening */
	double AMG_truncation_threshold; /**< truncation factor for interpolation */
	double AMG_max_row_sum; /**< maximal row sum */
} ASP_param;

/**
* \brief AFEM_param: parameters for adaptive finite element method.
*
* This is needed for adaptive finite element method.
*
*/
typedef struct {

	AMG_param *amgparam;                             /**< AMG params */
	ASP_param *aspparam;                             /**< AMG params */
	double tol;										/**< tolerance for AFEM */
	int max_iter;							        /**< max number of iterations for AFEM */
	int solver_type;   					        /**< type of solver */
	double mark_threshold;			        	/**< marking threshold for AFEM */
} AFEM_param;

/** 
 * \brief Dense three-dimensional matrix of double type.
 *
 * A dense three-dimensional double matrix
 */ 
typedef struct ddenmat3{
	/** number of pages */
	int pag;	  
	/** number of rows */
	int row;	  
	/** number of columns */
	int col;	
	/** actual matrix entries */
	double ***val;
}ddenmat3;


/** 
 * \brief Links
 */
struct linked_list
{
  //! data
	int                        data;
	//! next element
	struct linked_list *next_elt;
	//! previous element
	struct linked_list *prev_elt;
	//! starting position
	int                        head;
	//! ending position
	int                        tail;
};
typedef struct linked_list ListElement;

/** 
 * \brief List of links
 */
typedef ListElement  *LinkList;  



/*--- Function declaration ---*/

/* input.c */
int read_Input_data(char *filenm, Input_data *Input);

/* coarsening.c */
void coarsening(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param);
void generateS(dCSRmat *A, iCSRmat *S, AMG_param *param);
void generateSRS(dCSRmat *A, iCSRmat *S, double epsilon_str, int coarsening_type);
void getCoarseStiffMatrix(dCSRmat *PT, dCSRmat *A, dCSRmat *P, dCSRmat *Ac);

/* classicAMG.c */
void classicAMG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param);

/* classicAMG_setup.c */
void classicAMG_setup(dCSRmat *A, dCSRmat *P, dCSRmat *PT, int *levels, AMG_param *param);

/* classicAMG_solve.c */
void classicAMG_solve(dCSRmat *A, dvector *b, dvector *x, dCSRmat *P, dCSRmat *PT, int levelNum, AMG_param *param);

/* interpolation.c */
void interpP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg);
void interpP1toP2_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFp2, EDGE *edges);
void interpVecP1toDG2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFdg);
void interpVecP1toNcP1_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges);
void interpVecP1toMINI_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFmini);
void interpVecP1toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp1, ELEMENT_DOF *elementDOFcr, EDGE *edges);
void interpVecP2toCR_2d(dCSRmat *P, ELEMENT_DOF *elementDOFp2, ELEMENT_DOF *elementDOFcr);
void interpolation(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param);
void interpolationRS(dCSRmat *A, ivector *vertices, dCSRmat *Ptr, AMG_param *param);
int getiteval(dCSRmat *A, dCSRmat *it);

/* smoother.c */
void jacobi(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L);
void gs(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, int L);
void mulschwarz(dvector *u, int i_1, int i_n, int s, dCSRmat *A, dvector *b, dOBDmat *B, int L);
void getSchwarzblocks_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocks_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_vertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_edge(dOBDmat *B, dCSRmat *A, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_edgevertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, EDGE *edges, int nvertices, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_element(dOBDmat *B, dCSRmat *A, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF);
void getSchwarzblocksVec2_elementvertex(dOBDmat *B, dCSRmat *A, ELEMENT *elements, int nvertices, ELEMENT_DOF *elementDOF);

/* multigrid.c */
void multigrid(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P, int level, int levelNum, int smoother, int m1, int m2, int m0);
void multigridvar(dCSRmat *A, dvector *b, dvector *x, dCSRmat *R, dCSRmat *P, int level, int levelNum, int smoother, int *m1, int *m2, int m0);
double mgvVectorP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m);
double mgvP1_solve(dCSRmat *A, dvector *b, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *isInNode, ivector *nondirichlet, ivector *index, int levelNum, int smoother, int m);
void interpolationPvector2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe);
void interpolationP2d(EDGE *Cedges, ivector *CisInNode, ivector *Cindex, ivector *Fnondirichlet, ivector *nodeCEdge, int cnn, dvector *e, dvector *Pe);
void restrictionPTvector2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr);
void restrictionPT2d(iCSRmat *FedgesTran, ivector *Findex, ivector *Cnondirichlet, dvector *r, dvector *PTr);
double mgvVectorP1_solveOld(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, dvector *x, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge, ivector *nondirichlet, ivector *index, int levelNum, int m);
void interpolationPvector2dOld(EDGE *Cedges, ivector *nodeCEdge, dvector *e, dvector *Pe);
void restrictionPTvector2dOld(iCSRmat *edgesTran, dvector *r, dvector *PTr);

/* itsolver.c */
int standard_CG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level);
int diag_PCG(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, int print_level);
int classicAMG_PCG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level);
int DiagAsP1StokesCR_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1StokesMINI_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1StokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAMGStokesNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int AbfpAsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int Abfp2AsP1StokesNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1ElasCR_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1ElasNcP1P0_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int DiagAsP1ElasDG_MINRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int TriAsElasMINI_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type, int print_level);
int TriAsP1ElasMINI_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int TriAsElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type, int print_level);
int TriAsP1ElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int TriAsP2ElasCR_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int TriAsP1ElasNcP1P0_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int TriAsP1ElasDG_GMRES(dCSRmat *A, dvector *b, dvector *x, ASP_param *param, int print_level);
int classicAMG_GMRES(dCSRmat *A, dvector *b, dvector *x, AMG_param *param, int print_level);
int cg_pcg(dCSRmat *A, dvector *b, dvector *sigma, dvector *u, int MaxIt, double tol, AMG_param *param, int itsolver_type, int print_level);

/* asElasDG.c */
double mgvVectorP1As_solve(dCSRmat *A, dvector *b, dvector *x, void *data);

/* linklist.c */
void dispose_elt ( LinkList element_ptr );
void remove_point ( LinkList *LoL_head_ptr , LinkList *LoL_tail_ptr , int measure , int index , int *lists , int *where );
LinkList create_elt ( int Item );
void enter_on_lists ( LinkList *LoL_head_ptr , LinkList *LoL_tail_ptr , int measure , int index , int *lists , int *where );

/* sparse.c */
int get_block(dCSRmat *A, int m, int n, int *rows, int *cols, double *Aloc, int *mask);
int get_block_dden(dCSRmat *A, int m, int n, int *rows, int *cols, ddenmat *Aloc, int *mask);
iCSRmat getCSRfromIntT(idenmat *T, int nNode);
iCSRmat getTransposeOfA(iCSRmat *A);
void getTransposeOfSparse(dCSRmat *A,dCSRmat *AT);
void getTransposeOfiden(idenmat *A, iCSRmat *AT, int cols, int nNode);
void getTransposeOfELEMENT(ELEMENT *A, iCSRmat *AT, int cols, int nNode);
void getTransposeOfelementDoF(ELEMENT_DOF *A, iCSRmat *AT, int cols);
int sparseAddition(dCSRmat *A, dCSRmat *B, dCSRmat *C);
int sparseSubtraction(dCSRmat *A, dCSRmat *B, dCSRmat *C);
int dCSRPlusdDiagVector(dCSRmat *A, dvector *B, dCSRmat *C);
int dCSRMultiplydDiagVector(dCSRmat *A, dvector *D, dCSRmat *C);
int dDiagVectorMultiplydCSR(dvector *D, dCSRmat *A, dCSRmat *C);
int dCSRPlusdBD(dCSRmat *A, dBDmat *B, dCSRmat *C);
int dBDMinusdCSR(dBDmat *A, dCSRmat *B, dCSRmat *C);
void sparseMultiplication(dCSRmat *A, dCSRmat *B, dCSRmat *C);
void dBDMultiplydCSR(dBDmat *A, dCSRmat *B, dCSRmat *C);
void dBD2MultiplydCSR(dBDmat *A, int srow, int scol, dCSRmat *B, dCSRmat *C);
void dCSRMultiplydBD(dCSRmat *A, dBDmat *B, dCSRmat *C);
void sparseTripleMultiplication(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication1(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication2(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
void sparseTripleMultiplication3(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *B);
int dIJtoCSR(dIJmat *A, dCSRmat *B);
int dCSRtoIJ(dCSRmat *A, dIJmat *B);

/* kirchhoffplatefem.c */
void kirchhoffplatefem(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void kirchhoffplateMorley(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void kirchhoffplatefem_psp(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void kirchhoffplatefem_psprot(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void kirchhoffplatefem_pep(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void kirchhoffplateMorley_psp(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input);
void kirchhoffplateMorley_psprot(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input);
void kirchhoffplateMorley_pep(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOFm, Input_data *Input);

/* fem.c */
void stokesfem(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesfem_elas(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesfem_psv(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesNcP1P0(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesNcP1P0_elas(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesMINI(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesMINI_psv_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesMINI_psv_stressP1P0(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesCR(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesCR_elas_stressCompact(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesCR_elas_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesCR_psv_stressCompact(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);
void stokesCR_psv_stressP2(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, ELEMENT_DOF *elementDOF, Input_data *Input);

/* assemble.c */
int getmesh(int domain_num, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum);
void getElementDOF1D(ELEMENT_DOF *elementDOF, int ne, int dop);
void getElementDOF1D_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF(ELEMENT_DOF *elementDOF, int nt, int dop);
void getElementDOF_symTensor(ELEMENT_DOF *elementDOF, int nt, int dop);
void getElementDOF_Lagrange(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void getElementDOF_NoncfmP1(ELEMENT_DOF *elementDOF, idenmat *elementEdge, int ne);
void getElementDOF_Morley(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int ne, int nvertices);
void getElementDOF_MINI(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices);
void getElementDOF_MINIsymTensor(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices);
void getElementDOF_MINIsymTensorP1P0(ELEMENT_DOF *elementDOF, int nt);
void getElementDOF_CR(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int ne, int nvertices);
void getElementDOF_CRsymTensor(ELEMENT_DOF *elementDOF, int nt);
void getElementDOF_HuZhang(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop);
void assemble(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);

void assemble_kirchhoffplateMorley(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double nu);
void assembleStiffmatrixKirchhoffplateMorley2d(dCSRmat *A, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double nu);

void assemble_poissonMorley(dCSRmat *ptr_A1, dvector *ptr_A2, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assembleStiffmatrixPoissonMorley2d(dCSRmat *A1, dvector *A2, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran);
void assemble_kirchhoffplate_StokesNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu);
void assembleStiffmatrix_biharmonic_StokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm);
void assembleStiffmatrix_kirchhoffplate_StokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu);
void assemble_kirchhoffplate_StokesrotNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu);
void assembleStiffmatrix_biharmonic_StokesrotNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm);
void assembleStiffmatrix_kirchhoffplate_StokesrotNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu);
void assemble_kirchhoffplate_ElasNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu);
void assembleStiffmatrix_biharmonic_ElasNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm);
void assembleStiffmatrix_kirchhoffplate_ElasNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu);
void assembleRHScurl_biharmonic_PoissonMorley2d(dvector *ptr_b, dvector *uhv, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFv);
void assembleRHSgrad_biharmonic_PoissonMorley2d(dvector *ptr_b, dvector *uhv, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFv);
void assembleStressMassmatrixKirchhoffplateMorley2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu, double nu);
void assembleStressMassmatrixKirchhoffplateMorleydBD2d(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu, double nu);

void assemble_stokesCR(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrixStokesCR2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assemblePressMassmatrixStokesCR2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double t);
void assemblePressMassmatrixStokesCRdBD2d(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double t);

void assemble_stokesMINI(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrixStokesMINI2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assemblePressMassmatrixStokesMINI2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double t);

void assemble_stokesNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrixStokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assemblePressMassmatrixStokesNcP1P02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF);

void assemble_stokesCR_elas_stressCompact(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrix_stokesCR_elas_stressCompact(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStressMassmatrixStokesCR2d_stressCompact(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);
void assembleStressMassmatrixStokesCRdBD2d_stressCompact(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);

void assemble_stokesCR_elas_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrix_stokesCR_elas_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStressMassmatrixStokesCR2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);
void assembleStressMassmatrixStokesCRdBD2d_stressP2(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);

void assemble_stokesCR_psv_stressCompact(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrix_stokesCR_psv_stressCompact(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assemblePressStressMassmatrixStokesCR2d_stressCompact(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);
void assemblePressStressMassmatrixStokesCRdBD2d_stressCompact(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);

void assemble_stokesCR_psv_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrix_stokesCR_psv_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assemblePressStressMassmatrixStokesCR2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);
void assemblePressStressMassmatrixStokesCRdBD2d_stressP2(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);

void assemble_stokesMINI_psv_stressP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrix_stokesMINI_psv_stressP1P0(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assemblePressStressMassmatrixStokesMINI2d_stressP1P0(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu);
void assemblePressStressMassmatrixStokesMINIdBD2d_stressP1P0(dBDmat *M, dvector *diag, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);

void assemble_stokesMINI_psv_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrix_stokesMINI_psv_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assemblePressStressMassmatrixStokesMINI2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu);
void assemblePressStressMassmatrixStokesMINIdBD2d_stressP2(dBDmat *M, dvector *diag, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);

void assemble_stokesNcP1P0_elas(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStiffmatrix_stokesNcP1P0_elas(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t);
void assembleStressMassmatrixStokesNcP1P02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu);

void assembleStiffmatrixHuZhang2d(dCSRmat *A, dCSRmat *B, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu);
void assembleStiffmatrixHuZhangA11_2d(dCSRmat *A, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu);
void assembleStiffmatrixHuZhangIP2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu);
void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu);
void assembleStiffmatrixVecLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu);
void assembleStiffmatrixLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu);
void assembleStiffmatrixElasIPDG0(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void assembleStiffmatrixVecIPDG0(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void sumNormalDerivativeEta(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, ddenmat *etas, double *sum);
void jumpOperatorVectorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorVector(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void getinfo4jump(double lambda1, double lambda2, int *edgeNode, double elen, int element, ELEMENT *elements, int dop, int index, double *jump);
void TangentDerivative4Edge(double lambda1, double lambda2, int edge, int element, ELEMENT *elements, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *val);
void averageOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *average);
void getinfo4average(double lambda1, double lambda2, int *edgeNode, int element, ELEMENT *elements, int dop, int index, double *average);
void jumpOperatorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorNormal(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorTangent(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorDIVPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorTensorDIV(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump);
void jumpOperatorATensorTangent2(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump);
void jumpOperatorATensorTangent2Dirichlet(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump);
void jumpOperatorRotATensorTangentPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump);
void NeumannStress(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump);
void getElementEdgeGeoInfo(ELEMENT *elements, EDGE *edges, dennode *nodes);
void getFreenodesInfoHuZhang(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoLagrange(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoNcP1(EDGE *edges, ELEMENT_DOF *elementDOF);
void getFreenodesInfoMorley(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoNcP1Vector2(EDGE *edges, ELEMENT_DOF *elementDOF);
void getFreenodesInfoMINIVector2(dennode *nodes, ELEMENT_DOF *elementDOF);
void getFreenodesInfoCRVector2(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void getBoundaryInfoNormalTrace(EDGE *edges, ELEMENT_DOF *elementDOF, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void getBoundaryInfoEdge(EDGE *edges, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void getTensorBasis4HuZhang2d(EDGE *edges, iCSRmat *edgesTran, dennode *nodes);
int getlocaldofsEdge4HuZhang2d(int *patchnodes, int edgeindex, int dop);
void NeumannbdValHuZhang(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);


/* assembleInfo.c */
void getEdgeInfo(ELEMENT *elements, iCSRmat *elementsTran, EDGE *edges, iCSRmat *edgesTran, dennode *nodes);
void getEdgeBdflagFromVertexBdflag(EDGE *edges, dennode *nodes);
int getCoarseInfo(int domain_num, dennode *nodes, ELEMENT *elements, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge);
void uniformrefine(dennode *Cnodes, ELEMENT *Celements, EDGE *Cedges, dennode *Fnodes, ELEMENT *Felements, EDGE *Fedges, iCSRmat *FedgesTran, ivector *nodeCEdge);
void extractNondirichletMatrixVector(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index, dvector *uh);
void extractNondirichletMatrix11(dCSRmat *A, dCSRmat *A11, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index);
void extractNondirichletMatrix1r(dCSRmat *A, dCSRmat *A1, ivector *dirichlet, ivector *nondirichlet);
void extractNondirichletMatrix1c(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index);
void extractNondirichletMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index);
void extractNondirichletVector(dCSRmat *A, dvector *b, dvector *b1, ivector *dirichlet, ivector *nondirichlet, dvector *uh);
void extractFreenodesVector2StokesDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh);
void extractFreenodesVector2ElasDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh);
void extractFreenodesVector(dCSRmat *A, dvector *b, dvector *b1, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, dvector *uh);
void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
void extractFreenodesMatrix1r(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix1c(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF);
void extractFreenodesMatrix11cBlock(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc);
int getEdgeDOFsTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32);
int getEdgeDOFsScalarTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int rowstart, int row31, int row32);
int getEdgeDOFsVectorTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32);
int getEdgeDOFsVector(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row21);
void patchtogether22(dCSRmat *A11, dCSRmat *A12, dCSRmat *A21, dCSRmat *A22, dCSRmat *A);
void patchtogether33(dCSRmat *A11, dCSRmat *A12, dCSRmat *A13, dCSRmat *A21, dCSRmat *A22, dCSRmat *A23, dCSRmat *A31, dCSRmat *A32, dCSRmat *A33, dCSRmat *A);
void getPenaltyParameters(ddenmat *etas, idenmat *elementEdge, EDGE *edges, double *alpha, double *beta);

/* basicData.c */
void cart2pol(double x, double y, double *r, double *theta);
double f(double x, double y);
double u(double x, double y);
double u_x(double x, double y);
double u_y(double x, double y);
double u_xx(double x, double y);
double u_xy(double x, double y);
double u_yy(double x, double y);
double f1(double x, double y, double lambda, double mu);
double f2(double x, double y, double lambda, double mu);
double u1(double x, double y, double lambda, double mu);
double u2(double x, double y, double lambda, double mu);
double u1_x(double x, double y, double lambda, double mu);
double u1_y(double x, double y, double lambda, double mu);
double u2_x(double x, double y, double lambda, double mu);
double u2_y(double x, double y, double lambda, double mu);
double u1_xx(double x, double y, double lambda, double mu);
double u1_xy(double x, double y, double lambda, double mu);
double u1_yx(double x, double y, double lambda, double mu);
double u1_yy(double x, double y, double lambda, double mu);
double u2_xx(double x, double y, double lambda, double mu);
double u2_xy(double x, double y, double lambda, double mu);
double u2_yx(double x, double y, double lambda, double mu);
double u2_yy(double x, double y, double lambda, double mu);
double sigma11(double x, double y, double lambda, double mu);
double sigma22(double x, double y, double lambda, double mu);
double sigma12(double x, double y, double lambda, double mu);
double sigma21(double x, double y, double lambda, double mu);
double ut_t(double x, double y, double lambda, double mu, int edgeFlag);
double ut_t_1(double x, double y, double lambda, double mu);
double ut_t_2(double x, double y, double lambda, double mu);
double ut_t_3(double x, double y, double lambda, double mu);
double ut_t_4(double x, double y, double lambda, double mu);
double un_tt(double x, double y, double lambda, double mu, int edgeFlag);
double un_tt_1(double x, double y, double lambda, double mu);
double un_tt_2(double x, double y, double lambda, double mu);
double un_tt_3(double x, double y, double lambda, double mu);
double un_tt_4(double x, double y, double lambda, double mu);
void morley_basis(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double *phi);
void morley_basis1(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double phi[2]);
void morley_basis2(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double phi[3]);
void lagrange1D_basis(double lambda, int index, int dop, double *phi);
void lagrange1D_basis1(double lambda, int index, int dop, double h, double *phi);
void lagrange_basis(double *lambda, int index, int dop, double *phi);
void lagrange_basis1(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[2]);
void lagrange_basis2(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[3]);
void lagrangeSymTensor_basis(double *lambda, int index, int dop, double phi[3]);
void ncp1_basis(double *lambda, int index, double *phi);
void ncp1_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2]);
void mini_basis(double *lambda, int index, double *phi);
void mini_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2]);
void miniSymTensor_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3]);
void miniSymTensorP1P0_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3]);
void cr_basis(double *lambda, int index, double *phi);
void cr_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2]);
void crSymTensor_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3]);
void rt_basis(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[2]);
void rt_basis1(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[4]);
void arnoldwinther_basis(double *lambda, double *x, double *y, ddenmat3 *basisCoeffs, int element, int index, double *phi);
void arnoldwinther_basisDIV(double *lambda, ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi);
void arnoldwinther_basisDIV2(ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi);
void huzhang_basis(double *lambda, double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[3]);
void huzhang_basis1(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double(*phi)[2]);
void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2]);
void huzhang_basisROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2]);
void huzhang_basisCurlTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2]);
void huzhang_basisROTROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double *phi);
void huzhang_basisLaplaceTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double *phi);
double area(double x1,double x2,double x3,double y1,double y2,double y3);
void localb(double (*nodes)[2],double *b);

/* basiscoeff.c */
void generateBasisCoeffs(ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes);
void generateBasisCoeffsEdgeProj(dBDmat *B, dBDmat *bC, EDGE *edges, int dopl, int dopk);

/* quadrature.c */
int getNumQuadPoints(int dop, int dim);
void init_quadrature(int num_qp, int ncoor, double (*gauss)[3]);
void init_Gauss(int num_qp, int ncoor, double (*gauss)[3]);
void init_Gauss1D(int num_qp, int ncoor, double (*gauss)[2]);
void init_NewtonCotes1D(int num_qp, int ncoor, double (*newtoncotes)[2]);

/* lu.c */
int LU_Decomp(double *A, int pivot[], int n);
int LU_Solve(double *A, double b[], int pivot[], double x[], int n);

/* post_processing.c */
void projPiecewiseLagrangeDisplacement(dvector *Qhu, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);
void projPiecewiseLagrangeRHS(dvector *Qhf, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);
void postprocess2newDisplacement(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);

/* error.c */
void geterrors(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);
void getposteriorierrors(double *errors, dvector *sigmah, dvector *uh, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu);
void geterrors_kirchhoffpalte_morley(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);
void geterrors2discr_kirchhoffpalte_morley(double *errors, dvector *uh1, dvector *uh2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF);

/* adaptiveFEM.c */
double adaptiveFEM(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, iCSRmat *edgesTran, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uh, dvector *uhstar, double lambda, double mu, AFEM_param *afemparam);
int solve(dCSRmat *A, dvector *b, dvector *x, ASP_param *aspparam, int itsolver_type);

/* estimate.c */
double estimate(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *sigmah, dvector *uh, dvector *uhstar, double lambda, double mu, dvector *estimator);

/* refine.c */
double mark(idenmat *elementEdge, EDGE *edges, dennode *nodes, dvector *estimator, int *index, double theta, double totalError, int *marker);
void refine(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, iCSRmat *edgesTran, dennode *nodes, int *marker);

/* sort.c */
void quicksort(double *a, int low, int high);
void quicksortIndex(double *a, int *index, int low, int high);

/* output.c */
int getElementsNodes4Matlab(ELEMENT *elements, dennode *nodes, dvector *uh);


/* xuludmil.for (Xiangtan energy minimization code) */
//void get_p_xuludmil_(int *ia,int *ja,double *a, int *n, int *nc,int *ip,int *jp,double *pn,int *ipt,int *jpt,int *mf);
