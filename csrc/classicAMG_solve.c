/*
 *  classicAMG_solve.c
 *  classicAMG
 *
 *  Created by Xuehai Huang on 3/27/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */
 
/*! \file classicAMG_solve.c
 *  \brief Ruge-Stuben AMG solver
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "header.h"
#include "matvec.h"

/**
	* \fn void classicAMG_solve(dCSRmat *A, dvector *b, dvector *x, dCSRmat *P, dCSRmat *PT, 
											     int levelNum, AMG_param *param)
	* \brief solve phase of Ruge and Stuben's classic AMG
	* \param *A pointer to the coefficient matrices of all levels
	* \param *b pointer to the dvector of right hand side to input
	* \param *x pointer to the dvector of dofs
	* \param *P pointer to the prolongation operators of all levels
	* \param *PT pointer to the restriction operators of all levels
	* \param levelNum integer, the number of total levels
	* \param *param pointer to AMG parameters
	*
	* Solve Ax=b using multigrid method.
	* P, PT are obtained by Ruge and Stuben's classic AMG.
	*
*/

void classicAMG_solve(dCSRmat *A, dvector *b, dvector *x, dCSRmat *P, dCSRmat *PT, 
											int levelNum, AMG_param *param)
{
	int m=b->row;
	double r[A[0].row];
	double error;
	
	int MaxIt = param->max_iter; 
	double tol = param->tol;
	int print_level = param->print_level;
	int smoother = param->smoother;
	int pre = param->presmooth_iter;
	int post = param->postsmooth_iter;
	
	clock_t solve_start, solve_end;
	solve_start=clock();
	
	int iter=0;
	double sumb=dot_array(m,b->val,b->val); // (b,b);

	while ( ++iter <= MaxIt ) // MG solver here
	{		
		multigrid(A, b, x, PT, P, 0, levelNum, smoother, pre, post, 1);
		
		// r = b-A*u
		copy_array(m,b->val,r);
		sparse_mv(-1.0,&A[0],x->val,r);
		
		if (sumb > 1e-30)
			error=sqrt(dot_array(m,r,r)/sumb);
		else
			error=sqrt(dot_array(m,r,r));
		
		if (print_level>1)						
			printf("Iteration %3d: relative residual = %e\n",iter,error);
		
		if (error<tol) break;
	}
	
	solve_end=clock();
	
	double solveduration = (double)(solve_end - solve_start)/(double)(CLOCKS_PER_SEC);
	printf("Ruge-Stuben AMG solve costs %f seconds.\n", solveduration);

	if (iter>MaxIt)
		printf("Maximal iteration %d exceeded with relative residual %e.\n", MaxIt, error);
	else
		printf("Number of iterations = %d with relative residual %e.\n", iter, error);
}
