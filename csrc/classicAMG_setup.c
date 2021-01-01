/*
 *  classicAMG_setup.c
 *  classicAMG
 *
 *  Created by Xuehai Huang on 3/27/09.
 *  Modified by Chensong Zhang on 04/04/2009.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file classicAMG_setup.c
 *  \brief Ruge-Stuben AMG setup
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "header.h"
#include "checkmat.h"

/**
 * \fn void classicAMG_setup(dCSRmat *A, dCSRmat *P, dCSRmat *PT, int *levels, AMG_param *param
 * \brief Set up phase of Ruge and Stuben's classic AMG
 * \param *A pointer to the coefficient matrices of all levels
 * \param *P pointer to the prolongation operators of all levels
 * \param *PT pointer to the restriction operators of all levels
 * \param *levels pointer to the total number of levels
 * \param *param pointer to AMG parameters
 *
 * Setup A, P, PT, levels using Ruge and Stuben's classic AMG
 * concrete algorithm see Appendix A.7 in book
 * U. Trottenberg, C. W. Oosterlee, and A. Sch¨uller. Multigrid. Academic Press Inc., San Diego, CA, 2001. 
 * With contributions by A. Brandt, P. Oswald and K. St¨uben.
 */
void classicAMG_setup(dCSRmat *A, dCSRmat *P, dCSRmat *PT, int *levels, AMG_param *param)
{
	int critiDOF=param->coarse_dof;
	int print_level=param->print_level;
	int interp_type=param->interpolation_type;
	int max_levels=param->max_levels;
	
	int i,l;
	int levelNum=*levels;
	ivector vertices[max_levels]; // each level stores the information of the vertex is on fine (current level) or coarse grid (fine: 0; coarse: 1)
	
	for(i=0;i<max_levels;i++)
	{		
		vertices[i].row=0;
		vertices[i].val=NULL;
	}
	
	clock_t setup_start, setup_end;
	setup_start=clock();
	
	while((A[levelNum-1].row>critiDOF)&& (levelNum<max_levels))
	{
		/*-- Coarseing and form the structure of interpolation --*/
		coarsening(&A[levelNum-1], &vertices[levelNum-1], &P[levelNum-1], param);
		
		/*-- Form interpolation --*/
		interpolation(&A[levelNum-1], &vertices[levelNum-1], &P[levelNum-1], param);
		
		free(vertices[levelNum-1].val);
		
		/*-- Form coarse level stiffness matrix --*/
		
		getTransposeOfSparse(&P[levelNum-1],&PT[levelNum-1]);		
		getCoarseStiffMatrix(&PT[levelNum-1],&A[levelNum-1],&P[levelNum-1],&A[levelNum]);
		
		levelNum++;
	}
	
	*levels=levelNum;
	
	free(vertices[levelNum-1].val);
	
	if (print_level>0) {
		for(l=0;l<levelNum;l++)
			printf("Level =%3d, Arow =%10d, Acol =%10d\n", l, A[l].row, A[l].col);
	}
	
	setup_end=clock();
	double setupduration = (double)(setup_end - setup_start)/(double)(CLOCKS_PER_SEC);
	printf("Ruge-Stuben AMG setup costs %f seconds.\n", setupduration);	
}
