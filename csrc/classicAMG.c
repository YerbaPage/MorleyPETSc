/*
 *  classicAMG.c
 *  classicAMG
 *
 *  Created by Xuehai Huang on 3/27/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file classicAMG.c
 *  \brief Ruge-Stuben AMG (main file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "header.h"

/**
	* \fn void classicAMG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param)
	* \brief Solve Ax=b by Ruge and Stuben's classic AMG;
	* \param *A pointer to the coefficient matrix
	* \param *b pointer to the dvector of right hand side
	* \param *x pointer to the dvector of dofs
	* \param *param pointer to AMG parameters
	*
	* Concrete algorithm see Appendix A.7 in book
	* U. Trottenberg, C. W. Oosterlee, and A. Sch¨uller. Multigrid. Academic Press Inc., San Diego, CA, 2001. 
	* With contributions by A. Brandt, P. Oswald and K. St¨uben.
  */
void classicAMG(dCSRmat *A, dvector *b, dvector *x, AMG_param *param)
{
	int max_levels=param->max_levels;
	int levelNum=1, i, l;
		
	clock_t classicAMG_start, classicAMG_end;

	dCSRmat AA[max_levels];
	dCSRmat P[max_levels], PT[max_levels];	
		
	for(i=0;i<max_levels;i++)
	{
		AA[i].row=0;
		AA[i].col=0;
		AA[i].IA=NULL;
		AA[i].JA=NULL;
		AA[i].val=NULL;
		
		P[i].row=0;
		P[i].col=0;
		P[i].IA=NULL;
		P[i].JA=NULL;
		P[i].val=NULL;
		
		PT[i].row=0;
		PT[i].col=0;
		PT[i].IA=NULL;
		PT[i].JA=NULL;
		PT[i].val=NULL;
	}

	AA[0].row=A->row;
	AA[0].col=A->col;
	AA[0].IA=(int*)calloc(AA[0].row+1, sizeof(int));

	for(i=0;i<=AA[0].row;i++) AA[0].IA[i]=A->IA[i];
	AA[0].JA=(int*)calloc(AA[0].IA[AA[0].row], sizeof(int));
	AA[0].val=(double*)calloc(AA[0].IA[AA[0].row], sizeof(double));
	
	for(i=0;i<AA[0].IA[AA[0].row];i++) {
		AA[0].JA[i]=A->JA[i];
		AA[0].val[i]=A->val[i];
	}
				
	classicAMG_start=clock();

	classicAMG_setup(AA, P, PT, &levelNum, param); /**< classic AMG setup */
	classicAMG_solve(AA, b, x, P, PT, levelNum, param); /**< 	classic AMG solve */
	
	classicAMG_end=clock();
	
	double classicAMGduration = (double)(classicAMG_end - classicAMG_start)/(double)(CLOCKS_PER_SEC);
	
	printf("Ruge-Stuben AMG costs %f seconds.\n", classicAMGduration);
	
	for(l=0;l<levelNum;l++)
	{
		free(AA[l].IA);
		free(AA[l].JA);
		free(AA[l].val);
	}
	
	for(l=0;l<levelNum-1;l++)
	{
		free(P[l].IA);
		free(P[l].JA);
		free(P[l].val);
		
		free(PT[l].IA);
		free(PT[l].JA);
		free(PT[l].val);
	}
}
