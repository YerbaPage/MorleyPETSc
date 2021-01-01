/*
 *  coarsening.c
 *  classicAMG
 *
 *  Created by Xuehai Huang on 1/28/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file coarsening.c
 *  \brief Functions for coarsening 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "header.h"
#include "matvec.h"

/**
 * \fn void coarsening(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param)
 * \brief standard coarsening
 * \param *A pointer to the coefficient matrix, the index starts from zero
 * \param *vertices pointer to the indicator ivector of the CF splitting of the vertices, 0: fine (current level) gird or 1: coarse grid
 * \param *P pointer to the dCSRmat matrix of resulted interpolation, here just obtain the structure of resulted interpolation
 * \param *param pointer to AMG parameters
 *
 * Standard coarsening, refer to P475 in book
 * U. Trottenberg, C. W. Oosterlee, and A. Sch¨uller. Multigrid. Academic Press Inc., San Diego, CA, 2001. 
 * With contributions by A. Brandt, P. Oswald and K. St¨uben.
 */
void coarsening(dCSRmat *A, ivector *vertices, dCSRmat *P, AMG_param *param)
{	
	int i,j,k,l;
	int Cn=0;
	iCSRmat S; // strong n-couplings
	int row=A->row, col=A->col;
	int index;
	
	/** Step 1: generate S */
	generateS(A, &S, param);
	iCSRmat ST=getTransposeOfA(&S);
	
	/** Step 2: Standard coarsening algorithm */
	LinkList LoL_head;
	LinkList LoL_tail;
	int *lists, *where;
	LoL_head = NULL;
	LoL_tail = NULL;
	lists = (int*)calloc(A->row,sizeof(int));
	where = (int*)calloc(A->row,sizeof(int));
	
	/** first coarsening phase */
	/** 1. initialize lambda */
	int *lambda;
	lambda=(int *)calloc(A->row,sizeof(int));
	
	int maxlambda, maxnode;
	for(i=0;i<A->row;i++) lambda[i]=ST.IA[i+1]-ST.IA[i];
	
	vertices->row=A->row;
	vertices->val=(int*)calloc(vertices->row, sizeof(int));
	
	int num_left=0;
	/** 2. Before the following algorithm starts, filter out the variables which */
	/** have no connection at all and become special F-variables. */
	for(i=0;i<A->row;i++)
	{
		if((A->IA[i+1]-A->IA[i])==1)
		{
			vertices->val[i]=2; // set i as  special fine node
			lambda[i]=0;
		}
		else
		{
			vertices->val[i]=-1; // set i as  undecide node
			num_left++;
		}
	}
	
	/** 3. Set the variables with nonpositvie measure as F-variables */
	double measure, new_meas;
	for(i=0;i<A->row;i++)
	{
		measure=lambda[i];
		if(vertices->val[i]!=2)
		{
			if(measure>0)
			{
				enter_on_lists(&LoL_head, &LoL_tail, lambda[i], i, lists, where);
			}
			else
			{
				if(measure<0) printf("negative lambda!\n");
				vertices->val[i]=0; // set i as fine node
				for(k=S.IA[i];k<S.IA[i+1];k++)
				{
					j=S.JA[k];
					if(vertices->val[j]!=2)
					{
						if(j<i)
						{
							new_meas=lambda[j];
							if(new_meas>0)
							{
								remove_point(&LoL_head, &LoL_tail, new_meas, j, lists, where);
							}
							new_meas= ++(lambda[j]);
							enter_on_lists(&LoL_head, &LoL_tail,  new_meas, j, lists, where);
						}
						else
						{
							new_meas= ++(lambda[j]);
						}
					}
				}
				num_left--;
			}
		}
	}
	
	/** 4. Main loop */
	while(num_left>0)
	{
		/** pick $i\in U$ with $\max\lambda_i: C:=C\cup\{i\}, U:=U\\{i\}$ */
		maxnode=LoL_head->head;
		maxlambda=lambda[maxnode];
		
		vertices->val[maxnode]=1; /*! set maxnode as coarse node */
		lambda[maxnode]=0;
		num_left--;
		remove_point(&LoL_head, &LoL_tail, maxlambda, maxnode, lists, where);
		Cn++;
		
		/** for all $j\in S_i^T\cap U: F:=F\cup\{j\}, U:=U\backslash\{j\}$ */
		for(i=ST.IA[maxnode];i<ST.IA[maxnode+1];i++)
		{
			j=ST.JA[i];
			
			/** if j is unkown */
			if(vertices->val[j]==-1) 
			{
				vertices->val[j]=0;  /*! set j as fine node */
				remove_point(&LoL_head, &LoL_tail, lambda[j], j, lists, where);
				num_left--;
				
				for(l=S.IA[j];l<S.IA[j+1];l++)
				{
					k=S.JA[l];
					if(vertices->val[k]==-1) // k is unkown
					{
						remove_point(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
						new_meas= ++(lambda[k]);
						enter_on_lists(&LoL_head, &LoL_tail,new_meas, k, lists, where);
					}
				}
			} // if
		} // i
		
		for(i=S.IA[maxnode];i<S.IA[maxnode+1];i++)
		{
			j=S.JA[i];
			if(vertices->val[j]==-1) // j is unkown
			{
				measure=lambda[j];
				remove_point(&LoL_head, &LoL_tail,measure, j, lists, where);
				lambda[j]=--measure;
				if(measure>0)
				{
					enter_on_lists(&LoL_head, &LoL_tail,measure, j, lists, where);
				}
				else
				{
					vertices->val[j]=0; // set j as fine variable
					num_left--;
					
					for(l=S.IA[j];l<S.IA[j+1];l++)
					{
						k=S.JA[l];
						if(vertices->val[k]==-1) // k is unkown
						{
							remove_point(&LoL_head, &LoL_tail, lambda[k], k, lists, where);
							new_meas= ++(lambda[k]);
							enter_on_lists(&LoL_head, &LoL_tail,new_meas, k, lists, where);
						}
					}
				}
			}
		}
	} // while
	
	free(lists);
	free(where);
	LinkList   list_ptr;
	list_ptr =  LoL_head;
	if(LoL_head!=NULL)
	{
		list_ptr=LoL_head;
		LoL_head->prev_elt=NULL;
		LoL_head->next_elt=NULL;
		LoL_head = list_ptr -> next_elt;
		free(list_ptr);
	}
	
	/** second pass, check fine points for coarse neighbors */
	int		    ji, jj;
	int		    ci_tilde = -1;
	int		    ci_tilde_mark = -1;
	int		    set_empty = 1;
	int		    C_i_nonempty = 0;
	int       *graph_array;
	
	graph_array=(int*)calloc(row, sizeof(int));
	for (i = 0; i < row; i++) graph_array[i] = -1;
	
	for (i=0; i < row; i++)
	{
		if (ci_tilde_mark |= i) ci_tilde = -1;
		if (vertices->val[i] == 0) // F-variable
		{
			for (ji = S.IA[i]; ji < S.IA[i+1]; ji++)
			{
				j = S.JA[ji];
				if (vertices->val[j] ==1) // C-variable
					graph_array[j] = i;
			}
			for (ji = S.IA[i]; ji < S.IA[i+1]; ji++)
			{
				j = S.JA[ji];
				if (vertices->val[j] == 0)
				{
					set_empty = 1;
					for (jj = S.IA[j]; jj < S.IA[j+1]; jj++)
					{
						index = S.JA[jj];
						if (graph_array[index] == i)
						{
							set_empty = 0;
							break;
						}
					}
					if (set_empty)
					{
						if (C_i_nonempty)
						{
							vertices->val[i] = 1;
							Cn++;
							if (ci_tilde > -1)
							{
								vertices->val[ci_tilde] = 0;
								Cn--;
								ci_tilde = -1;
							}
							C_i_nonempty = 0;
							break;
						}
						else
						{
							ci_tilde = j;
							ci_tilde_mark = i;
							vertices->val[j] = 1;
							Cn++;
							C_i_nonempty = 1;
							i--;
							break;
						}
					}
				}
			}
		}
	}
	
	free(graph_array);
		
	/** Step 3: generate P which is the inersection of C and S */
	P->val=NULL;
	P->JA=NULL;
	P->row=A->row;
	P->col=Cn;
	P->IA=(int*)calloc(P->row+1, sizeof(int));	
	
	// step 1: Find first the structure IA of P		
	for(i=0;i<A->row;i++)
	{
		if(vertices->val[i]==0) // if node i is on fine grid
		{
			for(j=S.IA[i];j<S.IA[i+1];j++)
			{
				k=S.JA[j];
				if(vertices->val[k]==1)
				{
					P->IA[i+1]++;
				}
			}
		}
		else if(vertices->val[i]==2) // if node i is a special fine node
		{
			P->IA[i+1]=0;
		}
		else // if node i is on coarse grid 
		{
			P->IA[i+1]=1;
		}
	}	
	
	for(i=0;i<P->row;i++)
		P->IA[i+1]+=P->IA[i];
	
	P->nnz=P->IA[P->row];
	// step 2: Find the structure JA of P
	index=0;
	P->JA=(int*)calloc(P->nnz,sizeof(int));
	P->val=(double*)calloc(P->nnz,sizeof(double));
	for(i=0;i<A->row;i++)
	{
		if(vertices->val[i]==0) // fine node
		{
			for(j=S.IA[i];j<S.IA[i+1];j++)
			{
				k=S.JA[j];
				if(vertices->val[k]==1)
				{
					P->JA[index]=k;
					index++;
				}
			}
		}
		else if(vertices->val[i]==2) // if node i is a special fine node
		{
		}
		else // if node i is on coarse grid 
		{
			P->JA[index]=i;
			index++;
		}
	}	

	free(lambda);	
	free(S.IA);
	free(S.JA);
	free(ST.IA);
	free(ST.JA);
}

/**
 * \fn void generateS(dCSRmat *A, iCSRmat *S, AMG_param *param)
 * \brief generate the set of all strong couplings S
 * \param *A pointer to the coefficient matrix
 * \param *S pointer to the set of all strong couplings matrix
 * \param *param pointer to AMG parameters
 */
void generateS(dCSRmat *A, iCSRmat *S, AMG_param *param)
{
	double max_row_sum=param->max_row_sum;
	double epsilon_str=param->strong_threshold;
	int coarsening_type=param->coarsening_type;
	
	if(coarsening_type==2 || coarsening_type==3)
	{
		generateSRS(A, S, epsilon_str, coarsening_type);
		return;
	}
	
	int i,j,k,l;
	int row=A->row, col=A->col;
	int index;
	
	// get the diagnal entry of A
	dvector diag;
	getdiag(0, A, &diag);
	
	/** Step 1: generate S */
	double row_scale, row_sum;
	
	// copy the structure of A to S
	S->val=NULL;
	S->row=row;
	S->col=col;
	S->nnz=A->IA[row];
	S->IA=(int*)calloc(S->row+1, sizeof(int));
	S->JA=(int*)calloc(S->nnz, sizeof(int));
	for(i=0;i<=row;i++) S->IA[i]=A->IA[i];
	for(i=0;i<S->nnz;i++) S->JA[i]=A->JA[i];
	
	for(i=0;i<row;i++) {
		
		/** compute scaling factor and row sum */
		row_scale=0;
		row_sum=0;
		for(j=A->IA[i];j<A->IA[i+1];j++) {
			row_scale=min(row_scale, A->val[j]);
			row_sum+=A->val[j];
		}
		row_sum=fabs(row_sum/diag.val[i]);
		
		/* compute row entries of S */
		for(j=A->IA[i];j<A->IA[i+1];j++) {
			if(A->JA[j]==i) {
				S->JA[j]=-1;
				break;
			}
		}
		
		/** if $|\sum_{j=1}^n a_{ij}|> \theta_2 |a_{ii}|$ */
		if((row_sum>max_row_sum)&&(max_row_sum<1)) { 
			/** make all dependencies weak */
			for(j=A->IA[i];j<A->IA[i+1];j++) S->JA[j]=-1;
		}
		/** otherwise */
		else {
			for(j=A->IA[i];j<A->IA[i+1];j++) {
				if(A->val[j]>=epsilon_str*row_scale) S->JA[j]= -1; /*< if $a_{ij}>=\epsilon_{str}*\min a_{ij}$, the connection $a_{ij}$ is set to be weak connection */
			}
		}
	}
	
	free_dvector(&diag);
	
	/* Compress the strength matrix */
	index=0;
	for(i=0;i<row;i++) {
		S->IA[i]=index;
		for(j=A->IA[i];j<A->IA[i+1];j++) {
			if(S->JA[j]>-1) {
				S->JA[index]=S->JA[j];
				index++;
			}
		}
	}
	
	S->IA[row]=index;
	S->nnz=index;
	S->JA=(int*)realloc(S->JA, (S->IA[row])*sizeof(int));
}

/**
 * \fn void generateSRS(dCSRmat *A, iCSRmat *S, double epsilon_str, int coarsening_type)
 * \brief generate the set of all strong negative couplings(coarsening_type=2) or strong absolute couplings(coarsening_type=3) S
 * \param *A pointer to the coefficient matrix
 * \param *S pointer to the set of all strong couplings matrix
 * \param epsilon_str strong coupled ratio
 * \param coarsening_type coarsening type(2: strong negative couplings, 3: strong absolute couplings)
 */
void generateSRS(dCSRmat *A, iCSRmat *S, double epsilon_str, int coarsening_type)
{
	int i,j;
	double amax[A->row];
	
	// get the maximum absolute negative value data in each row of A
	if(coarsening_type==2) // original RS coarsening, just consider negative strong coupled
	{
		for(i=0;i<A->row;i++)
		{
			amax[i]=0;
			for(j=A->IA[i];j<A->IA[i+1];j++)
			{
				if((-A->val[j])>amax[i] && A->JA[j]!=i)
				{
					amax[i]=-A->val[j];
				}
			}
		}
	}
	
	if(coarsening_type==3) // original RS coarsening, consider absolute strong coupled
	{
		for(i=0;i<A->row;i++)
		{
			amax[i]=0;
			for(j=A->IA[i];j<A->IA[i+1];j++)
			{
				if(fabs(A->val[j])>amax[i] && A->JA[j]!=i)
				{
					amax[i]=fabs(A->val[j]);
				}
			}
		}
	}
	
	// step 1: Find first the structure IA of S
	S->val=NULL;
	S->JA=NULL;
	S->row=A->row;
	S->col=A->col;
	S->IA=(int*)calloc(S->row+1, sizeof(int));
	
	if(coarsening_type==2)
	{
		for(i=0;i<S->row;i++)
		{
			for(j=A->IA[i];j<A->IA[i+1];j++)
			{
				if((-A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i)
				{
					S->IA[i+1]++;
				}
			}	
		}
	}
	
	if(coarsening_type==3)
	{
		for(i=0;i<S->row;i++)
		{
			for(j=A->IA[i];j<A->IA[i+1];j++)
			{
				if(fabs(A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i)
				{
					S->IA[i+1]++;
				}
			}
		}
	}
	
	for(i=0;i<S->row;i++) S->IA[i+1]+=S->IA[i];
	
	// step 2: Find the structure JA of S
	int index=0;
	S->JA=(int*)calloc(S->IA[S->row],sizeof(int));
	if(coarsening_type==2)
	{
		for(i=0;i<S->row;i++)
		{
			for(j=A->IA[i];j<A->IA[i+1];j++)
			{
				if((-A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i)
				{
					S->JA[index]= A->JA[j];
					index++;
				}
			}
		}
	}
	
	if(coarsening_type==3)
	{
		for(i=0;i<S->row;i++)
		{
			for(j=A->IA[i];j<A->IA[i+1];j++)
			{
				if(fabs(A->val[j])>=epsilon_str*amax[i] && A->JA[j]!=i)
				{
					S->JA[index]= A->JA[j];
					index++;
				}
			}
		}
	}
}


// Ac=R*A*P: form coarse level stiffness matrix
void getCoarseStiffMatrix(dCSRmat *R, dCSRmat *A, dCSRmat *P, dCSRmat *Ac)
{
	sparseTripleMultiplication(R, A, P, Ac);
}
