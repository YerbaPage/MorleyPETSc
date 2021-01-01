/*
 *  refine.c
 *
 *  Created by Xuehai Huang on 08/27/2010.
 *  Copyright 2010 WZU. All rights reserved.
 *
 */

/*! \file refine.c
 *  \brief mark and refine
 */

#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "matvec.h"


/**
 * \fn double mark(idenmat *elementEdge, EDGE *edges, dennode *nodes, dvector *estimator, int *index, double theta, double totalError, int *marker)
 * \brief Dorfler marking strategy
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *estimator pointer to error estimator in ascending order
 * \param *index pointer to index of estimator
 * \param theta parameter in (0,1), usually chosen in [0.2,0.5]. 
 *           We mark minimal number of triangles M such that
 *             \sum_{T \in M} \eta_T^2 > \theta*\sum\eta_T^2
 * \param totalError total posterior error
 * \param *marker pointer to marking status of edges (posivite: marked, 0: unmarked)
 * \return total posterior error
 */
double mark(idenmat *elementEdge, EDGE *edges, dennode *nodes, dvector *estimator, int *index, double theta, double totalError, int *marker)
{
	int i, j, k;
	int Nt=estimator->row;
	int Ne=edges->row;
	int nvertices=nodes->row;
	int base, element;
	double current=0;
	int *isAdded;

	isAdded=(int*)calloc(Nt, sizeof(int));
	for(k=Nt-1;k>-1;k--)
	{
		if(current>theta*totalError)
			break;

		element=index[k];
		while(1)
		{
			base=elementEdge->val[element][0];
			if(marker[base]>0)
			{
				if(isAdded[element]==0)
				{
					current+=estimator->val[element];
					isAdded[element]=1;
				}
				break;
			}
			else
			{
				current+=estimator->val[element];
				isAdded[element]=1;
				marker[base]=nodes->row;
				nodes->row++;
				element=edges->val[base][2]+edges->val[base][3]-element;
				if(element==-1)
					break;
			}
		}
	}
	free(isAdded);

	// generate nodes information
	nodes->bdFlag = (int*)realloc(nodes->bdFlag, sizeof(int)*(nodes->row));
	nodes->val=(double**)realloc(nodes->val, sizeof(double *)*(nodes->row));
	nodes->tensorBasis = (double***)realloc(nodes->tensorBasis, sizeof(double **)*(nodes->row));
	for (i = nvertices; i < nodes->row; i++)
	{
		nodes->val[i] = (double*)calloc(nodes->col, sizeof(double));
		nodes->tensorBasis[i] = (double**)calloc(3, sizeof(double *));
		for (j = 0; j<3; j++)
			nodes->tensorBasis[i][j] = (double*)calloc(3, sizeof(double));
	}
	int point1, point2;
	for(base=0;base<Ne;base++)
	{
		point1=edges->val[base][0];
		point2=edges->val[base][1];
		k=marker[base];
		if(k>0)
		{
			nodes->val[k][0]=(nodes->val[point1][0]+nodes->val[point2][0])/2;
		    nodes->val[k][1]=(nodes->val[point1][1]+nodes->val[point2][1])/2;
			nodes->bdFlag[k] = edges->bdFlag[base];
		}
	}

	return current;
}


/**
 * \fn void refine(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, iCSRmat *edgesTran, dennode *nodes, int *marker)
 * \brief  refine the marked mesh using the newest vertex bisection
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                       the fourth column stores -1 if the edge is on boundary
 * \param *edgesTran the relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *marker pointer to marking status of edges (posivite: marked, 0: unmarked)
 * \return void
 */
void refine(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, iCSRmat *edgesTran, dennode *nodes, int *marker)
{
	int i,k;
	int Nt=elements->row;
//	int Ne=edges->row;
	int base, left, right, element;

	for(k=0;k<Nt;k++)
	{
		for(i=0;i<3;i++)
		{
			if(marker[elementEdge->val[k][i]]>0)
				elements->row++;
		}
	}
	
	elements->val=(int**)realloc(elements->val, sizeof(int *)*(elements->row));
	elements->parent = (int*)realloc(elements->parent, sizeof(int)*(elements->row));
	elements->xi=(double**)realloc(elements->xi, sizeof(double *)*(elements->row));
	elements->eta=(double**)realloc(elements->eta, sizeof(double *)*(elements->row));
	elements->edgeslength = (double**)realloc(elements->edgeslength, sizeof(double *)*(elements->row));
	elements->nvector = (double***)realloc(elements->nvector, sizeof(double **)*(elements->row));
	elements->tvector = (double***)realloc(elements->tvector, sizeof(double **)*(elements->row));
	for(k=Nt;k<elements->row;k++)
	{
		elements->val[k]=(int*)calloc(elements->col, sizeof(int));
		elements->xi[k]=(double*)calloc(elements->col, sizeof(double));
		elements->eta[k]=(double*)calloc(elements->col, sizeof(double));
		elements->edgeslength[k] = (double*)calloc(elements->col, sizeof(double));
		elements->nvector[k] = (double**)calloc(elements->col, sizeof(double *));
		elements->tvector[k] = (double**)calloc(elements->col, sizeof(double *));
		for (i = 0; i<elements->col; i++)
		{
			elements->nvector[k][i] = (double*)calloc(elements->col - 1, sizeof(double));
			elements->tvector[k][i] = (double*)calloc(elements->col - 1, sizeof(double));
		}

	}
	elements->vol=(double*)realloc(elements->vol, sizeof(double)*(elements->row));
	
	// newest vertex bisection
	int count=Nt;
	int vertices[3];
	for(k=0;k<Nt;k++)
	{
		for(i=0;i<3;i++)
			vertices[i]=elements->val[k][i];
		base=elementEdge->val[k][0];
		if(marker[base]>0) // base
		{
			right=elementEdge->val[k][1];
			left=elementEdge->val[k][2];
			if(marker[right]>0) // base,right
			{
				if(marker[left]>0) // base,right,left
				{
					elements->val[count][0]=marker[right];
					elements->val[count][1]=marker[base];
					elements->val[count][2]=vertices[2];
					elements->parent[count] = k;
					count++;
					elements->val[count][0]=marker[left];
					elements->val[count][1]=vertices[1];
					elements->val[count][2]=marker[base];
					elements->parent[count] = k;
					count++;
					elements->val[count][0]=marker[right];
					elements->val[count][1]=vertices[0];
					elements->val[count][2]=marker[base];
					elements->parent[count] = k;
					count++;
					elements->val[k][0]=marker[left];
					elements->val[k][1]=marker[base];
					elements->val[k][2]=vertices[0];
					elements->parent[k] = k;
				}
				else // base,right
				{
					elements->val[count][0]=marker[right];
					elements->val[count][1]=marker[base];
					elements->val[count][2]=vertices[2];
					elements->parent[count] = k;
					count++;
					elements->val[count][0]=marker[right];
					elements->val[count][1]=vertices[0];
					elements->val[count][2]=marker[base];
					elements->parent[count] = k;
					count++;
					elements->val[k][0]=marker[base];
					elements->val[k][1]=vertices[0];
					elements->val[k][2]=vertices[1];
					elements->parent[k] = k;
				}
			}
			else
			{
				if(marker[left]>0) // base,left
				{
					elements->val[count][0]=marker[base];
					elements->val[count][1]=vertices[2];
					elements->val[count][2]=vertices[0];
					elements->parent[count] = k;
					count++;
					elements->val[count][0]=marker[left];
					elements->val[count][1]=vertices[1];
					elements->val[count][2]=marker[base];
					elements->parent[count] = k;
					count++;
					elements->val[k][0]=marker[left];
					elements->val[k][1]=marker[base];
					elements->val[k][2]=vertices[0];
					elements->parent[k] = k;
				}
				else // base
				{
					elements->val[count][0]=marker[base];
					elements->val[count][1]=vertices[2];
					elements->val[count][2]=vertices[0];
					elements->parent[count] = k;
					count++;
					elements->val[k][0]=marker[base];
					elements->val[k][1]=vertices[0];
					elements->val[k][2]=vertices[1];
					elements->parent[k] = k;
				} // if
			} // if
		} // if
	} // for
	
	// get edge information
	iCSRmat elementsTran;
	getTransposeOfELEMENT(elements, &elementsTran, elements->col, nodes->row);
	getEdgeInfo(elements, &elementsTran, edges, edgesTran, nodes);
	free(elementsTran.IA);
	free(elementsTran.JA);

	// get elementEdge information
	int point1, point2, element1, element2;
	int l;
	free_iden_matrix(elementEdge);
	create_iden_matrix(elements->row, 3, elementEdge);
	for(i=0;i<edges->row;i++)
	{
		point1=edges->val[i][0];
		point2=edges->val[i][1];
		element1=edges->val[i][2];
		element2=edges->val[i][3];
		
		for(k=0;k<3;k++)
		{
			if(elements->val[element1][k]==point1)
				break;
		}
		for(l=0;l<3;l++)
		{
			if(elements->val[element1][l]==point2)
				break;
		}
		elementEdge->val[element1][3-k-l]=i;
		
		if(element2>-1)
		{
			for(k=0;k<3;k++)
			{
				if(elements->val[element2][k]==point1)
					break;
			}
			for(l=0;l<3;l++)
			{
				if(elements->val[element2][l]==point2)
					break;
			}
			elementEdge->val[element2][3-k-l]=i;
		}
	}
}