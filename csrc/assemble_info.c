/*
 *  assemble_info.c
 *
 *  Created by Xuehai Huang on 3/29/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file assemble_info.c
 *  \brief Subroutines for assembling purpose
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn void getEdgeInfo(ELEMENT *elements, iCSRmat *elementsTran, EDGE *edges, iCSRmat *edgesTran)
 * \brief get edges information with edgesTran from elements
 *			 ALgorithm: node-->element-->edge
 * \param *elements store 3 nodes corresponding to the element, the 4th column store the relation with coarse grid elements( store -1 if itself is the coarset grid)
 * \param *elementsTran pointer to the transpose of *elements
 * \param *edges the first two columns store the two vertice corresponding to the edge; 
 *				 the 3rd and 4th columns store the triangles which the edge belongs to;
 *				 if the edge is a boundary, the 4th column will stores -1;
 *				 the first column is in ascend order.
 * \param *edgesTran the relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodes pointer to the nodes location of the triangulation
 * \return void
 */
void getEdgeInfo(ELEMENT *elements, iCSRmat *elementsTran, EDGE *edges, iCSRmat *edgesTran, dennode *nodes)
{
	int i,j,k,l;
	int element, point1, point2, col;	
	
	int iSTART=0;
	int curEdgeRow = 0;

	// Euler's formula: F - E + V = 1
	create_EDGE(elements->row + elementsTran->row - 1, 4, edges);
	edgesTran->row=elementsTran->row;
	edgesTran->IA=(int*)calloc(edgesTran->row+1, sizeof(int));
	col=elements->col;
	for(i=0;i<elementsTran->row;i++)
	{
		for(j=elementsTran->IA[i];j<elementsTran->IA[i+1];j++)
		{
			element=elementsTran->JA[j]; // achieve the element

			for(k=0;k<col;k++)
			{
				if(elements->val[element][k]==i)
				{
					point1=elements->val[element][(k-1+col)%col];
					point2=elements->val[element][(k+1)%col];
					break;
				}
			} // k
			
			if(i<point1) // the case i>point1 already exists in edges
			{
				for(l=iSTART;l<curEdgeRow;l++)
				{					
					if((edges->val[l][0]==i) && (edges->val[l][1]==point1))
					{
						break;
					}
				}
				if(l == curEdgeRow) // the edge(i, point1) is new 
				{
					edges->val[l][0]=i;
					edges->val[l][1]=point1;
					edges->val[l][2]=element;
					edges->val[l][3]=-1;
					edgesTran->IA[i+1]++;
					edgesTran->IA[point1+1]++;
					curEdgeRow++;
				}
				else // the edge(i, point1) already exits
				{
					edges->val[l][3]=element;
				}				
			}
			
			if(i<point2) // the case i>point2 already exists in edges
			{
				for(l=iSTART;l<curEdgeRow;l++)
				{					
					if((edges->val[l][0]==i) && (edges->val[l][1]==point2))
					{
						break;
					}
				}
				if(l == curEdgeRow) // the edge(i, point2) is new
				{
					edges->val[l][0]=i;
					edges->val[l][1]=point2;
					edges->val[l][2]=element;
					edges->val[l][3]=-1;
					edgesTran->IA[i+1]++;
					edgesTran->IA[point2+1]++;
					curEdgeRow++;
				}
				else // the edge(i, point2) already exits
				{
					edges->val[l][3]=element;
				}				
			}
		} // j
		
		iSTART = curEdgeRow;
	} // i
	
	// generate egdesTran
	edgesTran->col=edges->row;
	for(i=0;i<edgesTran->row;i++)
	{
		edgesTran->IA[i+1]=edgesTran->IA[i+1]+edgesTran->IA[i];
	}
	edgesTran->nnz = edgesTran->IA[edgesTran->row];
	edgesTran->JA=(int*)calloc(edgesTran->nnz, sizeof(int));
	edgesTran->val=(int*)calloc(edgesTran->nnz, sizeof(int));
	
	int *count; // for counting in each row
	count=(int*)calloc(edgesTran->row, sizeof(int));
	for(i=0;i<edges->row;i++)
	{
		point1=edges->val[i][0];
		point2=edges->val[i][1];
		edgesTran->JA[edgesTran->IA[point1]+count[point1]]=i;
		edgesTran->val[edgesTran->IA[point1]+count[point1]]=point2;
		count[point1]++;
		edgesTran->JA[edgesTran->IA[point2]+count[point2]]=i;
		edgesTran->val[edgesTran->IA[point2]+count[point2]]=point1;
		count[point2]++;
	}
	free(count);

	getEdgeBdflagFromVertexBdflag(edges, nodes);
}

/**
* \fn void getEdgeBdflagFromVertexBdflag(EDGE *edges, dennode *nodes)
* \brief get edges boundary information from vertices boundary information
* \param *edges the first two columns store the two vertice corresponding to the edge;
*				 the 3rd and 4th columns store the triangles which the edge belongs to;
*				 if the edge is a boundary, the 4th column will stores -1;
*				 the first column is in ascend order.
* \param *nodes pointer to the nodes location of the triangulation
* \return void
************* boundary type of vertex *************
0 : non-boundary, i.e., an interior vertex.
1 : first type, i.e., a Dirichlet boundary vertex.
2 : second type, i.e., a Neumann boundary vertex.
3 : third type, i.e., a Robin boundary vertex.
12: a Dirichlet-Neumann boundary vertex.
22: a Neumann-Neumann boundary vertex.
************* boundary type of edge *************
* 0 : non - boundary, i.e., an interior edge or face.
* 1 : first type, i.e., a Dirichlet boundary edge or face.
* 2 : second type, i.e., a Neumann boundary edge or face.
* 3 : third type, i.e., a Robin boundary edge or face.
*/
void getEdgeBdflagFromVertexBdflag(EDGE *edges, dennode *nodes)
{
	int *vbdflag = nodes->bdFlag;
	int *ebdflag = edges->bdFlag;
	int i, j, k, temp;
	int ev[2];
	int lflag[2], rflag[2];

	for (i = 0; i < edges->row; i++)
	{
		if (edges->val[i][3] == -1)
		{
			ev[0] = edges->val[i][0];
			ev[1] = edges->val[i][1];
			lflag[0] = vbdflag[ev[0]] / 10;
			lflag[1] = vbdflag[ev[0]] % 10;
			rflag[0] = vbdflag[ev[1]] / 10;
			rflag[1] = vbdflag[ev[1]] % 10;
			if (lflag[0] > lflag[1])
			{
				temp = lflag[0];
				lflag[0] = lflag[1];
				lflag[1] = temp;
			}
			if (rflag[0] > rflag[1])
			{
				temp = rflag[0];
				rflag[0] = rflag[1];
				rflag[1] = temp;
			}
			
			if (lflag[0] > 0)
			{
				if (rflag[0] > 0)
				{
					if (lflag[0] == rflag[0] || lflag[0] == rflag[1])
						ebdflag[i] = lflag[0];
					else if (lflag[1] == rflag[0] || lflag[1] == rflag[1])
						ebdflag[i] = lflag[1];
					else
						ebdflag[i] = 0; // impossible
				}
				else
				{
					if (rflag[1] == lflag[0] || rflag[1] == lflag[1])
						ebdflag[i] = rflag[1];
					else
						ebdflag[i] = 0; // impossible
				}
			}
			else
			{
				if (lflag[1] == rflag[0] || lflag[1] == rflag[1])
					ebdflag[i] = lflag[1];
				else
					ebdflag[i] = 0; // impossible
			}
		/*	if (vbdflag[ev[0]] == 1 && vbdflag[ev[1]] == 1) // Dirichlet & Dirichlet
				ebdflag[i] = 1;

			else if (vbdflag[ev[0]] == 2 && vbdflag[ev[1]] == 2) // Neumann & Neumann
				ebdflag[i] = 2;
			
			else if (vbdflag[ev[0]] == 22 && vbdflag[ev[1]] == 22) // Neumann-Neumann & Neumann-Neumann
				ebdflag[i] = 2;
			
			else if (vbdflag[ev[0]] == 1 && vbdflag[ev[1]] == 12 || vbdflag[ev[0]] == 12 && vbdflag[ev[1]] == 1) // Dirichlet & Dirichlet-Neumann
				ebdflag[i] = 1;
			
			else if (vbdflag[ev[0]] == 2 && vbdflag[ev[1]] == 12 || vbdflag[ev[0]] == 12 && vbdflag[ev[1]] == 2) // Neumann & Dirichlet-Neumann
				ebdflag[i] = 2;
			
			else if (vbdflag[ev[0]] == 12 && vbdflag[ev[1]] == 22 || vbdflag[ev[0]] == 22 && vbdflag[ev[1]] == 12) // Dirichlet-Neumann & Neumann-Neumann
				ebdflag[i] = 2;
			
			else if (vbdflag[ev[0]] == 2 && vbdflag[ev[1]] == 22 || vbdflag[ev[0]] == 22 && vbdflag[ev[1]] == 2) // Neumann & Neumann-Neumann
				ebdflag[i] = 2;
			*/
		}
		else
			ebdflag[i] = 0;
	}

}

/**
 * \fn int getCoarseInfo(int domain_num, dennode *nodes, ELEMENT *elements, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge)
 * \brief generate the coarse grid information and store it into point, element, edge respectively
 * \param domain_num number of domain
 * \param *nodes the first column stores the x coordinate of points, the second column stores the y coordinate of points
 * \param *elements store 3 nodes corresponding to the element, the 4th column store the relation with coarse grid elements( store -1 if itself is the coarset grid)
 * \param *edges the first two columns store the two vertice, the third column stores the affiliated element
 * \param *edgesTran the relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \return 1 if succeed 0 if fail
 */
int getCoarseInfo(int domain_num, dennode *nodes, ELEMENT *elements, EDGE *edges, iCSRmat *edgesTran, ivector *nodeCEdge)
{
	// get data from inputFile
	char *str1 = "data/unitsquare.dat";
	char *str2 = "data/square.dat";
	char *str3 = "data/adapmesh.dat";
	char *str4 = "data/unstructuredmesh1.dat";
	char *str5 = "data/LDomain.dat";
	char *str6 = "data/rotLDomain.dat";
	char *str7 = "data/tempLDomain.dat";

	char *filename;
	switch (domain_num) {
	case 1:filename = str1; break;
	case 2:filename = str2; break;
	case 3:filename = str3; break;
	case 4:filename = str4; break;
	case 5:filename = str5; break;
	case 6:filename = str6; break;
	case 7:filename = str7; break;
	default:filename = str1;
	}

	FILE *inputFile;
	inputFile=fopen(filename, "r");
	if(inputFile==NULL)
	{
		printf("Opening file %s fails!\n", filename);
		return 0;
	}
	
	int NumNode, Ncoor, NumElem, Nnode;
	int i,j;
	
	// get the nodes' coordinates
	fscanf(inputFile, "%d %d", &NumNode, &Ncoor);
	create_dennode(NumNode, Ncoor, nodes);
	for(i=0;i<nodes->row;i++)
	{
		for(j=0;j<nodes->col;j++)
			fscanf(inputFile, "%lf", &nodes->val[i][j]);
		fscanf(inputFile, "%d", &nodes->bdFlag[i]);
	}
	
	// get triangular grid
	idenmat T;
	fscanf(inputFile, "%d %d", &NumElem, &Nnode);
	create_ELEMENT(NumElem, Nnode, elements);
		
	for(i=0;i<elements->row;i++)
	{
		for(j=0;j<elements->col;j++)
		{
			fscanf(inputFile, "%d", &elements->val[i][j]);
			elements->val[i][j]--; // the data from matlab start with 1 while from c start with 0
		}
		elements->parent[i]=-1;
	}	
	
	fclose(inputFile);
	
	iCSRmat elementsTran; // the transpose of *elements

	getTransposeOfELEMENT(elements, &elementsTran, elements->col, nodes->row);
	
	// get edge information
	getEdgeInfo(elements, &elementsTran, edges, edgesTran, nodes);

	free(elementsTran.IA);
	free(elementsTran.JA);
		
	// get nodeCEdge
	nodeCEdge->row=nodes->row;
	nodeCEdge->val=(int*)calloc(nodeCEdge->row, sizeof(int));
	for(i=0;i<nodeCEdge->row;i++)
		nodeCEdge->val[i]=-1;
	
	return 1;
}

/**
 * \fn void uniformrefine(dennode *Cnodes, ELEMENT *Celements, EDGE *Cedges, dennode *Fnodes, ELEMENT *Felements, EDGE *Fedges, iCSRmat *FedgesTran, ivector *nodeCEdge)
 * \brief generate fine grid using regular section
 * \param *Cnodes pointer to the nodes location of the triangulation on coasre grid
 * \param *Celements pointer to the structure of the triangulation on coasre grid
 * \param *Cedges pointer to the edge information of the triangulation on coasre grid
 * \param *Fnodes pointer to the nodes location of the triangulation on fine grid
 * \param *Felements pointer to the structure of the triangulation on fine grid
 * \param *Fedges pointer to the edge information of the triangulation on fine grid
 * \param *FedgesTran the relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \return void
 */
void uniformrefine(dennode *Cnodes, ELEMENT *Celements, EDGE *Cedges, dennode *Fnodes, ELEMENT *Felements, EDGE *Fedges, iCSRmat *FedgesTran, ivector *nodeCEdge)
{
	iCSRmat FelementsTran; // the transpose of *Felements
	int i,j,k,l;
	int NumCNodes;
	NumCNodes=Cnodes->row;
	
	// generate fine gird's nodes information
	int midElements[Celements->row][3]; // store the 3 middle point of each tirangle

	create_dennode(Cnodes->row + Cedges->row, Cnodes->col, Fnodes);
	nodeCEdge->row=Fnodes->row;
	nodeCEdge->val=(int*)realloc(nodeCEdge->val, sizeof(int)*(nodeCEdge->row));

	// copy Cnodes into Fnodes
	for (i = 0; i < Cnodes->row; i++)
	{
		for (j = 0; j < Cnodes->col; j++)
			Fnodes->val[i][j] = Cnodes->val[i][j];
		Fnodes->bdFlag[i] = Cnodes->bdFlag[i];
	}
	
	int col=Celements->col;
	int point1, point2,location;
	for(i=0;i<Cedges->row;i++)
	{
		point1=Cedges->val[i][0];
		point2=Cedges->val[i][1];
		Fnodes->val[NumCNodes+i][0]=(Fnodes->val[point1][0]+Fnodes->val[point2][0])/2;
		Fnodes->val[NumCNodes+i][1]=(Fnodes->val[point1][1]+Fnodes->val[point2][1])/2;
		Fnodes->bdFlag[NumCNodes + i] = Cedges->bdFlag[i];
		nodeCEdge->val[NumCNodes+i]=i;
		
		// get the relation between triangle and edge
		l=Cedges->val[i][2];
		for(j=0;j<col;j++)
		{
			if(Celements->val[l][j]==point1)
			{
				break;
			}
		}
		for(k=0;k<col;k++)
		{
			if(Celements->val[l][k]==point2)
			{
				break;
			}
		}
		location=3-j-k;
		midElements[l][location]=i;
		if(Cedges->val[i][3]>-1) // case this edge is not a boundary edge
		{
			l=Cedges->val[i][3];
			for(j=0;j<col;j++)
			{
				if(Celements->val[l][j]==point1)
				{
					break;
				}
			}
			for(k=0;k<col;k++)
			{
				if(Celements->val[l][k]==point2)
				{
					break;
				}
			}
			location=3-j-k;
			midElements[l][location]=i;
		}
	}
	
	// generate fine grid trianglues information
	create_ELEMENT(Celements->row * 4, Celements->col, Felements);
	
	for(i=0;i<Celements->row;i++)
	{
		// bisection
/*		 Felements->val[i*4][0] = midElements[i][2]+NumCNodes;
		 Felements->val[i*4][1] = midElements[i][0]+NumCNodes;
		 Felements->val[i*4][2] = Celements->val[i][0];
		 Felements->parent[i*4] = i;
		 Felements->val[i*4+1][0] = midElements[i][2]+NumCNodes;
		 Felements->val[i*4+1][1] = Celements->val[i][1];
		 Felements->val[i*4+1][2] = midElements[i][0]+NumCNodes;
		 Felements->parent[i*4+1] = i;
		 Felements->val[i*4+2][0] = midElements[i][1]+NumCNodes;
		 Felements->val[i*4+2][1] = midElements[i][0]+NumCNodes;
		 Felements->val[i*4+2][2] = Celements->val[i][2];
		 Felements->parent[i*4+2] = i;
		 Felements->val[i*4+3][0] = midElements[i][1]+NumCNodes;
		 Felements->val[i*4+3][1] = Celements->val[i][0];
		 Felements->val[i*4+3][2] = midElements[i][0]+NumCNodes;
		 Felements->parent[i*4+3] = i; */
		
		// regular section start
		Felements->val[i*4][0] = Celements->val[i][0];
		Felements->val[i*4][1] = midElements[i][2]+NumCNodes;
		Felements->val[i*4][2] = midElements[i][1]+NumCNodes;
		Felements->parent[i*4] = i;
		Felements->val[i*4+1][0] = midElements[i][2]+NumCNodes;
		Felements->val[i*4+1][1] = Celements->val[i][1];
		Felements->val[i*4+1][2] = midElements[i][0]+NumCNodes;
		Felements->parent[i*4+1] = i;
		Felements->val[i*4+2][0] = midElements[i][1]+NumCNodes;
		Felements->val[i*4+2][1] = midElements[i][0]+NumCNodes;
		Felements->val[i*4+2][2] = Celements->val[i][2];
		Felements->parent[i*4+2] = i;
		Felements->val[i*4+3][0] = midElements[i][0]+NumCNodes;
		Felements->val[i*4+3][1] = midElements[i][1]+NumCNodes;
		Felements->val[i*4+3][2] = midElements[i][2]+NumCNodes;  
		Felements->parent[i*4+3] = i;
		// regular section end
	}
	
	getTransposeOfELEMENT(Felements, &FelementsTran, Felements->col, Fnodes->row);

	getEdgeInfo(Felements, &FelementsTran, Fedges, FedgesTran, Fnodes);

	free(FelementsTran.IA);
	free(FelementsTran.JA);
}

/**
* \fn void extractNondirichletMatrixVector(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index, dvector *uh)
* \brief extract stiffness matrix by removing the corresponding dirichlet boundary condition
* \param *A pointer to the stiffness matrix with dirichelt boundary condition(without removed)
* \param *b pointer to the right hand side with dirichelt boundary condition(without removed)
* \param *A11 pointer to the stiffness matrix without dirichelt boundary condition(removed)
* \param *b1 pointer to the right hand side without dirichelt boundary condition(removed)
* \param *isInNode if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *dirichlet pointer to the indicator of the dirichlet boundary
* \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
* \param *index pointer to the transpose of *dirichlet and *nondirichlet
* \param *uh pointer to the dirichlet boundary value
* \return void
*/
void extractNondirichletMatrixVector(dCSRmat *A, dvector *b, dCSRmat *A11, dvector *b1, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index, dvector *uh)
{
	// achiveve b1 due to dirichlet boundary condition
	extractNondirichletVector(A, b, b1, dirichlet, nondirichlet, uh);

	// achiveve A11 due to dirichlet boundary condition
	extractNondirichletMatrix11(A, A11, isInNode, dirichlet, nondirichlet, index);
}

/**
* \fn void extractNondirichletMatrix11(dCSRmat *A, dCSRmat *A11, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
* \brief extract stiffness matrix by removing the corresponding dirichlet boundary condition
* \param *A pointer to the stiffness matrix with dirichelt boundary condition(without removed)
* \param *A11 pointer to the stiffness matrix without dirichelt boundary condition(removed)
* \param *isInNode if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *dirichlet pointer to the indicator of the dirichlet boundary
* \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
* \param *index pointer to the transpose of *dirichlet and *nondirichlet
* \return void
*/
void extractNondirichletMatrix11(dCSRmat *A, dCSRmat *A11, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
{
	// achiveve A11 due to dirichlet boundary condition
	int i, j, k, l, i1, j1;
	int count;
	A11->row = A->row - dirichlet->row;
	A11->col = A11->row;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	// form A11->IA
	for (i = 0; i<A11->row; i++)
	{
		l = nondirichlet->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (isInNode->val[j] != -1)
			{
				A11->IA[i + 1]++;
			}
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++)
	{
		l = nondirichlet->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (isInNode->val[j] != -1)
			{
				A11->JA[count] = index->val[j];
				A11->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
 * \fn void extractNondirichletMatrix1(dCSRmat *A, dCSRmat *A1, ivector *dirichlet, ivector *nondirichlet)
 * \brief extract stiffness matrix by removing the corresponding dirichlet boundary condition  
 * \param *A pointer to the stiffness matrix with dirichelt boundary condition(without removed)
 * \param *A1 pointer to the stiffness matrix without dirichelt boundary condition(removed)
 * \param *dirichlet pointer to the indicator of the dirichlet boundary
 * \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
 * \return void
 */
void extractNondirichletMatrix1r(dCSRmat *A, dCSRmat *A1, ivector *dirichlet, ivector *nondirichlet)
{
	// achiveve A1 due to dirichlet boundary condition
	int i,k,l;
	int count;
	A1->row=A->row-dirichlet->row;
	A1->col=A->col;
	A1->IA=(int*)calloc(A1->row+1, sizeof(int));
	A1->JA=NULL;
	A1->val=NULL;
	
	// form A1->IA
	for(i=0;i<A1->row;i++)
	{
		l=nondirichlet->val[i];
		A1->IA[i+1]=A->IA[l+1]-A->IA[l];
	}
	
	for(i=0;i<A1->row;i++)
		A1->IA[i+1]+=A1->IA[i];

	A1->nnz=A1->IA[A1->row];
	
	// form A1->JA, A1->val
	A1->JA=(int*)calloc(A1->nnz, sizeof(int));
	A1->val=(double*)calloc(A1->nnz, sizeof(double));
	count=0;
	for(i=0;i<A1->row;i++)
	{
		l=nondirichlet->val[i];
		for(k=A->IA[l];k<A->IA[l+1];k++)
		{
			A1->JA[count]=A->JA[k];
			A1->val[count]=A->val[k];
			count++;
		}
	}	
}

/**
* \fn void extractNondirichletMatrix1c(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index)
* \brief extract stiffness matrix by removing the corresponding dirichlet boundary condition
* \param *A pointer to the stiffness matrix with dirichelt boundary condition(without removed)
* \param *A1 pointer to the stiffness matrix without dirichelt boundary condition(removed)
* \param *isInNode if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *dirichlet pointer to the indicator of the dirichlet boundary
* \param *index pointer to the transpose of *dirichlet and *nondirichlet
* \return void
*/
void extractNondirichletMatrix1c(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index)
{
	// achiveve A1 due to dirichlet boundary condition
	int i, j, k;
	int count;
	A1->row = A->row;
	A1->col = A->col - dirichlet->row;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (isInNode->val[j] != -1)
			{
				A1->IA[i + 1]++;
			}
		}
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (isInNode->val[j] != -1)
			{
				A1->JA[count] = index->val[j];
				A1->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void extractNondirichletMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index)
* \brief extract stiffness matrix by removing the corresponding dirichlet boundary condition
* \param *A pointer to the stiffness matrix with dirichelt boundary condition(without removed)
* \param *A1 pointer to the stiffness matrix without dirichelt boundary condition(removed)
* \param *isInNode if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
* \param *dirichlet pointer to the indicator of the dirichlet boundary
* \param *index pointer to the transpose of *dirichlet and *nondirichlet
* \return void
*/
void extractNondirichletMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ivector *isInNode, ivector *dirichlet, ivector *index)
{
	// achiveve A1 due to dirichlet boundary condition
	int i, j, k;
	int count;
	A1->row = A->row;
	A1->col = A->col - dirichlet->row*2;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	int col = A->col / 2;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (isInNode->val[j%col] != -1)
			{
				A1->IA[i + 1]++;
			}
		}
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (isInNode->val[j%col] != -1)
			{
				A1->JA[count] = index->val[j%col] + (index->row - dirichlet->row)*(j / col);
				A1->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
 * \fn void extractNondirichletVector(dCSRmat *A, dvector *b, dvector *b1, ivector *dirichlet, ivector *nondirichlet, dvector *uh)
 * \brief extract vector by removing the corresponding dirichlet boundary condition  
 * \param *b pointer to the vector with dirichelt boundary condition(without removed)
 * \param *b1 pointer to the vector without dirichelt boundary condition(removed)
 * \param *dirichlet pointer to the indicator of the dirichlet boundary
 * \param *nondirichlet pointer to the indicator of the node which is not in the dirichlet boundary
 * \return void
 */
void extractNondirichletVector(dCSRmat *A, dvector *b, dvector *b1, ivector *dirichlet, ivector *nondirichlet, dvector *uh)
{
	int i,j,k,i1,j1;
	b1->row=b->row-dirichlet->row;
	b1->val=(double*)calloc(b1->row, sizeof(double));	
		
	// achieve b1 due to dirichlet boundary condition
	for (i1 = 0; i1<b1->row; i1++)
	{
		i = nondirichlet->val[i1];
		b1->val[i1] = b->val[i];

		if (uh != NULL)
		{
			for (j1 = 0; j1 < dirichlet->row; j1++)
			{
				j = dirichlet->val[j1];
				for (k = A->IA[i]; k < A->IA[i + 1]; k++)
				{
					if (A->JA[k] == j)
					{
						b1->val[i1] -= A->val[k] * uh->val[j];
						break;
					}
				}
			}
		}
	}
}

/**
* \fn void extractFreenodesVector2StokesDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
* \brief extract vector according to freenodes
* \param *A pointer to left-top stiffness submatrix
* \param *B pointer to left-bottom stiffness submatrix
* \param *b pointer to the orginal vector
* \param *b1 pointer to the vector related to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \param *uh pointer to the initial solution
* \return void
*/
void extractFreenodesVector2StokesDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
{
	int i, j, k, i1, j1;

	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *nfFlag = &elementDOF->nfFlag;

	create_dvector(freenodes->row, &b1[0]);

	if (uh == NULL)
		return;

	// achieve b1[0] due to freenodes
	for (i1 = 0; i1<b1[0].row; i1++)
	{
		i = freenodes->val[i1];
		b1[0].val[i1] = b->val[i];

		for (j1 = A->IA[i]; j1 < A->IA[i + 1]; j1++)
		{
			j = A->JA[j1];
			if (nfFlag->val[j] == 1)
			{
				b1[0].val[i1] -= A->val[j1] * uh->val[j];
			}
		}
	}

	// achieve b1[1] due to freenodes
	for (i = 0; i<b1[1].row; i++)
	{
		for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
		{
			j = B->JA[j1];
			if (nfFlag->val[j] == 1)
			{
				b1[1].val[i] -= B->val[j1] * uh->val[j];
			}
		}
	}
}


/**
* \fn void extractFreenodesVector2ElasDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
* \brief extract vector according to freenodes
* \param *A pointer to left-top stiffness submatrix
* \param *B pointer to left-bottom stiffness submatrix
* \param *b pointer to the orginal vector
* \param *b1 pointer to the vector related to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \param *uh pointer to the initial solution
* \return void
*/
void extractFreenodesVector2ElasDirichlet(dCSRmat *A, dCSRmat *B, dvector *b, dvector *b1, ELEMENT_DOF *elementDOF, dvector *uh)
{
	int i, j, k, i1, j1;

	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;

	create_dvector(freenodes->row, &b1[1]);

	// achieve b1[0] due to freenodes
	if (uh != NULL)
	{
		for (j1 = 0; j1 < nfreenodes->row; j1++)
		{
			j = nfreenodes->val[j1];
			for (k = B->IA[j]; k < B->IA[j + 1]; k++)
			{
				i = B->JA[k];
				b1[0].val[i] -= B->val[k] * uh->val[j];
			}
		}
	}

	// achieve b1[1] due to freenodes
	for (i1 = 0; i1<b1[1].row; i1++)
	{
		i = freenodes->val[i1];
		b1[1].val[i1] = b->val[i];
	}
}

/**
* \fn void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOF)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A11 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
/*void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOF)
{
	// achiveve A11 due to dirichlet boundary condition
	int i, j, k, l;
	int count;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;


	A11->row = freenodes->row;
	A11->col = A11->row;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	// form A11->IA
	for (i = 0; i<A11->row; i++)
	{
		l = freenodes->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j] == 0)
				A11->IA[i + 1]++;
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++)
	{
		l = freenodes->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j] == 0)
			{
				A11->JA[count] = index->val[j];
				A11->val[count] = A->val[k];
				count++;
			}
		}
	}
}*/

/**
* \fn void extractFreenodesVector(dCSRmat *A, dvector *b, dvector *b1, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, dvector *uh)
* \brief extract vector according to freenodes
* \param *A pointer to the original matrix
* \param *b pointer to the orginal vector
* \param *b1 pointer to the vector related to freenodes
* \param *elementDOFr pointer to relation between elements and DOFs for rows
* \param *elementDOFc pointer to relation between elements and DOFs for columns
* \param *uh pointer to the initial solution
* \return void
*/
void extractFreenodesVector(dCSRmat *A, dvector *b, dvector *b1, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc, dvector *uh)
{
	// achiveve A11 due to dirichlet boundary condition
	int i, j, k, l, i1, j1;
	int count;

	//	ivector *nfFlagr = &elementDOFr->nfFlag;
	ivector *freenodesr = &elementDOFr->freenodes;
	//	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	//	ivector *indexr = &elementDOFr->index;

	ivector *nfFlagc = &elementDOFc->nfFlag;
	ivector *freenodesc = &elementDOFc->freenodes;
	//	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	ivector *indexc = &elementDOFc->index;

	create_dvector(freenodesr->row, b1);

	if (uh == NULL)
		return;

	// achieve b1 due to freenodes
	for (i1 = 0; i1<b1->row; i1++)
	{
		i = freenodesr->val[i1];
		b1->val[i1] = b->val[i];
		for (j1 = A->IA[i]; j1 < A->IA[i + 1]; j1++)
		{
			j = A->JA[j1];
			if (nfFlagc->val[j] == 1)
			{
				b1->val[i1] -= A->val[j1] * uh->val[j];
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A11 pointer to the extracted matrix according to freenodes
* \param *elementDOFr pointer to relation between elements and DOFs for rows
* \param *elementDOFc pointer to relation between elements and DOFs for columns
* \return void
*/
void extractFreenodesMatrix11(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
{
	// achiveve A11 due to dirichlet boundary condition
	int i, j, k, l;
	int count;

	//	ivector *nfFlagr = &elementDOFr->nfFlag;
	ivector *freenodesr = &elementDOFr->freenodes;
	//	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	//	ivector *indexr = &elementDOFr->index;

	ivector *nfFlagc = &elementDOFc->nfFlag;
	ivector *freenodesc = &elementDOFc->freenodes;
	//	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	ivector *indexc = &elementDOFc->index;

	A11->row = freenodesr->row;
	A11->col = freenodesc->row;
	A11->col = A11->row;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	// form A11->IA
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j] == 0)
				A11->IA[i + 1]++;
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j] == 0)
			{
				A11->JA[count] = indexc->val[j];
				A11->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix1r(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
* \brief extract stiffness matrix A1 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void extractFreenodesMatrix1r(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
{
	// achiveve A1 due to dirichlet boundary condition
	int i, j, k, l;
	int count;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	A1->row = freenodes->row;
	A1->col = A->col;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		l = freenodes->val[i];
		A1->IA[i + 1] = A->IA[l + 1] - A->IA[l];
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		l = freenodes->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			A1->JA[count] = A->JA[k];
			A1->val[count] = A->val[k];
			count++;
		}
	}
}

/**
* \fn void extractFreenodesMatrix1c(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
* \brief extract stiffness matrix A1 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void extractFreenodesMatrix1c(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
{
	// achiveve A1 due to dirichlet boundary condition
	int i, j, k;
	int count;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	A1->row = A->row;
	A1->col = freenodes->row;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j] == 0)
				A1->IA[i + 1]++;
		}
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j] == 0)
			{
				A1->JA[count] = index->val[j];
				A1->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
* \brief extract stiffness matrix A1 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void extractFreenodesMatrix1cBlock(dCSRmat *A, dCSRmat *A1, ELEMENT_DOF *elementDOF)
{
	// achiveve A1 according to freenodes
	int i, j, k;
	int count;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	A1->row = A->row;
	A1->col = freenodes->row * 2;
	A1->IA = (int*)calloc(A1->row + 1, sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	int col = A->col / 2;

	// form A1->IA
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j%col] == 0)
				A1->IA[i + 1]++;
		}
	}

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// form A1->JA, A1->val
	A1->JA = (int*)calloc(A1->nnz, sizeof(int));
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A1->row; i++)
	{
		for (k = A->IA[i]; k<A->IA[i + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlag->val[j%col] == 0)
			{
				A1->JA[count] = index->val[j%col] + freenodes->row*(j / col);
				A1->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
* \fn void extractFreenodesMatrix11cBlock(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
* \brief extract stiffness matrix A11 according to freenodes
* \param *A pointer to the original matrix
* \param *A1 pointer to the extracted matrix according to freenodes
* \param *elementDOFr pointer to relation between elements and DOFs in rows
* \param *elementDOFc pointer to relation between elements and DOFs in columns, blockwise
* \return void
*/
void extractFreenodesMatrix11cBlock(dCSRmat *A, dCSRmat *A11, ELEMENT_DOF *elementDOFr, ELEMENT_DOF *elementDOFc)
{
	// achiveve A1 according to freenodes
	int i, j, k, l;
	int count;

	//	ivector *nfFlagr = &elementDOFr->nfFlag;
	ivector *freenodesr = &elementDOFr->freenodes;
	//	ivector *nfreenodesr = &elementDOFr->nfreenodes;
	//	ivector *indexr = &elementDOFr->index;

	ivector *nfFlagc = &elementDOFc->nfFlag;
	ivector *freenodesc = &elementDOFc->freenodes;
	//	ivector *nfreenodesc = &elementDOFc->nfreenodes;
	ivector *indexc = &elementDOFc->index;

	A11->row = freenodesr->row;
	A11->col = freenodesc->row * 2;
	A11->IA = (int*)calloc(A11->row + 1, sizeof(int));
	A11->JA = NULL;
	A11->val = NULL;

	int col = A->col / 2;

	// form A11->IA
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j%col] == 0)
				A11->IA[i + 1]++;
		}
	}

	for (i = 0; i<A11->row; i++)
		A11->IA[i + 1] += A11->IA[i];

	A11->nnz = A11->IA[A11->row];

	// form A11->JA, A11->val
	A11->JA = (int*)calloc(A11->nnz, sizeof(int));
	A11->val = (double*)calloc(A11->nnz, sizeof(double));
	count = 0;
	for (i = 0; i<A11->row; i++)
	{
		l = freenodesr->val[i];
		for (k = A->IA[l]; k<A->IA[l + 1]; k++)
		{
			j = A->JA[k];
			if (nfFlagc->val[j%col] == 0)
			{
				A11->JA[count] = indexc->val[j%col] + freenodesc->row*(j / col);
				A11->val[count] = A->val[k];
				count++;
			}
		}
	}
}

/**
 * \fn int getEdgeDOFsTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \param *row32 2/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
{
	int i,j,l;
	int node;
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	node=elementDOF->val[element][(l+1)%3];
	for(j=0;j<3;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	node=elementDOF->val[element][(l+2)%3];
	for(j=0;j<3;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node=elementDOF->val[element][3+l*(elementDOF->dop-1)+i];
		for(j=0;j<3;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
		}
		count++;
	}
	
	return count;
}

/**
 * \fn int getEdgeDOFsScalarTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int rowstart, int row31, int row32)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \param *row32 2/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsScalarTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int rowstart, int row31, int row32)
{
	int i,l;
	int node, taustart;

	if(elementDOF->dop==0)
	{
		taustart=3*element;
		A->JA[rowstart+count]=taustart;
		A->JA[rowstart+count+row31]=taustart+1;
		A->JA[rowstart+count+row32]=taustart+2;
		count++;
		return count;
	}
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	taustart=3*elementDOF->col*element;
	node=(l+1)%3;
	A->JA[rowstart+count]=taustart+node;
	A->JA[rowstart+count+row31]=taustart+node+elementDOF->col;
	A->JA[rowstart+count+row32]=taustart+node+elementDOF->col*2;
	count++;
	
	node=(l+2)%3;
	A->JA[rowstart+count]=taustart+node;
	A->JA[rowstart+count+row31]=taustart+node+elementDOF->col;
	A->JA[rowstart+count+row32]=taustart+node+elementDOF->col*2;
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node=3+l*(elementDOF->dop-1)+i;
		A->JA[rowstart+count]=taustart+node;
		A->JA[rowstart+count+row31]=taustart+node+elementDOF->col;
		A->JA[rowstart+count+row32]=taustart+node+elementDOF->col*2;
		count++;
	}
	
	return count;
}

/**
 * \fn int getEdgeDOFsVectorTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row31 1/3 of A->IA[i+1] - A->IA[i] 
 * \param *row32 2/3 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsVectorTensor(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row31, int *row32)
{
	int i,j,l;
	int node;

	if(elementDOF->dop==0)
	{
		node=elementDOF->val[element][0];
		for(j=0;j<2;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
		}
		count++;
		return count;
	}
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	node=elementDOF->val[element][(l+1)%3];
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	node=elementDOF->val[element][(l+2)%3];
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
		A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
	}
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node=elementDOF->val[element][3+l*(elementDOF->dop-1)+i];
		for(j=0;j<2;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row31[j]]=node+elementDOF->dof;
			A->JA[rowstart[j]+count+row32[j]]=node+elementDOF->dof*2;
		}
		count++;
	}
	
	return count;
}

/**
 * \fn int getEdgeDOFsVector(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row21)
 * \brief generate JA of A from the DOFs of edge in element 
 * \param *A pointer to stiffness matrix
 * \param count current index of A->JA
 * \param element index of current element
 * \param edge index of current edge
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *rowstart the starting location of i-th row, i.e. A->IA[i] 
 * \param *row21 1/2 of A->IA[i+1] - A->IA[i] 
 * \return count, the current index of A->JA
 */
int getEdgeDOFsVector(dCSRmat *A, int count, int element, int edge, idenmat *elementEdge, ELEMENT_DOF *elementDOF, int *rowstart, int *row21)
{
	int i,j,l;
	int node;
	
	for(l=0;l<elementEdge->col;l++)
	{
		if(elementEdge->val[element][l]==edge)
			break;
	}
	
	node = elementDOF->val[element][(l + 1) % 3];
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row21[j]]=node+elementDOF->dof;
	}
	count++;
	
	node = elementDOF->val[element][(l + 2) % 3];
	for(j=0;j<2;j++)
	{
		A->JA[rowstart[j]+count]=node;
		A->JA[rowstart[j]+count+row21[j]]=node+elementDOF->dof;
	}
	count++;
	
	for(i=0;i<elementDOF->dop-1;i++)
	{
		node = elementDOF->val[element][3 + l*(elementDOF->dop - 1) + i];
		for(j=0;j<2;j++)
		{
			A->JA[rowstart[j]+count]=node;
			A->JA[rowstart[j]+count+row21[j]]=node+elementDOF->dof;
		}
		count++;
	}
	
	return count;
}


/**
 * \fn void patchtogether22(dCSRmat *A11, dCSRmat *A12, dCSRmat *A21, dCSRmat *A22, dCSRmat *A)
 * \brief concatenate submatrices A11, A12, A21, A22 to matrix A by A=[A11 A12 
																	   A21 A22]
 * \param *A11 pointer to the submatrix located at the top-left corner
 * \param *A12 pointer to the submatrix located at the top-right corner
 * \param *A21 pointer to the submatrix located at the bottom-left corner
 * \param *A22 pointer to the submatrix located at the bottom-right corner
 * \param *A pointer to the resulted matrix
 * \return void
 */
void patchtogether22(dCSRmat *A11, dCSRmat *A12, dCSRmat *A21, dCSRmat *A22, dCSRmat *A)
{
	if((A11->row!=A12->row) || (A21->row!=A22->row) || (A11->col!=A21->col) || (A12->col!=A22->col))
	{
		printf("The dimensions of four submatrices do not match in function patchtogether!\n");
		exit(0);
	}

	int i,k,count;
	create_csr_matrix(A11->row+A21->row, A11->col+A12->col, A11->nnz+A12->nnz+A21->nnz+A22->nnz, A);
	
	// form A->IA
	for(i=0;i<A11->row;i++)
		A->IA[i+1]=A11->IA[i+1]-A11->IA[i]+A12->IA[i+1]-A12->IA[i];

	for(i=0;i<A21->row;i++)
		A->IA[A11->row+i+1]=A21->IA[i+1]-A21->IA[i]+A22->IA[i+1]-A22->IA[i];
	
	for(i=0;i<A->row;i++)
		A->IA[i+1]+=A->IA[i];

	// form A->JA, A->val
	count=0;
	for(i=0;i<A11->row;i++)
	{
		for(k=A11->IA[i];k<A11->IA[i+1];k++)
		{
			A->JA[count]=A11->JA[k];
			A->val[count]=A11->val[k];
			count++;
		}
		for(k=A12->IA[i];k<A12->IA[i+1];k++)
		{
			A->JA[count]=A12->JA[k]+A11->col;
			A->val[count]=A12->val[k];
			count++;
		}
	}

	for(i=0;i<A21->row;i++)
	{
		for(k=A21->IA[i];k<A21->IA[i+1];k++)
		{
			A->JA[count]=A21->JA[k];
			A->val[count]=A21->val[k];
			count++;
		}
		for(k=A22->IA[i];k<A22->IA[i+1];k++)
		{
			A->JA[count]=A22->JA[k]+A11->col;
			A->val[count]=A22->val[k];
			count++;
		}
	}	
}

/**
 * \fn void patchtogether33(dCSRmat *A11, dCSRmat *A12, dCSRmat *A13, dCSRmat *A21, dCSRmat *A22, dCSRmat *A23, dCSRmat *A31, dCSRmat *A32, dCSRmat *A33, dCSRmat *A)
 * \brief concatenate submatrices A11, A12, A13, A21, A22, A23, A31, A32, A33 to matrix A by 
 *                                                 A=[A11 A12 A13
 *  												  A21 A22 A23
 *											          A31 A32 A33]
 * \param *A11 pointer to the submatrix located at the top-left corner
 * \param *A12 pointer to the submatrix located at the top-middle corner
 * \param *A13 pointer to the submatrix located at the top-right corner
 * \param *A21 pointer to the submatrix located at the middle-left corner
 * \param *A22 pointer to the submatrix located at the middle-middle corner
 * \param *A23 pointer to the submatrix located at the middle-right corner
 * \param *A31 pointer to the submatrix located at the bottom-left corner
 * \param *A32 pointer to the submatrix located at the bottom-middle corner
 * \param *A33 pointer to the submatrix located at the bottom-right corner
 * \param *A pointer to the resulted matrix
 * \return void
 */
void patchtogether33(dCSRmat *A11, dCSRmat *A12, dCSRmat *A13, dCSRmat *A21, dCSRmat *A22, dCSRmat *A23, dCSRmat *A31, dCSRmat *A32, dCSRmat *A33, dCSRmat *A)
{
	if((A11->row!=A12->row) || (A11->row!=A13->row) || (A21->row!=A22->row) || (A21->row!=A23->row) || (A31->row!=A32->row) || (A31->row!=A33->row))
	{
		printf("The row dimensions of nine submatrices do not match in function patchtogether!\n");
		exit(0);
	}
	
	if((A11->col!=A21->col) || (A11->col!=A31->col) || (A12->col!=A22->col) || (A12->col!=A32->col) || (A13->col!=A23->col) || (A13->col!=A33->col))
	{
		printf("The column dimensions of nine submatrices do not match in function patchtogether!\n");
		exit(0);
	}

	int i,k,count;
	create_csr_matrix(A11->row+A21->row+A31->row, A11->col+A12->col+A13->col, A11->nnz+A12->nnz+A13->nnz+A21->nnz+A22->nnz+A23->nnz+A31->nnz+A32->nnz+A33->nnz, A);
	
	// form A->IA
	for(i=0;i<A11->row;i++)
		A->IA[i+1]=A11->IA[i+1]-A11->IA[i]+A12->IA[i+1]-A12->IA[i]+A13->IA[i+1]-A13->IA[i];

	for(i=0;i<A21->row;i++)
		A->IA[A11->row+i+1]=A21->IA[i+1]-A21->IA[i]+A22->IA[i+1]-A22->IA[i]+A23->IA[i+1]-A23->IA[i];
	
	for(i=0;i<A31->row;i++)
		A->IA[A11->row+A21->row+i+1]=A31->IA[i+1]-A31->IA[i]+A32->IA[i+1]-A32->IA[i]+A33->IA[i+1]-A33->IA[i];

	for(i=0;i<A->row;i++)
		A->IA[i+1]+=A->IA[i];

	// form A->JA, A->val
	count=0;
	for(i=0;i<A11->row;i++)
	{
		for(k=A11->IA[i];k<A11->IA[i+1];k++)
		{
			A->JA[count]=A11->JA[k];
			A->val[count]=A11->val[k];
			count++;
		}
		for(k=A12->IA[i];k<A12->IA[i+1];k++)
		{
			A->JA[count]=A12->JA[k]+A11->col;
			A->val[count]=A12->val[k];
			count++;
		}
		for(k=A13->IA[i];k<A13->IA[i+1];k++)
		{
			A->JA[count]=A13->JA[k]+A11->col+A12->col;
			A->val[count]=A13->val[k];
			count++;
		}
	}

	for(i=0;i<A21->row;i++)
	{
		for(k=A21->IA[i];k<A21->IA[i+1];k++)
		{
			A->JA[count]=A21->JA[k];
			A->val[count]=A21->val[k];
			count++;
		}
		for(k=A22->IA[i];k<A22->IA[i+1];k++)
		{
			A->JA[count]=A22->JA[k]+A11->col;
			A->val[count]=A22->val[k];
			count++;
		}
		for(k=A23->IA[i];k<A23->IA[i+1];k++)
		{
			A->JA[count]=A23->JA[k]+A11->col+A12->col;
			A->val[count]=A23->val[k];
			count++;
		}
	}

	for(i=0;i<A31->row;i++)
	{
		for(k=A31->IA[i];k<A31->IA[i+1];k++)
		{
			A->JA[count]=A31->JA[k];
			A->val[count]=A31->val[k];
			count++;
		}
		for(k=A32->IA[i];k<A32->IA[i+1];k++)
		{
			A->JA[count]=A32->JA[k]+A11->col;
			A->val[count]=A32->val[k];
			count++;
		}
		for(k=A33->IA[i];k<A33->IA[i+1];k++)
		{
			A->JA[count]=A33->JA[k]+A11->col+A12->col;
			A->val[count]=A33->val[k];
			count++;
		}
	}
}


/**
 * \fn void getPenaltyParameters(ddenmat *etas, idenmat *elementEdge, EDGE *edges, double *alpha, double *beta)
 * \brief get penalty parameters for HDG method 
 * \param *etas pointer to penalty parameter
  * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param alpha power of penalty parameter
 * \param beta coefficient of penalty parameter
 * \return void
 */
void getPenaltyParameters(ddenmat *etas, idenmat *elementEdge, EDGE *edges, double *alpha, double *beta)
{
	int i,j,k,l;
	
	for(i=0;i<3;i++)
		create_dden_matrix(elementEdge->row, elementEdge->col, etas+i);

	for(k=0;k<elementEdge->row;k++)
	{
		for(j=0;j<elementEdge->col;j++)
		{
			l=elementEdge->val[k][j];
			for(i=0;i<3;i++)
				etas[i].val[k][j] = pow(edges->length[l], alpha[i])*beta[i];
		}
	}	
}
