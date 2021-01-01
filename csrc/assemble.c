/*
 *  assemble.c
 *  DGFEM
 *
 *  Created by Xuehai Huang on 10/30/08.
 *  Modified by Xuehai Huang on 3/21/12.
 *  Copyright 2012 WZU. All rights reserved.
 *
 */

/*! \file assemble.c
 *  \brief Assembling for stiffness matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"

/**
 * \fn int getmesh(int domain_num, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum)
 * \brief generate mesh information
 * \param domain_num number of domain
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices, the last 3 columns store the indexes of edge's midpoints
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                       the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
                                       the third column stores the Dirichlet status of the points(0: nondirichlet, -1: dirichlet)
 * \param *edgesTran pointer to the tranpose of edges, used to get restriction operator. The relation between nodes and edges. JA stores edge index, A stores another vertex
 * \param *nodeCEdge record the index of coarse edge which the node belong to; if the node is located in the coarset grid, it will be set -1
 * \param levelNum total level number of grid
 * \return 1 if succeed 0 if fail
 */
int getmesh(int domain_num, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, iCSRmat *edgesTran, ivector *nodeCEdge, int levelNum)
{	
	int IsExist=getCoarseInfo(domain_num, &nodes[0], &elements[0], &edges[0], &edgesTran[0], nodeCEdge);
	if(IsExist==0)
	{
		printf("Constructing coarse grid fails!\n");
		return 0;
	}

	int i,j,k,l;
	int nvertices;

	for(l=0;l<levelNum-1;l++)
		uniformrefine(&nodes[l], &elements[l], &edges[l], &nodes[l+1], &elements[l+1], &edges[l+1], &edgesTran[l+1], nodeCEdge);
		
	// generate elementEdge	
	int point1, point2, element1, element2;
	
	for (l = 0; l < levelNum; l++)
	{
		create_iden_matrix(elements[l].row, 3, elementEdge + l);

		for (i = 0; i<edges[l].row; i++)
		{
			point1 = edges[l].val[i][0];
			point2 = edges[l].val[i][1];
			element1 = edges[l].val[i][2];
			element2 = edges[l].val[i][3];

			for (j = 0; j<3; j++)
			{
				if (elements[l].val[element1][j] == point1)
					break;
			}
			for (k = 0; k<3; k++)
			{
				if (elements[l].val[element1][k] == point2)
					break;
			}
			elementEdge[l].val[element1][3 - j - k] = i;

			if (element2>-1)
			{
				for (j = 0; j<3; j++)
				{
					if (elements[l].val[element2][j] == point1)
						break;
				}
				for (k = 0; k<3; k++)
				{
					if (elements[l].val[element2][k] == point2)
						break;
				}
				elementEdge[l].val[element2][3 - j - k] = i;
			}
		}
	}
	// end generate elementEdge
		
	return 1;
}


/**
 * \fn void getElementDOF1D(ELEMENT_DOF *elementDOF, int ne, int dop)
 * \brief get the degrees of freedom of piecewise Lagrange element on edges
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param ne number of edges
 * \param dop degree of polynomial
 */
void getElementDOF1D(ELEMENT_DOF *elementDOF, int ne, int dop)
{
	int i,j;
	int count;
	if(dop<0)
	{
		elementDOF=NULL;
		return;
	}

	create_elementDOF(dop, ne*(dop+1), ne, dop+1, elementDOF);

	for(j=0;j<ne;j++)
	{
		count=j*elementDOF->col;
		for(i=0;i<elementDOF->col;i++)
			elementDOF->val[j][i]=count+i;
	}

}

/**
 * \fn void getElementDOF1D_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
 * \brief get the degrees of freedom of Lagrange element on edges
 * \param *edgeDOF pointer to relation between edges and DOFs
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF1D_Continue(ELEMENT_DOF *edgeDOF, ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn + ne*(dop-1), ne, dop+1, edgeDOF);

	int node, edge;
	for(j=0;j<ne;j++)
	{
		for(i=0;i<2;i++)
		{
			node=edges->val[j][i];
			edgeDOF->val[j][i]=node;
		}
		
		for(i=0;i<dop-1;i++)
			edgeDOF->val[j][2+i] = nn + j*(dop-1) + i;
	}

	create_elementDOF(dop, nn + ne*(dop-1), nt, dop*3, elementDOF);

	int orient;
	for(k=0;k<nt;k++)
	{
		for(i=0;i<3;i++)
		{
			node=elements->val[k][i];
			elementDOF->val[k][i]=node;
		}
		
		for(j=0;j<3;j++)
		{
			edge=elementEdge->val[k][j];
			if(elements->val[k][(j+1)%3]== edges->val[edge][0])
				orient=1;
			else
				orient=0;
			
			if(orient==1)
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+i] = nn + edge*(dop-1) + i;
			}
			else
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+dop-2-i] = nn + edge*(dop-1) + i;
			}
		}
	}
}

/**
 * \fn void getElementDOF(ELEMENT_DOF *elementDOF, int nt, int dop)
 * \brief get the degrees of freedom of piecewise Lagrange element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param nt number of elements
 * \param dop degree of polynomial
 */
void getElementDOF(ELEMENT_DOF *elementDOF, int nt, int dop)
{
	int i,k;
	int count;
	if(dop<1)
		dop=0;

	create_elementDOF(dop, nt*(dop+1)*(dop+2)/2, nt, (dop+1)*(dop+2)/2, elementDOF);

	for(k=0;k<nt;k++)
	{
		count=k*elementDOF->col;
		for(i=0;i<elementDOF->col;i++)
			elementDOF->val[k][i]=count+i;
	}
}

/**
* \fn void getElementDOF_symTensor(ELEMENT_DOF *elementDOF, int nt, int dop)
* \brief get the degrees of freedom of symmtric tensor of piecewise Lagrange element
* \param *elementDOF pointer to relation between elements and DOFs
* \param nt number of elements
* \param dop degree of polynomial
*/
void getElementDOF_symTensor(ELEMENT_DOF *elementDOF, int nt, int dop)
{
	int i, k;
	int count;
	if (dop<1)
		dop = 0;

	create_elementDOF(dop, 3 * nt*(dop + 1)*(dop + 2) / 2, nt, 3 * (dop + 1)*(dop + 2) / 2, elementDOF);

	for (k = 0; k<nt; k++)
	{
		count = k*elementDOF->col;
		for (i = 0; i<elementDOF->col; i++)
			elementDOF->val[k][i] = count + i;
	}
}

/**
 * \fn void getElementDOF_Lagrange(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int nedges, int nvertices, int dop)
 * \brief get the degrees of freedom of Lagrange element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_Lagrange(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn + ne*(dop-1) + nt*(dop-1)*(dop-2)/2, nt, (dop+1)*(dop+2)/2, elementDOF);

	int node, edge;
	int orient;
	for(k=0;k<nt;k++)
	{
		for(i=0;i<3;i++)
		{
			node=elements->val[k][i];
			elementDOF->val[k][i]=node;
		}
		
		for(j=0;j<3;j++)
		{
			edge=elementEdge->val[k][j];
			if(elements->val[k][(j+1)%3]== edges->val[edge][0])
				orient=1;
			else
				orient=0;
			
			if(orient==1)
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+i] = nn + edge*(dop-1) + i;
			}
			else
			{
				for(i=0;i<dop-1;i++)
					elementDOF->val[k][3+(dop-1)*j+dop-2-i] = nn + edge*(dop-1) + i;
			}
		}

		for(i=0;i<(dop-1)*(dop-2)/2;i++)
		{
			elementDOF->val[k][3*dop+i] = nn + ne*(dop-1) + k*(dop-1)*(dop-2)/2 + i;
		}
	}
}

/**
* \fn void getElementDOF_NoncfmP1(ELEMENT_DOF *elementDOF, idenmat *elementEdge, int ne)
* \brief get the degrees of freedom of nonconforming P1 element
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param ne number of edges
*/
void getElementDOF_NoncfmP1(ELEMENT_DOF *elementDOF, idenmat *elementEdge, int ne)
{
	int i, k;
	int nt = elementEdge->row;

	create_elementDOF(1, ne, nt, 3, elementDOF);

	for (k = 0; k<nt; k++)
	{
		for (i = 0; i<3; i++)
			elementDOF->val[k][i] = elementEdge->val[k][i];
	}
}

/**
* \fn void getElementDOF_Morley(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int ne, int nvertices)
* \brief get the degrees of freedom of Morley element
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param ne number of edges
* \param nvertices number of vertices
*/
void getElementDOF_Morley(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int ne, int nvertices)
{
	int i, j, k;
	int nt = elements->row;
	int nn = nvertices;

	create_elementDOF(2, nn + ne, nt, 6, elementDOF);

	int node, edge;
	for (k = 0; k<nt; k++)
	{
		for (i = 0; i<3; i++)
		{
			node = elements->val[k][i];
			elementDOF->val[k][i] = node;
		}

		for (j = 0; j<3; j++)
		{
			edge = elementEdge->val[k][j];
			elementDOF->val[k][3 + j] = nn + edge;
		}
	}
}

/**
* \fn void getElementDOF_MINI(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices)
* \brief get the degrees of freedom of MINI element for Stokes equation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
* \param nvertices number of vertices
*/
void getElementDOF_MINI(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices)
{
	int i, j, k;
	int nt = elements->row;
	int nn = nvertices;

	create_elementDOF(1, nn + nt, nt, 4, elementDOF);

	int node, edge;
	for (k = 0; k<nt; k++)
	{
		for (i = 0; i<3; i++)
		{
			node = elements->val[k][i];
			elementDOF->val[k][i] = node;
		}

		elementDOF->val[k][3] = nn + k;
	}
}

/**
* \fn void getElementDOF_MINIsymTensor(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices)
* \brief get the degrees of freedom of symmtric stress of MINI element for Stokes equation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
* \param nvertices number of vertices
*/
void getElementDOF_MINIsymTensor(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices)
{
	int i, k;
	int nt = elements->row;
	int nn = nvertices;

	create_elementDOF(1, nn + nt * 5, nt, 8, elementDOF);

	int node;
	for (k = 0; k<nt; k++)
	{
		for (i = 0; i<3; i++)
		{
			node = elements->val[k][i];
			elementDOF->val[k][i] = node;
		}

		for (i = 0; i<5; i++)
			elementDOF->val[k][3 + i] = nn + k * 5 + i;
	}
}

/**
* \fn void getElementDOF_MINIsymTensorP1P0(ELEMENT_DOF *elementDOF, ELEMENT *elements, int nvertices)
* \brief get the degrees of freedom of symmtric stress of MINI element for Stokes equation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
* \param nvertices number of vertices
*/
void getElementDOF_MINIsymTensorP1P0(ELEMENT_DOF *elementDOF, int nt)
{
	int i, k;
	int count;

	create_elementDOF(1, nt * 9, nt, 9, elementDOF);

	int node;
	for (k = 0; k<nt; k++)
	{
		count = k*elementDOF->col;
		for (i = 0; i<elementDOF->col; i++)
			elementDOF->val[k][i] = count + i;
	}
}

/**
* \fn void getElementDOF_CR(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int ne, int nvertices)
* \brief get the degrees of freedom of Crouzeix每Raviart element for Stokes equation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param ne number of edges
* \param nvertices number of vertices
*/
void getElementDOF_CR(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, int ne, int nvertices)
{
	int i, j, k;
	int nt = elements->row;
	int nn = nvertices;

	create_elementDOF(2, nn + ne + nt, nt, 7, elementDOF);

	int node, edge;
	for (k = 0; k<nt; k++)
	{
		for (i = 0; i<3; i++)
		{
			node = elements->val[k][i];
			elementDOF->val[k][i] = node;
		}

		for (j = 0; j<3; j++)
		{
			edge = elementEdge->val[k][j];
			elementDOF->val[k][3 + j] = nn + edge;
		}

		elementDOF->val[k][6] = nn + ne + k;
	}
}

/**
* \fn void getElementDOF_CRsymTensor(ELEMENT_DOF *elementDOF, int nt)
* \brief get the degrees of freedom of symmtric stress of Crouzeix每Raviart element for Stokes equation
* \param *elementDOF pointer to relation between elements and DOFs
* \param nt number of elements
* \param dop degree of polynomial
*/
void getElementDOF_CRsymTensor(ELEMENT_DOF *elementDOF, int nt)
{
	int i, k;
	int count;

	create_elementDOF(1, nt * 11, nt, 11, elementDOF);

	for (k = 0; k<nt; k++)
	{
		count = k*elementDOF->col;
		for (i = 0; i<elementDOF->col; i++)
			elementDOF->val[k][i] = count + i;
	}
}


/**
 * \fn void getElementDOF_HuZhang(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
 * \brief get the degrees of freedom of Hu-Zhang element
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param nvertices number of vertices
 * \param dop degree of polynomial
 */
void getElementDOF_HuZhang(ELEMENT_DOF *elementDOF, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, int nvertices, int dop)
{
	int i,j,k;
	int nt=elements->row;
	int ne=edges->row;
	int nn=nvertices;
	if(dop<2)
		dop=1;

	create_elementDOF(dop, nn*3 + ne*(dop-1)*2 + nt*((dop-1)*3 + (dop-1)*(dop-2)/2*3), nt, (dop+1)*(dop+2)/2*3, elementDOF);

	int node, edge;
	int orient;
	for(k=0;k<nt;k++)
	{
		for(i=0;i<3;i++)
		{
			node=elements->val[k][i];
			elementDOF->val[k][i]=node;
			elementDOF->val[k][i+3]=node+nn;
			elementDOF->val[k][i+6]=node+nn*2;
		}
		
		for(j=0;j<3;j++)
		{
			edge=elementEdge->val[k][j];
			if(elements->val[k][(j+1)%3]== edges->val[edge][0])
				orient=1;
			else
				orient=0;
			
			if(orient==1)
			{
				for(i=0;i<dop-1;i++)
				{
					elementDOF->val[k][9+(dop-1)*j+i] = nn*3 + edge*(dop-1) + i;
					elementDOF->val[k][9+(dop-1)*3+(dop-1)*j+i] = nn*3 + ne*(dop-1) + edge*(dop-1) + i;
				}
			}
			else
			{
				for(i=0;i<dop-1;i++)
				{
					elementDOF->val[k][9+(dop-1)*j+dop-2-i] = nn*3 + edge*(dop-1) + i;
					elementDOF->val[k][9+(dop-1)*3+(dop-1)*j+dop-2-i] = nn*3 + ne*(dop-1) + edge*(dop-1) + i;
				}
			}

			for(i=0;i<dop-1;i++)
				elementDOF->val[k][9+(dop-1)*6+(dop-1)*j+i] = nn*3 + ne*(dop-1)*2 + k*((dop-1)*3 + (dop-1)*(dop-2)/2*3) + (dop-1)*j + i;
	//			elementDOF->val[k][9+(dop-1)*6+(dop-1)*j+i] = nn*3 + ne*(dop-1)*2 + k*(dop-1)*3 + (dop-1)*j + i;
		} // j

		for(i=0;i<(dop-1)*(dop-2)/2*3;i++)
			elementDOF->val[k][9*dop+i] = nn*3 + ne*(dop-1)*2 + k*((dop-1)*3 + (dop-1)*(dop-2)/2*3) + (dop-1)*3 + i;
	
/*		for(i=0;i<(dop-1)*(dop-2)/2;i++)
		{
			elementDOF->val[k][9*dop+i] = nn*3 + ne*(dop-1)*2 + nt*(dop-1)*3 + k*(dop-1)*(dop-2)/2 + i;
			elementDOF->val[k][9*dop+(dop-1)*(dop-2)/2+i] = nn*3 + ne*(dop-1)*2 + nt*((dop-1)*3 + (dop-1)*(dop-2)/2) + k*(dop-1)*(dop-2)/2 + i;
			elementDOF->val[k][9*dop+(dop-1)*(dop-2)+i] = nn*3 + ne*(dop-1)*2 + nt*((dop-1)*3 + (dop-1)*(dop-2)) + k*(dop-1)*(dop-2)/2 + i;
		}*/
	} // k
}

void assemble(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	// void for adaptive fem to run correctly
}

/**
* \fn void assemble_kirchhoffplateMorley(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double nu)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemble_kirchhoffplateMorley(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double nu)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = b
	Bx1 + Cx2   = 0
	where x1: u_h, x2: p_h
	**/
	dCSRmat A;
	dvector b;

	assembleStiffmatrixKirchhoffplateMorley2d(&A, &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, nu);

	// initial solution
	create_dvector(b.row, uh);

    //extract
	extractFreenodesVector(&A, &b, ptr_b, elementDOF, elementDOF, uh);
	free_dvector(&b);
	extractFreenodesMatrix11(&A, ptr_A, elementDOF, elementDOF);
	free_csr_matrix(&A);
}

/**
* \fn void assembleStiffmatrixKirchhoffplateMorley2d(dCSRmat *A, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double nu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrixKirchhoffplateMorley2d(dCSRmat *A, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double nu)
{
	int i, j, k, l;

	A->row = elementDOF->dof;
	A->col = A->row;
	A->IA = (int*)malloc((A->row + 1) * sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node, curnode[2];

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints(0, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	A->IA[0] = 0;
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)malloc(A->nnz * sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i+1]; j++)
		{
	/*		A->JA[A->IA[i]+A->IA[i + 1]-j-1] = istart;
			istart = index[istart];
			index[A->JA[A->IA[i] + A->IA[i + 1] - j - 1]] = -1;*/
			A->JA[j] = istart;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							morley_basis2(lambdas, s, elen, eta, xi, nv, nve, k1, phi1);
							morley_basis2(lambdas, s, elen, eta, xi, nv, nve, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * ((1 - nu)*(phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) + nu*(phi1[0] + phi1[1])*(phi2[0] + phi2[1]));
						}
						break;
					}
				} // j1
			} // k2		
		} // k1
	} // k

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				morley_basis(lambdas, s, elen, eta, xi, nv, nve, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b->val[i] += 2 * s*gauss0[i1][2] * f(x, y)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemble_poissonMorley(dCSRmat *ptr_A1, dvector *ptr_A2,  dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemble_poissonMorley(dCSRmat *ptr_A1, dvector *ptr_A2,  dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	A1x1      = b1
	    A2x2  = b2
	where x1: u_h, x2: p_h
	**/
	dCSRmat A1;
	dvector A2, b[2];

	assembleStiffmatrixPoissonMorley2d(&A1, &A2, b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran);

	// initial solution
	create_dvector(b[0].row, &uh[0]);
	create_dvector(b[1].row, &uh[1]);

	//extract
	extractFreenodesVector(&A1, &b[0], &ptr_b[0], &elementDOF[0], &elementDOF[0], &uh[0]);
	free_dvector(&b[0]);
	extractFreenodesMatrix11(&A1, ptr_A1, &elementDOF[0], &elementDOF[0]);
	free_csr_matrix(&A1);
	int i1, i;
	ivector *freenodes = &elementDOF[1].freenodes;
	create_dvector(freenodes->row, ptr_A2);
	create_dvector(freenodes->row, &ptr_b[1]);
	for (i1 = 0; i1 < ptr_b[1].row; i1++)
	{
		i = freenodes->val[i1];
		ptr_A2->val[i1] = A2.val[i];
		ptr_b[1].val[i1] = b[1].val[i];
	}
	free_dvector(&A2);
	free_dvector(&b[1]);
}

/**
* \fn void assembleStiffmatrixPoissonMorley2d(dCSRmat *A1, dvector *A2, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrixPoissonMorley2d(dCSRmat *A1, dvector *A2, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran)
{
	int i, j, k, l;

	A1->row = elementDOF[0].dof;
	A1->col = A1->row;
	A1->IA = (int*)malloc((A1->row + 1) * sizeof(int));
	A1->JA = NULL;
	A1->val = NULL;

	create_dvector(elementDOF[1].dof, A2);
	create_dvector(A1->row, &b[0]);
	create_dvector(A2->row, &b[1]);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node, curnode[2];

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A1->col * sizeof(int));
	for (i = 0; i<A1->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	A1->IA[0] = 0;
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A1->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A1->row; i++)
		A1->IA[i + 1] += A1->IA[i];

	A1->nnz = A1->IA[A1->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A1->JA = (int*)malloc(A1->nnz * sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A1->IA[i]; j<A1->IA[i + 1]; j++)
		{
			/*		A->JA[A->IA[i]+A->IA[i + 1]-j-1] = istart;
			istart = index[istart];
			index[A->JA[A->IA[i] + A->IA[i + 1] - j - 1]] = -1;*/
			A1->JA[j] = istart;
			istart = index[istart];
			index[A1->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A1->val = (double*)calloc(A1->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		// A1
		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A1->IA[i]; j1<A1->IA[i + 1]; j1++)
				{
					if (A1->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							morley_basis1(lambdas, s, elen, eta, xi, nv, nve, k1, phi1);
							morley_basis1(lambdas, s, elen, eta, xi, nv, nve, k2, phi2);
							A1->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2	
		} // k1

		// A2
		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp21; i1++)
			{
				lambdas[0] = gauss21[i1][0];
				lambdas[1] = gauss21[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + k1, phi1);
				A2->val[i] += 2 * s*gauss21[i1][2] * (phi1[0] * phi1[0] + phi1[1] * phi1[1]);
			}
		}
	} // k

	  /************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				morley_basis(lambdas, s, elen, eta, xi, nv, nve, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b[0].val[i] += 2 * s*gauss0[i1][2] * f(x, y)*phi;
			} // i1
		} // k1 

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				morley_basis(lambdas, s, elen, eta, xi, nv, nve, 3 + k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b[1].val[i] += 2 * s*gauss0[i1][2] * f(x, y)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemble_kirchhoffplate_StokesNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemble_kirchhoffplate_StokesNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = b
	Bx1 + Cx2   = 0
	where x1: u_h, x2: p_h
	**/
	dCSRmat A, B;
	dvector b;

	if (problem_num == 1)
		assembleStiffmatrix_biharmonic_StokesNcP1P02d(&A, &B, &ptr_A[3], &b, &ptr_b[1], uhpm, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, elementDOFpm);
	else
	{
		assembleStiffmatrix_kirchhoffplate_StokesNcP1P02d(&A, &B, &ptr_A[3], &b, &ptr_b[1], uhpm, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, elementDOFpm, nu);
		// p / (1-nu) as the new stress for the consideration of solver
		int i;
/*		for (i = 0; i < A.nnz; i++)
			A.val[i] /= (1 - nu);
		for (i = 0; i < b.row; i++)
			b.val[i] /= (1 - nu);*/
	}
	// initial solution
	create_dvector(b.row, &uh[0]);
	create_dvector(ptr_b[1].row, &uh[1]);
	
	//extract
	extractFreenodesVector2StokesDirichlet(&A, &B, &b, ptr_b, &elementDOF[0], &uh[0]);
	free_dvector(&b);
	extractFreenodesMatrix11(&A, &ptr_A[0], &elementDOF[0], &elementDOF[0]);
	free_csr_matrix(&A);
	extractFreenodesMatrix1c(&B, &ptr_A[2], &elementDOF[0]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrix_biharmonic_StokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_biharmonic_StokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	/*BT->row = A->row;
	BT->col = elementDOF[1].dof;
	BT->IA = (int*)calloc(BT->row + 1, sizeof(int));
	BT->JA = NULL;
	BT->val = NULL; */

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints((elementDOF[0].dop - 1) * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp23 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
	double gauss23[num_qp23][3];
	init_Gauss(num_qp23, 2, gauss23); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;
		A->IA[i + 1 + elementDOF[0].dof] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
			A->JA[j - A->IA[i] + A->IA[i + elementDOF[0].dof]] = istart + elementDOF[0].dof;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	for (i = 0; i < A->IA[elementDOF[0].dof]; i++)
		A->val[i + A->IA[elementDOF[0].dof]] = A->val[i];

	/************************************************** stiffness matrix B *****************************************************************/
	// step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i < elementDOF[1].dof; i++)
		B->IA[i + 1] += elementDOF[0].col * 2;

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		rowstart[0] = B->IA[elementDOF[1].val[k][0]];
		for (j = 0; j < elementDOF[0].col; j++)
		{
			B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			B->JA[rowstart[0] + j + elementDOF[0].col] = elementDOF[0].val[k][j] + elementDOF[0].dof;
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		i = elementDOF[1].val[k][0];
		row21[0] = elementDOF[0].col;
		for (k2 = 0; k2<elementDOF[0].col; k2++)
		{
			j = elementDOF[0].val[k][k2];
			// b11
			for (j1 = B->IA[i]; j1<B->IA[i] + row21[0]; j1++)
			{
				if (B->JA[j1] == j)
				{
					for (i1 = 0; i1<num_qp22; i1++)
					{
						lambdas[0] = gauss22[i1][0];
						lambdas[1] = gauss22[i1][1];
						lambdas[2] = 1 - lambdas[0] - lambdas[1];
						ncp1_basis1(lambdas, s, eta, xi, k2, phi1);
						B->val[j1] += 2 * s*gauss22[i1][2] * phi1[0];
						B->val[j1 + row21[0]] += 2 * s*gauss22[i1][2] * phi1[1];
					}
					break;
				}
				if (j1 == (B->IA[i] + row21[0]))
				{
					printf("There is something wrong in constructing b11 of B\n");
					exit(1);
				}
			}
		} // k1
	} // k

	/************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				ncp1_basis(lambdas, k1, &phi);
				for (j = 0; j<elementDOFpm[0].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, j, phi1);
					j1 = elementDOFpm[0].val[k][j];
					temp[0] += phi1[1] * uhpm[0].val[j1];
					temp[1] -= phi1[0] * uhpm[0].val[j1];
				}
				for (j = 0; j<elementDOFpm[1].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3+j, phi1);
					j1 = elementDOFpm[1].val[k][j];
					temp[0] += phi1[1] * uhpm[1].val[j1];
					temp[1] -= phi1[0] * uhpm[1].val[j1];
				}
				b1->val[i] += 2 * s*gauss0[i1][2] * temp[0] * phi;
				b1->val[i + elementDOF[0].dof] += 2 * s*gauss0[i1][2] * temp[1] * phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStiffmatrix_kirchhoffplate_StokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_kirchhoffplate_StokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	/*BT->row = A->row;
	BT->col = elementDOF[1].dof;
	BT->IA = (int*)calloc(BT->row + 1, sizeof(int));
	BT->JA = NULL;
	BT->val = NULL; */

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count, curnode[2];

	int num_qp21 = getNumQuadPoints((elementDOF[0].dop - 1) * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp23 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
	double gauss23[num_qp23][3];
	init_Gauss(num_qp23, 2, gauss23); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[i + 1 + elementDOF[0].dof] = count * 2;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		row21[0] = (A->IA[i + 1] - A->IA[i]) / 2;
		rowstart[1] = A->IA[i + elementDOF[0].dof];
		for (j = A->IA[i]; j<A->IA[i] + row21[0]; j++)
		{
			A->JA[j] = istart;
			A->JA[j + row21[0]] = istart + elementDOF[0].dof;
			A->JA[j - A->IA[i] + rowstart[1]] = A->JA[j];
			A->JA[j - A->IA[i] + rowstart[1] + row21[0]] = A->JA[j + row21[0]];
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			curnode[0] = elementDOF[0].val[k][k1];
			curnode[1] = curnode[0] + elementDOF[0].dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * ((1 - nu)*phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
				  // a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF[0].dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi1[1] * phi2[0] * nu;
						}
						break;
					}
				} // j1
				  // a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] * nu;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF[0].dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + (1 - nu)*phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k


	/************************************************** stiffness matrix B *****************************************************************/
	// step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i < elementDOF[1].dof; i++)
		B->IA[i + 1] += elementDOF[0].col * 2;

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		rowstart[0] = B->IA[elementDOF[1].val[k][0]];
		for (j = 0; j < elementDOF[0].col; j++)
		{
			B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			B->JA[rowstart[0] + j + elementDOF[0].col] = elementDOF[0].val[k][j] + elementDOF[0].dof;
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		i = elementDOF[1].val[k][0];
		row21[0] = elementDOF[0].col;
		for (k2 = 0; k2<elementDOF[0].col; k2++)
		{
			j = elementDOF[0].val[k][k2];
			// b11
			for (j1 = B->IA[i]; j1<B->IA[i] + row21[0]; j1++)
			{
				if (B->JA[j1] == j)
				{
					for (i1 = 0; i1<num_qp22; i1++)
					{
						lambdas[0] = gauss22[i1][0];
						lambdas[1] = gauss22[i1][1];
						lambdas[2] = 1 - lambdas[0] - lambdas[1];
						ncp1_basis1(lambdas, s, eta, xi, k2, phi1);
						B->val[j1] += 2 * s*gauss22[i1][2] * phi1[0];
						B->val[j1 + row21[0]] += 2 * s*gauss22[i1][2] * phi1[1];
					}
					break;
				}
				if (j1 == (B->IA[i] + row21[0]))
				{
					printf("There is something wrong in constructing b11 of B\n");
					exit(1);
				}
			}
		} // k1
	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				ncp1_basis(lambdas, k1, &phi);
				for (j = 0; j<elementDOFpm[0].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, j, phi1);
					j1 = elementDOFpm[0].val[k][j];
					temp[0] += phi1[1] * uhpm[0].val[j1];
					temp[1] -= phi1[0] * uhpm[0].val[j1];
				}
				for (j = 0; j<elementDOFpm[1].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + j, phi1);
					j1 = elementDOFpm[1].val[k][j];
					temp[0] += phi1[1] * uhpm[1].val[j1];
					temp[1] -= phi1[0] * uhpm[1].val[j1];
				}
				b1->val[i] += 2 * s*gauss0[i1][2] * temp[0] * phi;
				b1->val[i + elementDOF[0].dof] += 2 * s*gauss0[i1][2] * temp[1] * phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemble_kirchhoffplate_StokesrotNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemble_kirchhoffplate_StokesrotNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = b
	Bx1 + Cx2   = 0
	where x1: u_h, x2: p_h
	**/
	dCSRmat A, B;
	dvector b;

	if (problem_num == 1)
		assembleStiffmatrix_biharmonic_StokesrotNcP1P02d(&A, &B, &ptr_A[3], &b, &ptr_b[1], uhpm, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, elementDOFpm);
	else
	{
		assembleStiffmatrix_kirchhoffplate_StokesrotNcP1P02d(&A, &B, &ptr_A[3], &b, &ptr_b[1], uhpm, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, elementDOFpm, nu);
		// p / (1-nu) as the new stress for the consideration of solver
		int i;
		/*		for (i = 0; i < A.nnz; i++)
		A.val[i] /= (1 - nu);
		for (i = 0; i < b.row; i++)
		b.val[i] /= (1 - nu);*/
	}
	// initial solution
	create_dvector(b.row, &uh[0]);
	create_dvector(ptr_b[1].row, &uh[1]);

	//extract
	extractFreenodesVector2StokesDirichlet(&A, &B, &b, ptr_b, &elementDOF[0], &uh[0]);
	free_dvector(&b);
	extractFreenodesMatrix11(&A, &ptr_A[0], &elementDOF[0], &elementDOF[0]);
	free_csr_matrix(&A);
	extractFreenodesMatrix1c(&B, &ptr_A[2], &elementDOF[0]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrix_biharmonic_StokesrotNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_biharmonic_StokesrotNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	/*BT->row = A->row;
	BT->col = elementDOF[1].dof;
	BT->IA = (int*)calloc(BT->row + 1, sizeof(int));
	BT->JA = NULL;
	BT->val = NULL; */

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints((elementDOF[0].dop - 1) * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp23 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
	double gauss23[num_qp23][3];
	init_Gauss(num_qp23, 2, gauss23); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;
		A->IA[i + 1 + elementDOF[0].dof] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
			A->JA[j - A->IA[i] + A->IA[i + elementDOF[0].dof]] = istart + elementDOF[0].dof;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	for (i = 0; i < A->IA[elementDOF[0].dof]; i++)
		A->val[i + A->IA[elementDOF[0].dof]] = A->val[i];

	/************************************************** stiffness matrix B *****************************************************************/
	// step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i < elementDOF[1].dof; i++)
		B->IA[i + 1] += elementDOF[0].col * 2;

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		rowstart[0] = B->IA[elementDOF[1].val[k][0]];
		for (j = 0; j < elementDOF[0].col; j++)
		{
			B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			B->JA[rowstart[0] + j + elementDOF[0].col] = elementDOF[0].val[k][j] + elementDOF[0].dof;
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		i = elementDOF[1].val[k][0];
		row21[0] = elementDOF[0].col;
		for (k2 = 0; k2<elementDOF[0].col; k2++)
		{
			j = elementDOF[0].val[k][k2];
			// b11
			for (j1 = B->IA[i]; j1<B->IA[i] + row21[0]; j1++)
			{
				if (B->JA[j1] == j)
				{
					for (i1 = 0; i1<num_qp22; i1++)
					{
						lambdas[0] = gauss22[i1][0];
						lambdas[1] = gauss22[i1][1];
						lambdas[2] = 1 - lambdas[0] - lambdas[1];
						ncp1_basis1(lambdas, s, eta, xi, k2, phi1);
						B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[1];
						B->val[j1 + row21[0]] += 2 * s*gauss22[i1][2] * phi1[0];
					}
					break;
				}
				if (j1 == (B->IA[i] + row21[0]))
				{
					printf("There is something wrong in constructing b11 of B\n");
					exit(1);
				}
			}
		} // k1
	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				ncp1_basis(lambdas, k1, &phi);
				for (j = 0; j<elementDOFpm[0].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, j, phi1);
					j1 = elementDOFpm[0].val[k][j];
					temp[0] += phi1[0] * uhpm[0].val[j1];
					temp[1] += phi1[1] * uhpm[0].val[j1];
				}
				for (j = 0; j<elementDOFpm[1].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + j, phi1);
					j1 = elementDOFpm[1].val[k][j];
					temp[0] += phi1[0] * uhpm[1].val[j1];
					temp[1] += phi1[1] * uhpm[1].val[j1];
				}
				b1->val[i] += 2 * s*gauss0[i1][2] * temp[0] * phi;
				b1->val[i + elementDOF[0].dof] += 2 * s*gauss0[i1][2] * temp[1] * phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStiffmatrix_kirchhoffplate_StokesrotNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_kirchhoffplate_StokesrotNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	/*BT->row = A->row;
	BT->col = elementDOF[1].dof;
	BT->IA = (int*)calloc(BT->row + 1, sizeof(int));
	BT->JA = NULL;
	BT->val = NULL; */

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count, curnode[2];

	int num_qp21 = getNumQuadPoints((elementDOF[0].dop - 1) * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp23 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
	double gauss23[num_qp23][3];
	init_Gauss(num_qp23, 2, gauss23); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[i + 1 + elementDOF[0].dof] = count * 2;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		row21[0] = (A->IA[i + 1] - A->IA[i]) / 2;
		rowstart[1] = A->IA[i + elementDOF[0].dof];
		for (j = A->IA[i]; j<A->IA[i] + row21[0]; j++)
		{
			A->JA[j] = istart;
			A->JA[j + row21[0]] = istart + elementDOF[0].dof;
			A->JA[j - A->IA[i] + rowstart[1]] = A->JA[j];
			A->JA[j - A->IA[i] + rowstart[1] + row21[0]] = A->JA[j + row21[0]];
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			curnode[0] = elementDOF[0].val[k][k1];
			curnode[1] = curnode[0] + elementDOF[0].dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + (1 - nu)*phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
				  // a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF[0].dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] * nu;
						}
						break;
					}
				} // j1
				  // a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[1] * phi2[0] * nu;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF[0].dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * ((1 - nu)*phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k


	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i < elementDOF[1].dof; i++)
		B->IA[i + 1] += elementDOF[0].col * 2;

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		rowstart[0] = B->IA[elementDOF[1].val[k][0]];
		for (j = 0; j < elementDOF[0].col; j++)
		{
			B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			B->JA[rowstart[0] + j + elementDOF[0].col] = elementDOF[0].val[k][j] + elementDOF[0].dof;
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		i = elementDOF[1].val[k][0];
		row21[0] = elementDOF[0].col;
		for (k2 = 0; k2<elementDOF[0].col; k2++)
		{
			j = elementDOF[0].val[k][k2];
			// b11
			for (j1 = B->IA[i]; j1<B->IA[i] + row21[0]; j1++)
			{
				if (B->JA[j1] == j)
				{
					for (i1 = 0; i1<num_qp22; i1++)
					{
						lambdas[0] = gauss22[i1][0];
						lambdas[1] = gauss22[i1][1];
						lambdas[2] = 1 - lambdas[0] - lambdas[1];
						ncp1_basis1(lambdas, s, eta, xi, k2, phi1);
						B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[1];
						B->val[j1 + row21[0]] += 2 * s*gauss22[i1][2] * phi1[0];
					}
					break;
				}
				if (j1 == (B->IA[i] + row21[0]))
				{
					printf("There is something wrong in constructing b11 of B\n");
					exit(1);
				}
			}
		} // k1
	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				ncp1_basis(lambdas, k1, &phi);
				for (j = 0; j<elementDOFpm[0].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, j, phi1);
					j1 = elementDOFpm[0].val[k][j];
					temp[0] += phi1[0] * uhpm[0].val[j1];
					temp[1] += phi1[1] * uhpm[0].val[j1];
				}
				for (j = 0; j<elementDOFpm[1].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + j, phi1);
					j1 = elementDOFpm[1].val[k][j];
					temp[0] += phi1[0] * uhpm[1].val[j1];
					temp[1] += phi1[1] * uhpm[1].val[j1];
				}
				b1->val[i] += 2 * s*gauss0[i1][2] * temp[0] * phi;
				b1->val[i + elementDOF[0].dof] += 2 * s*gauss0[i1][2] * temp[1] * phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemble_kirchhoffplate_ElasNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_kirchhoffplate_ElasNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, dvector *uhpm, ELEMENT_DOF *elementDOFpm, int problem_num, double nu)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	if (problem_num == 1)
		assembleStiffmatrix_biharmonic_ElasNcP1P02d(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, uhpm, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, elementDOFpm);
	else
	{
		assembleStiffmatrix_kirchhoffplate_ElasNcP1P02d(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, uhpm, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, elementDOFpm, nu);
		// sigma / (1-nu) as the new stress for the consideration of solver
		int i;
		for (i = 0; i < ptr_A[0].nnz; i++)
			ptr_A[0].val[i] *= 1 - nu;
		for (i = 0; i < b.row; i++)
			b.val[i] /= 1 - nu;
	}

	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);

	//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[1], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[1]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrix_biharmonic_ElasNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_biharmonic_ElasNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 4;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof * 2;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart, colstart;
	int count;

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
									  //	int *index;
	int istart;
	//	index = (int*)calloc(A->col, sizeof(int));
	//	for (i = 0; i<A->col; i++)
	//		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (k = 0; k < elements->row; k++)
	{
		A->IA[4 * k + 1] = 2;
		A->IA[4 * k + 2] = 2;
		A->IA[4 * k + 3] = 1;
		A->IA[4 * k + 4] = 1;
	}

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		j = A->IA[4 * k];
		A->JA[j] = 4 * k;
		A->JA[j + 1] = 4 * k + 1;

		j = A->IA[4 * k + 1];
		A->JA[j] = 4 * k;
		A->JA[j + 1] = 4 * k + 1;

		j = A->IA[4 * k + 2];
		A->JA[j] = 4 * k + 2;

		j = A->IA[4 * k + 3];
		A->JA[j] = 4 * k + 3;
	}
	//	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		// end set parameters

		j = A->IA[4 * k];
		A->val[j] = s / 2;
		A->val[j + 1] = -s / 2;

		j = A->IA[4 * k + 1];
		A->val[j] = -s / 2;
		A->val[j + 1] = s / 2;

		j = A->IA[4 * k + 2];
		A->val[j] = s;

		j = A->IA[4 * k + 3];
		A->val[j] = s;
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i])*elementDOF[0].col * 2;
		B->IA[i + elementDOF[1].dof + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i])*elementDOF[0].col * 2;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			j1 = (j - elementdofTran->IA[i])*elementDOF[0].col * 2;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[1].dof] + j1;
			colstart = element[0] * elementDOF[0].col * 4;
			for (i1 = 0; i1 < elementDOF[0].col; i1++)
			{
				B->JA[rowstart[0] + i1] = colstart + i1;
				B->JA[rowstart[0] + elementDOF[0].col + i1] = colstart + elementDOF[0].col * 2 + i1;
				B->JA[rowstart[1] + i1] = colstart + elementDOF[0].col + i1;
				B->JA[rowstart[1] + elementDOF[0].col + i1] = colstart + elementDOF[0].col * 3 + i1;
			}
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				//				j = elementDOF[0].val[k][k2];
				// b11
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == k * elementDOF[0].col * 4 + k2)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[0];
							B->val[j1 + elementDOF[0].col] -= 2 * s*gauss22[i1][2] * phi1[1];
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
				// b21
				for (j1 = B->IA[i + elementDOF[1].dof]; j1 < B->IA[i + elementDOF[1].dof + 1]; j1++)
				{
					if (B->JA[j1] == k * elementDOF[0].col * 4 + elementDOF[0].col + k2)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[1];
							B->val[j1 + elementDOF[0].col] -= 2 * s*gauss22[i1][2] * phi1[0];
						}
						break;
					}
					if (j1 == B->IA[i + elementDOF[1].dof + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				ncp1_basis(lambdas, k1, &phi);
				for (j = 0; j<elementDOFpm[0].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, j, phi1);
					j1 = elementDOFpm[0].val[k][j];
					temp[0] += phi1[1] * uhpm[0].val[j1];
					temp[1] -= phi1[0] * uhpm[0].val[j1];
				}
				for (j = 0; j<elementDOFpm[1].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + j, phi1);
					j1 = elementDOFpm[1].val[k][j];
					temp[0] += phi1[1] * uhpm[1].val[j1];
					temp[1] -= phi1[0] * uhpm[1].val[j1];
				}
				b2->val[i] -= 2 * s*gauss0[i1][2] * temp[0] * phi;
				b2->val[i + elementDOF[1].dof] -= 2 * s*gauss0[i1][2] * temp[1] * phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStiffmatrix_kirchhoffplate_ElasNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_kirchhoffplate_ElasNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, dvector *uhpm, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, ELEMENT_DOF *elementDOFpm, double nu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 4;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof * 2;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart, colstart;
	int count;

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
									  //	int *index;
	int istart;
	//	index = (int*)calloc(A->col, sizeof(int));
	//	for (i = 0; i<A->col; i++)
	//		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (k = 0; k < elements->row; k++)
	{
		A->IA[4 * k + 1] = 2;
		A->IA[4 * k + 2] = 2;
		A->IA[4 * k + 3] = 2;
		A->IA[4 * k + 4] = 2;
	}

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		j = A->IA[4 * k];
		A->JA[j] = 4 * k;
		A->JA[j + 1] = 4 * k + 1;

		j = A->IA[4 * k + 1];
		A->JA[j] = 4 * k;
		A->JA[j + 1] = 4 * k + 1;

		j = A->IA[4 * k + 2];
		A->JA[j] = 4 * k + 2;
		A->JA[j + 1] = 4 * k + 3;

		j = A->IA[4 * k + 3];
		A->JA[j] = 4 * k + 2;
		A->JA[j+1] = 4 * k + 3;
	}
	//	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		// end set parameters

		j = A->IA[4 * k];
		A->val[j] = s / 2.0 / (1 - nu);
		A->val[j + 1] = -s / 2.0 / (1 - nu);

		j = A->IA[4 * k + 1];
		A->val[j] = -s / 2.0 / (1 - nu);
		A->val[j + 1] = s / 2.0 / (1 - nu);

		j = A->IA[4 * k + 2];
		A->val[j] = s / (1 - nu*nu);
		A->val[j + 1] = s*nu / (1 - nu*nu);

		j = A->IA[4 * k + 3];
		A->val[j] = s*nu / (1 - nu*nu);
		A->val[j + 1] = s / (1 - nu*nu);
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i])*elementDOF[0].col * 2;
		B->IA[i + elementDOF[1].dof + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i])*elementDOF[0].col * 2;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			j1 = (j - elementdofTran->IA[i])*elementDOF[0].col * 2;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[1].dof] + j1;
			colstart = element[0] * elementDOF[0].col * 4;
			for (i1 = 0; i1 < elementDOF[0].col; i1++)
			{
				B->JA[rowstart[0] + i1] = colstart + i1;
				B->JA[rowstart[0] + elementDOF[0].col + i1] = colstart + elementDOF[0].col * 2 + i1;
				B->JA[rowstart[1] + i1] = colstart + elementDOF[0].col + i1;
				B->JA[rowstart[1] + elementDOF[0].col + i1] = colstart + elementDOF[0].col * 3 + i1;
			}
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				//				j = elementDOF[0].val[k][k2];
				// b11
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == k * elementDOF[0].col * 4 + k2)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[0];
							B->val[j1 + elementDOF[0].col] -= 2 * s*gauss22[i1][2] * phi1[1];
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
				// b21
				for (j1 = B->IA[i + elementDOF[1].dof]; j1 < B->IA[i + elementDOF[1].dof + 1]; j1++)
				{
					if (B->JA[j1] == k * elementDOF[0].col * 4 + elementDOF[0].col + k2)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[1];
							B->val[j1 + elementDOF[0].col] -= 2 * s*gauss22[i1][2] * phi1[0];
						}
						break;
					}
					if (j1 == B->IA[i + elementDOF[1].dof + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				ncp1_basis(lambdas, k1, &phi);
				for (j = 0; j<elementDOFpm[0].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, j, phi1);
					j1 = elementDOFpm[0].val[k][j];
					temp[0] += phi1[1] * uhpm[0].val[j1];
					temp[1] -= phi1[0] * uhpm[0].val[j1];
				}
				for (j = 0; j<elementDOFpm[1].col; j++)
				{
					morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + j, phi1);
					j1 = elementDOFpm[1].val[k][j];
					temp[0] += phi1[1] * uhpm[1].val[j1];
					temp[1] -= phi1[0] * uhpm[1].val[j1];
				}
				b2->val[i] -= 2 * s*gauss0[i1][2] * temp[0] * phi;
				b2->val[i + elementDOF[1].dof] -= 2 * s*gauss0[i1][2] * temp[1] * phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleRHScurl_biharmonic_PoissonMorley2d(dvector *ptr_b, dvector *uhv, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFv)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleRHScurl_biharmonic_PoissonMorley2d(dvector *ptr_b, dvector *uhv, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFv)
{
	int i, j, k, l;
	dvector b[2];

	create_dvector(elementDOF[0].dof, &b[0]);
	create_dvector(elementDOF[1].dof, &b[1]);

//	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count;

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				morley_basis1(lambdas, s, elen, eta, xi, nv, nve, k1, phi1);
				for (j = 0; j<elementDOFv[0].col; j++)
				{
					ncp1_basis(lambdas, j, &phi);
					j1 = elementDOFv[0].val[k][j];
					temp[0] += phi * uhv[0].val[j1];
					temp[1] += phi * uhv[0].val[j1 + elementDOFv[0].dof];
				}
				b[0].val[i] += 2 * s*gauss0[i1][2] * (temp[0] * phi1[1] - temp[1] * phi1[0]);
			} // i1
		} // k1 
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + k1, phi1);
				for (j = 0; j<elementDOFv[0].col; j++)
				{
					ncp1_basis(lambdas, j, &phi);
					j1 = elementDOFv[0].val[k][j];
					temp[0] += phi * uhv[0].val[j1];
					temp[1] += phi * uhv[0].val[j1 + elementDOFv[0].dof];
				}
				b[1].val[i] += 2 * s*gauss0[i1][2] * (temp[0] * phi1[1] - temp[1] * phi1[0]);
			} // i1
		} // k1 
	} // k

	ivector *freenodes = &elementDOF[0].freenodes;
	create_dvector(freenodes->row, &ptr_b[0]);
	for (i1 = 0; i1 < ptr_b[0].row; i1++)
	{
		i = freenodes->val[i1];
		ptr_b[0].val[i1] = b[0].val[i];
	}
	free_dvector(&b[0]);
	freenodes = &elementDOF[1].freenodes;
	create_dvector(freenodes->row, &ptr_b[1]);
	for (i1 = 0; i1 < ptr_b[1].row; i1++)
	{
		i = freenodes->val[i1];
		ptr_b[1].val[i1] = b[1].val[i];
	}
	free_dvector(&b[1]);
}

/**
* \fn void assembleRHSgrad_biharmonic_PoissonMorley2d(dvector *ptr_b, dvector *uhv, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFv)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleRHSgrad_biharmonic_PoissonMorley2d(dvector *ptr_b, dvector *uhv, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, ELEMENT_DOF *elementDOFv)
{
	int i, j, k, l;
	dvector b[2];

	create_dvector(elementDOF[0].dof, &b[0]);
	create_dvector(elementDOF[1].dof, &b[1]);

	//	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3], temp[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	int rowstart[2], row21[2], taustart;
	int count;

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = getNumQuadPoints(2 + 1 - 1, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		nv = elements->nvector[k];
		elen = elements->edgeslength[k];
		for (i = 0; i < elementEdge->col; i++)
			nve[i] = edges->nvector[elementEdge->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				morley_basis1(lambdas, s, elen, eta, xi, nv, nve, k1, phi1);
				for (j = 0; j<elementDOFv[0].col; j++)
				{
					ncp1_basis(lambdas, j, &phi);
					j1 = elementDOFv[0].val[k][j];
					temp[0] += phi * uhv[0].val[j1];
					temp[1] += phi * uhv[0].val[j1 + elementDOFv[0].dof];
				}
				b[0].val[i] += 2 * s*gauss0[i1][2] * (temp[0] * phi1[0] + temp[1] * phi1[1]);
			} // i1
		} // k1 
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				temp[0] = 0;
				temp[1] = 0;
				morley_basis1(lambdas, s, elen, eta, xi, nv, nve, 3 + k1, phi1);
				for (j = 0; j<elementDOFv[0].col; j++)
				{
					ncp1_basis(lambdas, j, &phi);
					j1 = elementDOFv[0].val[k][j];
					temp[0] += phi * uhv[0].val[j1];
					temp[1] += phi * uhv[0].val[j1 + elementDOFv[0].dof];
				}
				b[1].val[i] += 2 * s*gauss0[i1][2] * (temp[0] * phi1[0] + temp[1] * phi1[1]);
			} // i1
		} // k1 
	} // k

	ivector *freenodes = &elementDOF[0].freenodes;
	create_dvector(freenodes->row, &ptr_b[0]);
	for (i1 = 0; i1 < ptr_b[0].row; i1++)
	{
		i = freenodes->val[i1];
		ptr_b[0].val[i1] = b[0].val[i];
	}
	free_dvector(&b[0]);
	freenodes = &elementDOF[1].freenodes;
	create_dvector(freenodes->row, &ptr_b[1]);
	for (i1 = 0; i1 < ptr_b[1].row; i1++)
	{
		i = freenodes->val[i1];
		ptr_b[1].val[i1] = b[1].val[i];
	}
	free_dvector(&b[1]);
}

/**
* \fn void assembleStressMassmatrixKirchhoffplateMorley2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu, double nu)
* \brief assemble stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStressMassmatrixKirchhoffplateMorley2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu, double nu)
{
	int i, j, k;
	double s;
	M->row = elementDOF[0].dof * 4;
	M->col = M->row;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;

	M->IA[0] = 0;
	for (k = 0; k < elements->row; k++)
	{
		M->IA[4 * k + 1] = 1;
		M->IA[4 * k + 2] = 1;
		M->IA[4 * k + 3] = 2;
		M->IA[4 * k + 4] = 2;
	}
	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];
	M->nnz = M->IA[M->row];

	M->JA = (int*)malloc(M->nnz*sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		j = M->IA[4 * k];
		M->JA[j] = 4 * k;

		j = M->IA[4 * k + 1];
		M->JA[j] = 4 * k + 1;

		j = M->IA[4 * k + 2];
		M->JA[j] = 4 * k + 2;
		M->JA[j + 1] = 4 * k + 3;

		j = M->IA[4 * k + 3];
		M->JA[j] = 4 * k + 2;
		M->JA[j + 1] = 4 * k + 3;
	}

	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		j = M->IA[4 * k];
		M->val[j] = s / (2 * mu);
		
		j = M->IA[4 * k + 1];
		M->val[j] = s / (2 * mu);
		
		j = M->IA[4 * k + 2];
		M->val[j] = s / (2 * mu) / (1 + nu);
		M->val[j + 1] = s*nu / (2 * mu) / (1 + nu);

		j = M->IA[4 * k + 3];
		M->val[j] = s*nu / (2 * mu) / (1 + nu);
		M->val[j + 1] = s / (2 * mu) / (1 + nu);
	} // k
}

/**
* \fn void assembleStressMassmatrixKirchhoffplateMorleydBD2d(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu, double nu)
* \brief assemble stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStressMassmatrixKirchhoffplateMorleydBD2d(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu, double nu)
{
	int i, j, k;
	create_dbd_matrix(elementDOF->dof * 4, elementDOF->dof * 4, elements->row * 3, M);
	for (i = 0; i < elements->row; i++)
	{
		create_dden_matrix(1, 1, M->blk + 3 * i);
		create_dden_matrix(1, 1, M->blk + 3 * i + 1);
		create_dden_matrix(2, 2, M->blk + 3 * i + 2);
	}

	double s;

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters
		M->blk[3 * k].val[0][0] = s / (2 * mu);
		M->blk[3 * k + 1].val[0][0] = s / (2 * mu);
		M->blk[3 * k + 2].val[0][0] = s / (2 * mu) / (1 + nu);
		M->blk[3 * k + 2].val[0][1] = s*nu / (2 * mu) / (1 + nu);
		M->blk[3 * k + 2].val[1][0] = s*nu / (2 * mu) / (1 + nu);
		M->blk[3 * k + 2].val[1][1] = s / (2 * mu) / (1 + nu);
	} // k
}

/**
* \fn void assemble_stokesCR(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemble_stokesCR(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = b
	Bx1 + Cx2   = 0
	where x1: u_h, x2: p_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrixStokesCR2d(&A, &B, &ptr_A[3], &b, &ptr_b[1], elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);
	
	// initial solution
	create_dvector(b.row, &uh[0]);
	create_dvector(ptr_b[1].row, &uh[1]);
	init_dvector(&uh[0], 1);///////////////////////////////////

	//extract
	extractFreenodesVector2StokesDirichlet(&A, &B, &b, ptr_b, &elementDOF[0], &uh[0]);
	free_dvector(&b);
	extractFreenodesMatrix11(&A, &ptr_A[0], &elementDOF[0], &elementDOF[0]);
	free_csr_matrix(&A);
	extractFreenodesMatrix1c(&B, &ptr_A[2], &elementDOF[0]);
	free_csr_matrix(&B);

	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrixStokesCR2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrixStokesCR2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node, curnode[2];

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints((3 - 1) * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(3 + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  //	int num_qp23 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
									  //	double gauss23[num_qp23][3];
									  //	init_Gauss(num_qp23, 2, gauss23); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[i + 1 + elementDOF[0].dof] = count * 2;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		row21[0] = (A->IA[i + 1] - A->IA[i]) / 2;
		rowstart[1] = A->IA[i + elementDOF[0].dof];
		for (j = A->IA[i]; j<A->IA[i] + row21[0]; j++)
		{
			A->JA[j] = istart;
			A->JA[j + row21[0]] = istart + elementDOF[0].dof;
			A->JA[j - A->IA[i] + rowstart[1]] = A->JA[j];
			A->JA[j - A->IA[i] + rowstart[1] + row21[0]] = A->JA[j + row21[0]];
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			curnode[0] = elementDOF[0].val[k][k1];
			curnode[1] = curnode[0] + elementDOF[0].dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							cr_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2) * 2 * mu;
						}
						break;
					}
				} // j1
				  // a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF[0].dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							cr_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[1] * phi2[0] / 2 * 2 * mu;
						}
						break;
					}
				} // j1
				  // a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							cr_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] / 2 * 2 * mu;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF[0].dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							cr_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2		
		} // k1
	} // k


	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i < elementDOF[1].dof; i++)
		B->IA[i + 1] += elementDOF[0].col * 2;

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = B->IA[elementDOF[1].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
			{
				B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
				B->JA[rowstart[0] + j + elementDOF[0].col] = elementDOF[0].val[k][j] + elementDOF[0].dof;
			}
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			row21[0] = (B->IA[i + 1] - B->IA[i]) / 2;
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// b11
				for (j1 = B->IA[i]; j1 < B->IA[i] + row21[0]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, &phi);
							cr_basis1(lambdas, s, eta, xi, k2, phi1);
							B->val[j1] += 2 * s*gauss22[i1][2] * phi * phi1[0];
							B->val[j1 + row21[0]] += 2 * s*gauss22[i1][2] * phi * phi1[1];
						}
						break;
					}
					if (j1 == (B->IA[i] + row21[0]))
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k


	  /************************************************** stiffness matrix C *****************************************************************/
	if (t > 0)
	{
		assemblePressMassmatrixStokesCR2d(C, elements, &elementDOF[1], -t);
	}
	else
	{
		C->row = 0;
		C->col = 0;
		C->IA = NULL;
		C->JA = NULL;
		C->val = NULL;
	}

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				cr_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b1->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b1->val[i + elementDOF[0].dof] += 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemblePressMassmatrixStokesCR2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double t)
* \brief assemble mass matrix
* \param *M pointer to mass matrix
* \param *elements pointer to the structure of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressMassmatrixStokesCR2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double t)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int count, rowstart[2];

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

									  /************************************************** stiffness matrix M *****************************************************************/
	M->row = elementDOF->dof;
	M->col = elementDOF->dof;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;
	// step 1M: Find first the structure IA of the stiffness matrix M
	M->IA[0] = 0;
	for (i = 0; i<elementDOF->dof; i++)
		M->IA[i + 1] = elementDOF->col;

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];
	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)calloc(M->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF->col; i++)
		{
			rowstart[0] = M->IA[elementDOF->val[k][i]];
			for (j = 0; j < elementDOF->col; j++)
				M->JA[rowstart[0] + j] = elementDOF->val[k][j];
		}
	}

	// step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//			xi = elements->xi[k];
		//			eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			for (k2 = 0; k2 < elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * t;
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemblePressMassmatrixStokesCRdBD2d(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double t)
* \brief assemble pressure-stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressMassmatrixStokesCRdBD2d(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double t)
{
	int i, j, k, l;
	int dop = elementDOF->dop;
	int col = (dop + 1)*(dop + 2) / 2;
	create_dbd_matrix(elementDOF->dof, elementDOF->dof, elements->row, M);
	for (i = 0; i<elements->row; i++)
		create_dden_matrix(col, col, M->blk + i);

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;

	int num_qp21 = getNumQuadPoints(dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<col; k1++)
		{
			for (k2 = 0; k2<col; k2++)
			{
				for (i1 = 0; i1<num_qp21; i1++)
				{
					lambdas[0] = gauss21[i1][0];
					lambdas[1] = gauss21[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, dop, &phi1[0]);
					lagrange_basis(lambdas, k2, dop, &phi2[0]);
					M->blk[k].val[k1][k2] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * t;
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesMINI(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemble_stokesMINI(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = b
	Bx1 + Cx2   = 0
	where x1: u_h, x2: p_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrixStokesMINI2d(&A, &B, &ptr_A[3], &b, &ptr_b[1], elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(b.row, &uh[0]);
	create_dvector(ptr_b[1].row, &uh[1]);
	init_dvector(&uh[0], 1);///////////////////////////////////

							//extract
	extractFreenodesVector2StokesDirichlet(&A, &B, &b, ptr_b, &elementDOF[0], &uh[0]);
	free_dvector(&b);
	extractFreenodesMatrix11(&A, &ptr_A[0], &elementDOF[0], &elementDOF[0]);
	free_csr_matrix(&A);
	extractFreenodesMatrix1c(&B, &ptr_A[2], &elementDOF[0]);
	free_csr_matrix(&B);

	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrixStokesMINI2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrixStokesMINI2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node, curnode[2];

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints((3 - 1) * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(3 + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  //	int num_qp23 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
									  //	double gauss23[num_qp23][3];
									  //	init_Gauss(num_qp23, 2, gauss23); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[i + 1 + elementDOF[0].dof] = count * 2;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		row21[0] = (A->IA[i + 1] - A->IA[i]) / 2;
		rowstart[1] = A->IA[i + elementDOF[0].dof];
		for (j = A->IA[i]; j<A->IA[i] + row21[0]; j++)
		{
			A->JA[j] = istart;
			A->JA[j + row21[0]] = istart + elementDOF[0].dof;
			A->JA[j - A->IA[i] + rowstart[1]] = A->JA[j];
			A->JA[j - A->IA[i] + rowstart[1] + row21[0]] = A->JA[j + row21[0]];
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	  //	free(index);

	  // step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			curnode[0] = elementDOF[0].val[k][k1];
			curnode[1] = curnode[0] + elementDOF[0].dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							mini_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2) * 2 * mu;
						}
						break;
					}
				} // j1
				  // a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF[0].dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							mini_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[1] * phi2[0] / 2 * 2 * mu;
						}
						break;
					}
				} // j1
				  // a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							mini_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] / 2 * 2 * mu;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF[0].dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							mini_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2		
		} // k1
	} // k


	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran[1].IA[i]; j<elementdofTran[1].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[1].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		B->IA[i + 1] = count * 2;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran[1].IA[i]; j<elementdofTran[1].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[1].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		row21[0] = (B->IA[i + 1] - B->IA[i]) / 2;
		for (j = B->IA[i]; j<B->IA[i] + row21[0]; j++)
		{
			B->JA[j] = istart;
			B->JA[j + row21[0]] = istart + elementDOF[0].dof;
			istart = index[istart];
			index[B->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			row21[0] = (B->IA[i + 1] - B->IA[i]) / 2;
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// b11
				for (j1 = B->IA[i]; j1 < B->IA[i] + row21[0]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, &phi);
							mini_basis1(lambdas, s, eta, xi, k2, phi1);
							B->val[j1] += 2 * s*gauss22[i1][2] * phi * phi1[0];
							B->val[j1 + row21[0]] += 2 * s*gauss22[i1][2] * phi * phi1[1];
						}
						break;
					}
					if (j1 == (B->IA[i] + row21[0]))
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k


	  /************************************************** stiffness matrix C *****************************************************************/
	if (t > 0)
	{
		assemblePressMassmatrixStokesMINI2d(C, elements, &elementDOF[1], &elementdofTran[1], -t);
	}
	else
	{
		C->row = 0;
		C->col = 0;
		C->IA = NULL;
		C->JA = NULL;
		C->val = NULL;
	}

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				mini_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b1->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b1->val[i + elementDOF[0].dof] += 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemblePressMassmatrixStokesMINI2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double t)
* \brief assemble mass matrix
* \param *M pointer to mass matrix
* \param *elements pointer to the structure of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressMassmatrixStokesMINI2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double t)
{
	int i, j, k, l;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

									  /************************************************** stiffness matrix M *****************************************************************/
	M->row = elementDOF->dof;
	M->col = elementDOF->dof;
	M->IA = (int*)calloc(M->row + 1, sizeof(int));
	M->JA = NULL;
	M->val = NULL;
	int *index;
	int istart;
	index = (int*)malloc(M->col * sizeof(int));
	for (i = 0; i<M->col; i++)
		index[i] = -1;
	// step 1M: Find first the structure IA of the stiffness matrix M
	for (i = 0; i<elementDOF->dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		M->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];
	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)calloc(M->nnz, sizeof(int));
	for (i = 0; i<elementDOF->dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = M->IA[i]; j<M->IA[i + 1]; j++)
		{
			M->JA[j] = istart;
			istart = index[istart];
			index[M->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//			xi = elements->xi[k];
		//			eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			for (k2 = 0; k2 < elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * t;
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemble_stokesNcP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = b
	Bx1 + Cx2   = 0
	where x1: u_h, x2: p_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrixStokesNcP1P02d(&A, &B, &ptr_A[3], &b, &ptr_b[1], elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(b.row, &uh[0]);
	create_dvector(ptr_b[1].row, &uh[1]);
	init_dvector(&uh[0], 1);///////////////////////////////////

							//extract
	extractFreenodesVector2StokesDirichlet(&A, &B, &b, ptr_b, &elementDOF[0], &uh[0]);
	free_dvector(&b);
	extractFreenodesMatrix11(&A, &ptr_A[0], &elementDOF[0], &elementDOF[0]);
	free_csr_matrix(&A);
	extractFreenodesMatrix1c(&B, &ptr_A[2], &elementDOF[0]);
	free_csr_matrix(&B);

	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrixStokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrixStokesNcP1P02d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	/*BT->row = A->row;
	BT->col = elementDOF[1].dof;
	BT->IA = (int*)calloc(BT->row + 1, sizeof(int));
	BT->JA = NULL;
	BT->val = NULL; */

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;

	int num_qp21 = getNumQuadPoints((elementDOF[0].dop - 1) * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp23 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
	double gauss23[num_qp23][3];
	init_Gauss(num_qp23, 2, gauss23); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;
		A->IA[i + 1 + elementDOF[0].dof] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
			A->JA[j - A->IA[i] + A->IA[i + elementDOF[0].dof]] = istart + elementDOF[0].dof;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							ncp1_basis1(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	for (i = 0; i < A->IA[elementDOF[0].dof]; i++)
		A->val[i + A->IA[elementDOF[0].dof]] = A->val[i];

	/************************************************** stiffness matrix B *****************************************************************/
	// step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i < elementDOF[1].dof; i++)
		B->IA[i + 1] += elementDOF[0].col * 2;

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		rowstart[0] = B->IA[elementDOF[1].val[k][0]];
		for (j = 0; j < elementDOF[0].col; j++)
		{
			B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			B->JA[rowstart[0] + j + elementDOF[0].col] = elementDOF[0].val[k][j] + elementDOF[0].dof;
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		i = elementDOF[1].val[k][0];
		row21[0] = elementDOF[0].col;
		for (k2 = 0; k2<elementDOF[0].col; k2++)
		{
			j = elementDOF[0].val[k][k2];
			// b11
			for (j1 = B->IA[i]; j1<B->IA[i] + row21[0]; j1++)
			{
				if (B->JA[j1] == j)
				{
					for (i1 = 0; i1<num_qp22; i1++)
					{
						lambdas[0] = gauss22[i1][0];
						lambdas[1] = gauss22[i1][1];
						lambdas[2] = 1 - lambdas[0] - lambdas[1];
						ncp1_basis1(lambdas, s, eta, xi, k2, phi1);
						B->val[j1] += 2 * s*gauss22[i1][2] * phi1[0];
						B->val[j1 + row21[0]] += 2 * s*gauss22[i1][2] * phi1[1];
					}
					break;
				}
				if (j1 == (B->IA[i] + row21[0]))
				{
					printf("There is something wrong in constructing b11 of B\n");
					exit(1);
				}
			}
		} // k1
	} // k

	  /************************************************** stiffness matrix BT *****************************************************************/
	  /*
	  // step 1BT: Find first the structure IA of the stiffness matrix BT
	  for (i = 0; i<elementDOF[0].dof; i++)
	  {
	  BT->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]);
	  BT->IA[i + 1 + elementDOF[0].dof] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]);

	  } // i

	  for (i = 0; i<BT->row; i++)
	  BT->IA[i + 1] += BT->IA[i];

	  BT->nnz = BT->IA[BT->row];

	  // step 2BT: Find the structure JA of the stiffness matrix BT
	  BT->JA = (int*)calloc(BT->nnz, sizeof(int));
	  for (i = 0; i<elementDOF[0].dof; i++)
	  {
	  for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
	  {
	  element[0] = elementdofTran->JA[j];
	  BT->JA[BT->IA[i] + j - elementdofTran->IA[i]] = element[0];
	  BT->JA[BT->IA[i + elementDOF[0].dof] + j - elementdofTran->IA[i]] = element[0];
	  } // j
	  } // i

	  // step 3BT: Loop element by element and compute the actual entries storing them in BT
	  BT->val = (double*)calloc(BT->nnz, sizeof(double));
	  for (k = 0; k<elements->row; k++)
	  {
	  // set parameters
	  s = elements->vol[k];
	  xi = elements->xi[k];
	  eta = elements->eta[k];
	  // end set parameters

	  for (k1 = 0; k1<elementDOF[0].col; k1++)
	  {
	  i = elementDOF[0].val[k][k1];
	  j = elementDOF[1].val[k][0];
	  // b11
	  for (j1 = BT->IA[i]; j1<BT->IA[i + 1]; j1++)
	  {
	  if (BT->JA[j1] == j)
	  {
	  for (i1 = 0; i1<num_qp22; i1++)
	  {
	  lambdas[0] = gauss22[i1][0];
	  lambdas[1] = gauss22[i1][1];
	  lambdas[2] = 1 - lambdas[0] - lambdas[1];
	  ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
	  BT->val[j1] += 2 * s*gauss22[i1][2] * phi1[0];
	  }
	  break;
	  }
	  if (j1 == (BT->IA[i + 1]))
	  {
	  printf("There is something wrong in constructing b11 of BT\n");
	  exit(1);
	  }
	  }
	  // b21
	  for (j1 = BT->IA[i + elementDOF[0].dof]; j1<BT->IA[i + elementDOF[0].dof + 1]; j1++)
	  {
	  if (BT->JA[j1] == j)
	  {
	  for (i1 = 0; i1<num_qp22; i1++)
	  {
	  lambdas[0] = gauss22[i1][0];
	  lambdas[1] = gauss22[i1][1];
	  lambdas[2] = 1 - lambdas[0] - lambdas[1];
	  ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
	  BT->val[j1] += 2 * s*gauss22[i1][2] * phi1[1];
	  }
	  break;
	  }
	  if (j1 == (BT->IA[i + elementDOF[0].dof + 1]))
	  {
	  printf("There is something wrong in constructing b21 of BT\n");
	  exit(1);
	  }
	  }
	  } // k1
	  } // k
	  */

	  /************************************************** stiffness matrix C *****************************************************************/
	if (t > 0)
	{
		C->row = elementDOF[1].dof;
		C->col = elementDOF[1].dof;
		C->IA = (int*)calloc(C->row + 1, sizeof(int));
		C->JA = NULL;
		C->val = NULL;
		for (i = 0; i<elementDOF[1].dof; i++)
			C->IA[i + 1] = 1;
		for (i = 0; i<C->row; i++)
			C->IA[i + 1] += C->IA[i];
		C->nnz = C->IA[C->row];
		C->JA = (int*)calloc(C->nnz, sizeof(int));
		for (i = 0; i<C->nnz; i++)
			C->JA[i] = i;
		C->val = (double*)calloc(C->nnz, sizeof(double));
		for (k = 0; k < elements->row; k++)
			C->val[k] = -elements->vol[k] * t;
	}
	else
	{
		C->row = 0;
		C->col = 0;
		C->IA = NULL;
		C->JA = NULL;
		C->val = NULL;
	}

	/************************************************** right hand side b1 *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				ncp1_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b1->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b1->val[i + elementDOF[0].dof] += 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemblePressMassmatrixStokesNcP1P02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF)
* \brief assemble mass matrix
* \param *M pointer to mass matrix
* \param *elements pointer to the structure of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assemblePressMassmatrixStokesNcP1P02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF)
{
	int i, j, k;

	M->row = elementDOF->dof;
	M->col = elementDOF->dof;
	M->IA = (int*)calloc(M->row + 1, sizeof(int));
	M->JA = NULL;
	M->val = NULL;
	for (i = 0; i<M->row; i++)
		M->IA[i + 1] = 1;
	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];
	M->nnz = M->IA[M->row];
	M->JA = (int*)calloc(M->nnz, sizeof(int));
	for (i = 0; i<M->nnz; i++)
		M->JA[i] = i;
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k < elements->row; k++)
		M->val[k] = elements->vol[k];
}

/**
* \fn void assemble_stokesCR_elas_stressCompact(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_stokesCR_elas_stressCompact(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrix_stokesCR_elas_stressCompact(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);
	init_dvector(&uh[1], 1);///////////////////////////////////

	//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[2], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[2]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);

	/////////////////////
	/*	int i,j;
	for (i = 0; i < ptr_A[1].row; i++)
	{
	for (j = ptr_A[1].IA[i]; j < ptr_A[1].IA[i + 1]; j++)
	{
	printf("( %d, %d, %e); ", i, ptr_A[1].JA[j], ptr_A[1].val[j]);
	}
	printf("\n");
	}*/
	//	exit(1);
}

/**
* \fn void assembleStiffmatrix_stokesCR_elas_stressCompact(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_stokesCR_elas_stressCompact(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	if (t < 0)
		t = 0;

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	int num_qp20 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss20[num_qp20][3];
	init_Gauss(num_qp20, 2, gauss20); // gauss intergation initial

	int num_qp21 = getNumQuadPoints(3, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	dBDmat Mp, Mpinv;
	dCSRmat Msp, MspT, Ms, tempA, tempB;
	// Mp, Mpinv
	create_dbd_matrix(elementDOF[0].dof, elementDOF[0].dof, elements->row, &Mp);
	for (i = 0; i<Mp.nb; i++)
		create_dden_matrix(elementDOF[0].col, elementDOF[0].col, Mp.blk + i);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				for (i1 = 0; i1<num_qp20; i1++)
				{
					lambdas[0] = gauss20[i1][0];
					lambdas[1] = gauss20[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, elementDOF[0].dop, &phi1[0]);
					lagrange_basis(lambdas, k2, elementDOF[0].dop, &phi2[0]);
					Mp.blk[k].val[k1][k2] += 2 * s*gauss20[i1][2] * phi1[0] * phi2[0] * (1.0 / mu + t);
				}
			} // k2
		} // k1
	} // k
	inverse_dBDmat(&Mp, &Mpinv);
	free_dbd_matrix(&Mp);

	// Msp, MspT
	Msp.row = elementDOF[1].dof;
	Msp.col = elementDOF[0].dof;
	Msp.IA = (int*)malloc((Msp.row + 1) * sizeof(int));
	Msp.JA = NULL;
	Msp.val = NULL;
	// step 1Msp: Find first the structure IA of the mass matrix Msp
	Msp.IA[0] = 0;
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			if (k1 < 6)
				Msp.IA[i + 1] = elementDOF[0].col;
			else if (k1<9)
				Msp.IA[i + 1] = 0;
			else
				Msp.IA[i + 1] = elementDOF[0].col;
		}
	}
	for (i = 0; i<Msp.row; i++)
		Msp.IA[i + 1] += Msp.IA[i];
	Msp.nnz = Msp.IA[Msp.row];
	// step 2Msp: Find the structure JA of the mass matrix Msp
	Msp.JA = (int*)malloc(Msp.nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			if (i >= 6 && i<9)
				continue;

			rowstart[0] = Msp.IA[elementDOF[1].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				Msp.JA[rowstart[0] + j] = elementDOF[0].val[k][j];
		} // i
	} // k
	  // step 3Msp: Loop element by element and compute the actual entries storing them in Msp
	Msp.val = (double*)calloc(Msp.nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters
		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			if (k1 >= 6 && k1<9)
				continue;

			i = elementDOF[1].val[k][k1];
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = Msp.IA[i]; j1 < Msp.IA[i + 1]; j1++)
				{
					if (Msp.JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							//	lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							Msp.val[j1] += 2 * s*gauss21[i1][2] * phi2[0] * (phi1[0] + phi1[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == Msp.IA[i + 1])
					{
						printf("There is something wrong in constructing Msp\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
	getTransposeOfSparse(&Msp, &MspT);
	assembleStressMassmatrixStokesCR2d_stressCompact(&Ms, elements, elementDOF + 1, mu);
	/************ form Shur complement A = Ms - MspMp^{-1}MspT  *********/
	dBDMultiplydCSR(&Mpinv, &MspT, &tempA);
	free_dbd_matrix(&Mpinv);
	free_csr_matrix(&MspT);
	sparseMultiplication(&Msp, &tempA, &tempB);
	free_csr_matrix(&Msp);
	free_csr_matrix(&tempA);
	sparseSubtraction(&Ms, &tempB, A);
	free_csr_matrix(&Ms);
	free_csr_matrix(&tempB);

	/************************************************** stiffness matrix B *****************************************************************/
	B->row = elementDOF[2].dof * 2;
	B->col = A->col;
	B->IA = (int*)malloc((B->row + 1) * sizeof(int));
	B->JA = NULL;
	B->val = NULL;
	// step 1B: Find first the structure IA of the stiffness matrix B
	B->IA[0] = 0;
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 8;
		B->IA[i + elementDOF[2].dof + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 8;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)malloc(B->nnz * sizeof(int));
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			j1 = (j - elementdofTran->IA[i]) * 8;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[2].dof] + j1;
			for (i1 = 0; i1 < 3; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1];
			for (i1 = 3; i1 < 8; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][3 + i1];
			for (i1 = 0; i1 < 8; i1++)
				B->JA[rowstart[1] + i1] = elementDOF[1].val[element[0]][3 + i1];
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		// b11
		count = 8;
		for (i = 0; i < 3; i++)
			patchnodes[i] = i;
		for (i = 3; i < 8; i++)
			patchnodes[i] = 3 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2];
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							//	lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		  // b21
		count = 8;
		for (i = 0; i < 8; i++)
			patchnodes[i] = 3 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1] + elementDOF[2].dof;
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2];
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							//	lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[1] * phi2[1] + phi1[0] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				cr_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] -= 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b2->val[i + elementDOF[2].dof] -= 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStressMassmatrixStokesCR2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStressMassmatrixStokesCR2d_stressCompact(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	M->row = elementDOF->dof;
	M->col = M->row;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial


									  /************************************************** stiffness matrix M *****************************************************************/
									  // step 1M: Find first the structure IA of the stiffness matrix M
	M->IA[0] = 0;
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			if (k1 < 6)
				M->IA[i + 1] = 4;
			else if (k1<9)
				M->IA[i + 1] = 5;
			else
				M->IA[i + 1] = 8;
		}
	}

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];

	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)malloc(M->nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF->col; i++)
		{
			rowstart[0] = M->IA[elementDOF->val[k][i]];
			if (i < 3)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][j];
				M->JA[rowstart[0] + 3] = elementDOF->val[k][9];
			}
			else if (i < 6)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][3 + j];
				M->JA[rowstart[0] + 3] = elementDOF->val[k][10];
			}
			else if (i < 9)
			{
				for (j = 0; j < 5; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][6 + j];
			}
			else if (i == 9)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][j];
				for (j = 3; j < 8; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][3 + j];
			}
			else // (i == 10)
			{
				for (j = 0; j < 8; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][3 + j];
			}
		} // i
	} // k

	  // step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			if (k1 < 3)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				patchnodes[3] = 9;
			}
			else if (k1 < 6)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = 3 + j;
				patchnodes[3] = 10;
			}
			else if (k1 < 9)
			{
				count = 5;
				for (j = 0; j < 5; j++)
					patchnodes[j] = 6 + j;
			}
			else if (k1 == 9)
			{
				count = 8;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				for (j = 3; j < 8; j++)
					patchnodes[j] = 3 + j;
			}
			else // k1==10 
			{
				count = 8;
				for (j = 0; j < 8; j++)
					patchnodes[j] = 3 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF->val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							//		lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							//		lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assembleStressMassmatrixStokesCRdBD2d(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStressMassmatrixStokesCRdBD2d_stressCompact(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;
	//	int dop = elementDOF->dop;
	int col = 11;
	create_dbd_matrix(elementDOF->dof, elementDOF->dof, elements->row, M);
	for (i = 0; i<M->nb; i++)
		create_dden_matrix(col, col, M->blk + i);

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;

	int num_qp21 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<col; k1++)
		{
			for (k2 = 0; k2<col; k2++)
			{
				for (i1 = 0; i1<num_qp21; i1++)
				{
					lambdas[0] = gauss21[i1][0];
					lambdas[1] = gauss21[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
					crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
					M->blk[k].val[k1][k2] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesCR_elas_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_stokesCR_elas_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrix_stokesCR_elas_stressP2(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);
	init_dvector(&uh[1], 1);///////////////////////////////////

	//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[2], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[2]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);

	/////////////////////
/*	int i,j;
	for (i = 0; i < ptr_A[1].row; i++)
	{
		for (j = ptr_A[1].IA[i]; j < ptr_A[1].IA[i + 1]; j++)
		{
			printf("( %d, %d, %e); ", i, ptr_A[1].JA[j], ptr_A[1].val[j]);
		}
		printf("\n");
	}*/
//	exit(1);
}

/**
* \fn void assembleStiffmatrix_stokesCR_elas_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_stokesCR_elas_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	if (t < 0)
		t = 0;

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	int num_qp20 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss20[num_qp20][3];
	init_Gauss(num_qp20, 2, gauss20); // gauss intergation initial

	int num_qp21 = getNumQuadPoints(3, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	dBDmat Mp, Mpinv;
	dCSRmat Msp, MspT, Ms, tempA, tempB;
	// Mp, Mpinv
	create_dbd_matrix(elementDOF[0].dof, elementDOF[0].dof, elements->row, &Mp);
	for (i = 0; i<Mp.nb; i++)
		create_dden_matrix(elementDOF[0].col, elementDOF[0].col, Mp.blk + i);
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				for (i1 = 0; i1<num_qp20; i1++)
				{
					lambdas[0] = gauss20[i1][0];
					lambdas[1] = gauss20[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, elementDOF[0].dop, &phi1[0]);
					lagrange_basis(lambdas, k2, elementDOF[0].dop, &phi2[0]);
					Mp.blk[k].val[k1][k2] += 2 * s*gauss20[i1][2] * phi1[0] * phi2[0] * (1.0 / mu + t);
				}
			} // k2
		} // k1
	} // k
	inverse_dBDmat(&Mp, &Mpinv);
	free_dbd_matrix(&Mp);

	// Msp, MspT
	Msp.row = elementDOF[1].dof;
	Msp.col = elementDOF[0].dof;
	Msp.IA = (int*)malloc((Msp.row + 1) * sizeof(int));
	Msp.JA = NULL;
	Msp.val = NULL;
	// step 1Msp: Find first the structure IA of the mass matrix Msp
	Msp.IA[0] = 0;
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			if (k1 < 6 * 2)
				Msp.IA[i + 1] = elementDOF[0].col;
			else
				Msp.IA[i + 1] = 0;
		}
	}
	for (i = 0; i<Msp.row; i++)
		Msp.IA[i + 1] += Msp.IA[i];
	Msp.nnz = Msp.IA[Msp.row];
	// step 2Msp: Find the structure JA of the mass matrix Msp
	Msp.JA = (int*)malloc(Msp.nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < 12; i++)
		{
			rowstart[0] = Msp.IA[elementDOF[1].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				Msp.JA[rowstart[0] + j] = elementDOF[0].val[k][j];
		} // i
	} // k
	  // step 3Msp: Loop element by element and compute the actual entries storing them in Msp
	Msp.val = (double*)calloc(Msp.nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
//		xi = elements->xi[k];
//		eta = elements->eta[k];
		// end set parameters
		for (k1 = 0; k1 < 12; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = Msp.IA[i]; j1 < Msp.IA[i + 1]; j1++)
				{
					if (Msp.JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							//				crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							Msp.val[j1] += 2 * s*gauss21[i1][2] * phi2[0] * (phi1[0] + phi1[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == Msp.IA[i+1])
					{
						printf("There is something wrong in constructing Msp\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
	getTransposeOfSparse(&Msp, &MspT);
	assembleStressMassmatrixStokesCR2d_stressP2(&Ms, elements, elementDOF+1, mu);
	/************ form Shur complement A = Ms - MspMp^{-1}MspT  *********/
	dBDMultiplydCSR(&Mpinv, &MspT, &tempA);
	free_dbd_matrix(&Mpinv);
	free_csr_matrix(&MspT);
	sparseMultiplication(&Msp, &tempA, &tempB);
	free_csr_matrix(&Msp);
	free_csr_matrix(&tempA);
	sparseSubtraction(&Ms, &tempB, A);
	free_csr_matrix(&Ms);
	free_csr_matrix(&tempB);

	  /************************************************** stiffness matrix B *****************************************************************/
	B->row = elementDOF[2].dof * 2;
	B->col = A->col;
	B->IA = (int*)malloc((B->row + 1) * sizeof(int));
	B->JA = NULL;
	B->val = NULL;
	  // step 1B: Find first the structure IA of the stiffness matrix B
	B->IA[0] = 0;
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 12;
		B->IA[i + elementDOF[2].dof + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 12;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)malloc(B->nnz*sizeof(int));
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			j1 = (j - elementdofTran->IA[i]) * 12;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[2].dof] + j1;
			for (i1 = 0; i1 < 6; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1];
			for (i1 = 6; i1 < 12; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][6 + i1];
			for (i1 = 0; i1 < 12; i1++)
				B->JA[rowstart[1] + i1] = elementDOF[1].val[element[0]][6 + i1];
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		// b11
		count = 12;
		for (i = 0; i < 6; i++)
			patchnodes[i] = i;
		for (i = 6; i < 12; i++)
			patchnodes[i] = 6 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2];
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
				//			crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		// b21
		count = 12;
		for (i = 0; i < 12; i++)
			patchnodes[i] = 6 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1] + elementDOF[2].dof;
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2];
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
					//		crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[1] * phi2[1] + phi1[0] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				cr_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] -= 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b2->val[i + elementDOF[2].dof] -= 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStressMassmatrixStokesCR2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStressMassmatrixStokesCR2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	M->row = elementDOF->dof;
	M->col = M->row;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial


									  /************************************************** stiffness matrix M *****************************************************************/
	// step 1M: Find first the structure IA of the stiffness matrix M
	M->IA[0] = 0;
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			M->IA[i + 1] = 6;
		}
	}

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];

	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)malloc(M->nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF->col; i++)
		{
			rowstart[0] = M->IA[elementDOF->val[k][i]];
			if (i < 6)
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][j];
			}
			else if (i < 12)
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][6 + j];
			}
			else
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF->val[k][12 + j];
			}
		} // i
	} // k

	  // step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
//		xi = elements->xi[k];
//		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			if (k1 < 6)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = j;
			}
			else if (k1 < 12)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 6 + j;
			}
			else
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 12 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF->val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//	crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							//	crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assembleStressMassmatrixStokesCRdBD2d_stressP2(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStressMassmatrixStokesCRdBD2d_stressP2(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;
	int dop = elementDOF->dop;
	int col = (dop + 1)*(dop + 2) / 2;
	create_dbd_matrix(elementDOF->dof, elementDOF->dof, elements->row*3, M);
	for (i = 0; i<M->nb; i++)
		create_dden_matrix(col, col, M->blk + i);

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;

	int num_qp21 = getNumQuadPoints(dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<col; k1++)
		{
			for (k2 = 0; k2<col; k2++)
			{
				for (i1 = 0; i1<num_qp21; i1++)
				{
					lambdas[0] = gauss21[i1][0];
					lambdas[1] = gauss21[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, elementDOF[0].dop, &phi1[0]);
					lagrange_basis(lambdas, k2, elementDOF[0].dop, &phi2[0]);
					M->blk[3 * k].val[k1][k2] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] / (2 * mu);
				}
				M->blk[3 * k + 1].val[k1][k2] = M->blk[3 * k].val[k1][k2];
				M->blk[3 * k + 2].val[k1][k2] = M->blk[3 * k].val[k1][k2] * 2;
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesCR_psv_stressCompact(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_stokesCR_psv_stressCompact(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrix_stokesCR_psv_stressCompact(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);
	init_dvector(&uh[1], 1);///////////////////////////////////

	//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[2], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[2]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);

	/////////////////////
	/*	int i,j;
	for (i = 0; i < ptr_A[1].row; i++)
	{
	for (j = ptr_A[1].IA[i]; j < ptr_A[1].IA[i + 1]; j++)
	{
	printf("( %d, %d, %e); ", i, ptr_A[1].JA[j], ptr_A[1].val[j]);
	}
	printf("\n");
	}*/
	//	exit(1);
}

/**
* \fn void assembleStiffmatrix_stokesCR_psv_stressCompact(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_stokesCR_psv_stressCompact(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	if (t < 0)
		t = 0;

	A->row = elementDOF[0].dof + elementDOF[1].dof;
	A->col = A->row;
	A->IA = (int*)malloc((A->row + 1) * sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[2].dof * 2;
	B->col = A->col;
	B->IA = (int*)malloc((B->row + 1) * sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
									  //	int *index;
	int istart;
	//	index = (int*)calloc(A->col, sizeof(int));
	//	for (i = 0; i<A->col; i++)
	//		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	A->IA[0] = 0;
	for (i = 0; i < elementDOF[0].dof; i++)
	{
		A->IA[i + 1] = 11; // elementDOF[0].col + 8;
	}
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			if (k1 < 6)
				A->IA[i + 1] = 7;
			else if (k1<9)
				A->IA[i + 1] = 5;
			else
				A->IA[i + 1] = 11;

		}
	}

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)malloc(A->nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		// pressure p
		for (i = 0; i < elementDOF[0].col; i++)
		{
			rowstart[0] = A->IA[elementDOF[0].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			for (j = 0; j < 6; j++)
				A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			A->JA[rowstart[0] + elementDOF[0].col + 6] = elementDOF[1].val[k][9] + elementDOF[0].dof;
			A->JA[rowstart[0] + elementDOF[0].col + 7] = elementDOF[1].val[k][10] + elementDOF[0].dof;
		}

		// stress sigma
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = A->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 3)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 0; j < 3; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				A->JA[rowstart[0] + elementDOF[0].col + 3] = elementDOF[1].val[k][9] + elementDOF[0].dof;
			}
			else if (i < 6)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 0; j < 3; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
				A->JA[rowstart[0] + elementDOF[0].col + 3] = elementDOF[1].val[k][10] + elementDOF[0].dof;
			}
			else if (i < 9)
			{
				for (j = 0; j < 5; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else if (i == 9)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 0; j < 3; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				for (j = 6; j < 11; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
			else // (i == 10)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 3; j < 11; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
		} // i
	} // k
	  //	free(index);

	  // step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			// A11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + elementDOF[0].col; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * (1.0 / mu + t);
						}
						break;
					}
					if (j1 == A->IA[i] + elementDOF[0].col)
					{
						printf("There is something wrong in constructing A11\n");
						exit(1);
					}
				}
			} // k2
			  // A12
			for (k2 = 0; k2 < elementDOF[1].col; k2++)
			{
				if (k2 >= 6 && k2<9)
					continue;

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + elementDOF[0].col; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi1[0] * (phi2[0] + phi2[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A12\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			row21[1] = 3;
			if (k1 >= 6 && k1<9)
				row21[1] = 0;

			// A21
			if (k1 < 6 || k1 >= 9)
			{
				count = 3;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
			}
			else
				count = 0;

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi2[0] * (phi1[0] + phi1[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i] + row21[1])
					{
						printf("There is something wrong in constructing A21\n");
						exit(1);
					}
				}
			} // k2

			  // A22
			if (k1 < 3)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				patchnodes[3] = 9;
			}
			else if (k1 < 6)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = 3 + j;
				patchnodes[3] = 10;
			}
			else if (k1 < 9)
			{
				count = 5;
				for (j = 0; j < 5; j++)
					patchnodes[j] = 6 + j;
			}
			else if (k1 == 9)
			{
				count = 8;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				for (j = 3; j < 8; j++)
					patchnodes[j] = 3 + j;
			}
			else // k1==10 
			{
				count = 8;
				for (j = 0; j < 8; j++)
					patchnodes[j] = 3 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + row21[1]; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	B->IA[0] = 0;
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 8;
		B->IA[i + elementDOF[2].dof + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 8;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)malloc(B->nnz * sizeof(int));
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			j1 = (j - elementdofTran->IA[i]) * 8;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[2].dof] + j1;
			for (i1 = 0; i1 < 3; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1] + elementDOF[0].dof;
			for (i1 = 3; i1 < 8; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][3 + i1] + elementDOF[0].dof;
			for (i1 = 0; i1 < 8; i1++)
				B->JA[rowstart[1] + i1] = elementDOF[1].val[element[0]][3 + i1] + elementDOF[0].dof;
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		// b11
		count = 8;
		for (i = 0; i < 3; i++)
			patchnodes[i] = i;
		for (i = 3; i < 8; i++)
			patchnodes[i] = 3 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		  // b21
		for (i = 0; i < 8; i++)
			patchnodes[i] = 3 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1] + elementDOF[2].dof;
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[1] * phi2[1] + phi1[0] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				cr_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] -= 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b2->val[i + elementDOF[2].dof] -= 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesCR2d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesCR2d_stressCompact(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	M->row = elementDOF[0].dof + elementDOF[1].dof;
	M->col = M->row;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial


									  /************************************************** stiffness matrix A *****************************************************************/
									  //	int *index;
	int istart;
	//	index = (int*)calloc(A->col, sizeof(int));
	//	for (i = 0; i<A->col; i++)
	//		index[i] = -1;
	// step 1M: Find first the structure IA of the stiffness matrix M
	M->IA[0] = 0;
	for (i = 0; i < elementDOF[0].dof; i++)
		M->IA[i + 1] = 3; // elementDOF[0].col;

	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			if (k1 < 6)
				M->IA[i + 1] = 4;
			else if (k1<9)
				M->IA[i + 1] = 5;
			else
				M->IA[i + 1] = 8;
		}
	}

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];

	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)malloc(M->nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		// pressure p
		for (i = 0; i < elementDOF[0].col; i++)
		{
			rowstart[0] = M->IA[elementDOF[0].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				M->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
		}

		// stress sigma
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = M->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 3)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				M->JA[rowstart[0] + 3] = elementDOF[1].val[k][9] + elementDOF[0].dof;
			}
			else if (i < 6)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
				M->JA[rowstart[0] + 3] = elementDOF[1].val[k][10] + elementDOF[0].dof;
			}
			else if (i < 9)
			{
				for (j = 0; j < 5; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else if (i == 9)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				for (j = 3; j < 8; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
			}
			else // (i == 10)
			{
				for (j = 0; j < 8; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
			}
		} // i
	} // k
	  //	free(index);

	  // step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			// M11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * (1.0 / mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M11\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			// M22
			if (k1 < 3)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				patchnodes[3] = 9;
			}
			else if (k1 < 6)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = 3 + j;
				patchnodes[3] = 10;
			}
			else if (k1 < 9)
			{
				count = 5;
				for (j = 0; j < 5; j++)
					patchnodes[j] = 6 + j;
			}
			else if (k1 == 9)
			{
				count = 8;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				for (j = 3; j < 8; j++)
					patchnodes[j] = 3 + j;
			}
			else // k1==10 
			{
				count = 8;
				for (j = 0; j < 8; j++)
					patchnodes[j] = 3 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesCRdBD2d_stressCompact(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesCRdBD2d_stressCompact(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;
	int dop1 = elementDOF[0].dop;
//	int dop2 = elementDOF[1].dop;
	int col1 = (dop1 + 1)*(dop1 + 2) / 2;
	int col2 = 11;
	create_dbd_matrix(elementDOF[0].dof + elementDOF[1].dof, elementDOF[0].dof + elementDOF[1].dof, elements->row * 2, M);
	for (i = 0; i<elements->row; i++)
		create_dden_matrix(col1, col1, M->blk + i);
	for (i = elements->row; i<M->nb; i++)
		create_dden_matrix(col2, col2, M->blk + i);

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;

	int num_qp21 = getNumQuadPoints(dop1 * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<col1; k1++)
		{
			for (k2 = 0; k2<col1; k2++)
			{
				for (i1 = 0; i1<num_qp21; i1++)
				{
					lambdas[0] = gauss21[i1][0];
					lambdas[1] = gauss21[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, elementDOF[0].dop, &phi1[0]);
					lagrange_basis(lambdas, k2, elementDOF[0].dop, &phi2[0]);
					M->blk[k].val[k1][k2] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] / mu;
				}
			} // k2
		} // k1

		for (k1 = 0; k1<col2; k1++)
		{
			for (k2 = 0; k2<col2; k2++)
			{
				for (i1 = 0; i1<num_qp22; i1++)
				{
					lambdas[0] = gauss22[i1][0];
					lambdas[1] = gauss22[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
					crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
					M->blk[elements->row + k].val[k1][k2] += 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesCR_psv_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_stokesCR_psv_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrix_stokesCR_psv_stressP2(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);
	
	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);
	init_dvector(&uh[1], 1);///////////////////////////////////

	//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[2], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[2]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);

	/////////////////////
	/*	int i,j;
	for (i = 0; i < ptr_A[1].row; i++)
	{
	for (j = ptr_A[1].IA[i]; j < ptr_A[1].IA[i + 1]; j++)
	{
	printf("( %d, %d, %e); ", i, ptr_A[1].JA[j], ptr_A[1].val[j]);
	}
	printf("\n");
	}*/
	//	exit(1);
}

/**
* \fn void assembleStiffmatrix_stokesCR_psv_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_stokesCR_psv_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	if (t < 0)
		t = 0;

	A->row = elementDOF[0].dof + elementDOF[1].dof;
	A->col = A->row;
	A->IA = (int*)malloc((A->row + 1) * sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[2].dof * 2;
	B->col = A->col;
	B->IA = (int*)malloc((B->row + 1) * sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp20 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss20[num_qp20][3];
	init_Gauss(num_qp20, 2, gauss20); // gauss intergation initial

	int num_qp21 = getNumQuadPoints(3, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
									  //	int *index;
	int istart;
	//	index = (int*)calloc(A->col, sizeof(int));
	//	for (i = 0; i<A->col; i++)
	//		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	A->IA[0] = 0;
	for (i = 0; i < elementDOF[0].dof; i++)
	{
		A->IA[i + 1] = 3 + 6 * 2; // elementDOF[0].col + 8;
	}
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			if (k1 < 6 * 2)
				A->IA[i + 1] = 3 + 6;
			else
				A->IA[i + 1] = 6;

		}
	}

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)malloc(A->nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		// pressure p
		for (i = 0; i < elementDOF[0].col; i++)
		{
			rowstart[0] = A->IA[elementDOF[0].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			for (j = 0; j < 6 * 2; j++)
				A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
		}

		// stress sigma
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = A->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 6)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 0; j < 6; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
			else if (i < 12)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 0; j < 6; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else
			{
				for (j = 0; j < 6; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][12 + j] + elementDOF[0].dof;
			}
		} // i
	} // k
	  //	free(index);

	  // step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			// A11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + elementDOF[0].col; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp20; i1++)
						{
							lambdas[0] = gauss20[i1][0];
							lambdas[1] = gauss20[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] += 2 * s*gauss20[i1][2] * phi1[0] * phi2[0] * (1.0 / mu + t);
						}
						break;
					}
					if (j1 == A->IA[i] + elementDOF[0].col)
					{
						printf("There is something wrong in constructing A11\n");
						exit(1);
					}
				}
			} // k2
			  // A12
			for (k2 = 0; k2 < 6 * 2; k2++)
			{
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + elementDOF[0].col; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//						crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi1[0] * (phi2[0] + phi2[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A12\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			row21[1] = 3;
			if (k1 >= 12)
				row21[1] = 0;

			// A21
			if (k1 < 12)
			{
				count = 3;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
			}
			else
				count = 0;

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							//				crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi2[0] * (phi1[0] + phi1[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i] + row21[1])
					{
						printf("There is something wrong in constructing A21\n");
						exit(1);
					}
				}
			} // k2

			  // A22
			if (k1 < 6)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = j;
			}
			else if (k1 < 12)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 6 + j;
			}
			else
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 12 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + row21[1]; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//			crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							//			crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	B->IA[0] = 0;
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 12;
		B->IA[i + elementDOF[2].dof + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i]) * 12;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)malloc(B->nnz * sizeof(int));
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			j1 = (j - elementdofTran->IA[i]) * 12;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[2].dof] + j1;
			for (i1 = 0; i1 < 6; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1] + elementDOF[0].dof;
			for (i1 = 6; i1 < 12; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][6 + i1] + elementDOF[0].dof;
			for (i1 = 0; i1 < 12; i1++)
				B->JA[rowstart[1] + i1] = elementDOF[1].val[element[0]][6 + i1] + elementDOF[0].dof;
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		// b11
		count = 12;
		for (i = 0; i < 6; i++)
			patchnodes[i] = i;
		for (i = 6; i < 12; i++)
			patchnodes[i] = 6 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//			crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		  // b21
		for (i = 0; i < 12; i++)
			patchnodes[i] = 6 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1] + elementDOF[2].dof;
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							cr_basis1(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//		crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[1] * phi2[1] + phi1[0] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				cr_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] -= 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b2->val[i + elementDOF[2].dof] -= 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesCR2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesCR2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;

	M->row = elementDOF[0].dof + elementDOF[1].dof;
	M->col = M->row;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
									  //	int *index;
	int istart;
	//	index = (int*)calloc(A->col, sizeof(int));
	//	for (i = 0; i<A->col; i++)
	//		index[i] = -1;
	// step 1M: Find first the structure IA of the stiffness matrix M
	M->IA[0] = 0;
	for (i = 0; i < elementDOF[0].dof; i++)
		M->IA[i + 1] = 3; // elementDOF[0].col;

	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			M->IA[i + 1] = 6;
		}
	}

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];

	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)malloc(M->nnz * sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		// pressure p
		for (i = 0; i < elementDOF[0].col; i++)
		{
			rowstart[0] = M->IA[elementDOF[0].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				M->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
		}

		// stress sigma
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = M->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 6)
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
			else if (i < 12)
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][12 + j] + elementDOF[0].dof;
			}
		} // i
	} // k
	  //	free(index);

	  // step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			// M11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * (1.0 / mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M11\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			// M22
			if (k1 < 6)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = j;
			}
			else if (k1 < 12)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 6 + j;
			}
			else
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 12 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//	crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							//	crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							M->val[j1] += 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesCRdBD2d_stressP2(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of Crouzeix每Raviart element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesCRdBD2d_stressP2(dBDmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;
	int dop1 = elementDOF[0].dop;
	int dop2 = elementDOF[1].dop;
	int col1 = (dop1 + 1)*(dop1 + 2) / 2;
	int col2 = (dop2 + 1)*(dop2 + 2) / 2;
	create_dbd_matrix(elementDOF[0].dof + elementDOF[1].dof, elementDOF[0].dof + elementDOF[1].dof, elements->row * 4, M);
	for (i = 0; i<elements->row; i++)
		create_dden_matrix(col1, col1, M->blk + i);
	for (i = elements->row; i<M->nb; i++)
		create_dden_matrix(col2, col2, M->blk + i);

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;

	int num_qp21 = getNumQuadPoints(dop1 * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(dop2 * 2, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<col1; k1++)
		{
			for (k2 = 0; k2<col1; k2++)
			{
				for (i1 = 0; i1<num_qp21; i1++)
				{
					lambdas[0] = gauss21[i1][0];
					lambdas[1] = gauss21[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, elementDOF[0].dop, &phi1[0]);
					lagrange_basis(lambdas, k2, elementDOF[0].dop, &phi2[0]);
					M->blk[k].val[k1][k2] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] / mu;
				}
			} // k2
		} // k1

		for (k1 = 0; k1<col2; k1++)
		{
			for (k2 = 0; k2<col2; k2++)
			{
				for (i1 = 0; i1<num_qp22; i1++)
				{
					lambdas[0] = gauss22[i1][0];
					lambdas[1] = gauss22[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi1[0]);
					lagrange_basis(lambdas, k2, elementDOF[1].dop, &phi2[0]);
					M->blk[elements->row + 3 * k].val[k1][k2] += 2 * s*gauss22[i1][2] * phi1[0] * phi2[0] / (2 * mu);
				}
				M->blk[elements->row + 3 * k + 1].val[k1][k2] = M->blk[elements->row + 3 * k].val[k1][k2];
				M->blk[elements->row + 3 * k + 2].val[k1][k2] = M->blk[elements->row + 3 * k].val[k1][k2] * 2;
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesMINI_psv_stressP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_stokesMINI_psv_stressP1P0(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrix_stokesMINI_psv_stressP1P0(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);
	init_dvector(&uh[1], 1);///////////////////////////////////

							//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[2], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[2]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrix_stokesMINI_psv_stressP1P0(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix of MINI element for Stokes equation with stress discretized by piecewise P1-P0 element with bubble
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_stokesMINI_psv_stressP1P0(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	if (t < 0)
		t = 0;

	A->row = elementDOF[0].dof + elementDOF[1].dof;
	A->col = A->row;
	A->IA = (int*)malloc((A->row + 1) * sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[2].dof * 2;
	B->col = A->col;
	B->IA = (int*)malloc((B->row + 1) * sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp20 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss20[num_qp20][3];
	init_Gauss(num_qp20, 2, gauss20); // gauss intergation initial

	int num_qp21 = getNumQuadPoints(3, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	A->IA[0] = 0;
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;
		A->IA[i + 1] += 8 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]);

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			if (k1 < 6)
				A->IA[i + 1] = 3 + 4;
			else if (k1 < 7)
				A->IA[i + 1] = 0 + 3;
			else
				A->IA[i + 1] = 3 + 6;

		}
	}

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)malloc(A->nnz * sizeof(int));
	// pressure p
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}

			//			j1 = (j - elementdofTran[0].IA[i]) * 8;
			//			rowstart[0] = A->IA[i + 1] - 8 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]) + j1;
			rowstart[0] = A->IA[i + 1] + 8 * (j - elementdofTran[0].IA[i + 1]);
			for (i1 = 0; i1 < 6; i1++)
				A->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1] + elementDOF[0].dof;
			A->JA[rowstart[0] + 6] = elementDOF[1].val[element[0]][7] + elementDOF[0].dof;
			A->JA[rowstart[0] + 7] = elementDOF[1].val[element[0]][8] + elementDOF[0].dof;
		}

		for (j = A->IA[i]; j<A->IA[i + 1] - 8 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]); j++)
		{
			A->JA[j] = istart;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);
	// stress sigma
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = A->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 3)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
				for (j = 0; j < 3; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				A->JA[rowstart[0] + elementDOF[0].col + 3] = elementDOF[1].val[k][7] + elementDOF[0].dof;
			}
			else if (i < 6)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
				for (j = 0; j < 3; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
				A->JA[rowstart[0] + elementDOF[0].col + 3] = elementDOF[1].val[k][8] + elementDOF[0].dof;
			}
			else if (i < 7)
			{
				for (j = 0; j < 3; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else if (i == 7)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 0; j < 3; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				for (j = 6; j < 9; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
			else // (i == 8)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 3; j < 9; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
		} // i
	} // k
	  //	free(index);

	  // step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			row21[0] = A->IA[i + 1] - A->IA[i] - 8 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]);
			// A11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp20; i1++)
						{
							lambdas[0] = gauss20[i1][0];
							lambdas[1] = gauss20[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] += 2 * s*gauss20[i1][2] * phi1[0] * phi2[0] * (1.0 / mu + t);
						}
						break;
					}
					if (j1 == A->IA[i] + row21[0])
					{
						printf("There is something wrong in constructing A11\n");
						exit(1);
					}
				}
			} // k2
			  // A12
			for (k2 = 0; k2 < 9; k2++)
			{
				if (k2 == 6)
					continue;

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + row21[0]; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi1[0] * (phi2[0] + phi2[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A12\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			row21[1] = 3;
			if (k1 == 6)
				row21[1] = 0;

			// A21
			if (k1 != 6)
			{
				count = 3;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
			}
			else
				count = 0;

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi2[0] * (phi1[0] + phi1[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i] + row21[1])
					{
						printf("There is something wrong in constructing A21\n");
						exit(1);
					}
				}
			} // k2

			  // A22
			if (k1 < 3)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				patchnodes[3] = 7;
			}
			else if (k1 < 6)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = 3 + j;
				patchnodes[3] = 8;
			}
			else if (k1 < 7)
			{
				count = 3;
				for (j = 0; j < 3; j++)
					patchnodes[j] = 6 + j;
			}
			else if (k1 == 7)
			{
				count = 6;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				for (j = 3; j < 6; j++)
					patchnodes[j] = 3 + j;
			}
			else // k1==8 
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 3 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + row21[1]; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k1, phi1);
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	B->IA[0] = 0;
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran[1].IA[i + 1] - elementdofTran[1].IA[i]) * 6;
		B->IA[i + elementDOF[2].dof + 1] = (elementdofTran[1].IA[i + 1] - elementdofTran[1].IA[i]) * 6;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)malloc(B->nnz * sizeof(int));
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		for (j = elementdofTran[1].IA[i]; j<elementdofTran[1].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[1].JA[j];
			j1 = (j - elementdofTran[1].IA[i]) * 6;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[2].dof] + j1;
			for (i1 = 0; i1 < 3; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1] + elementDOF[0].dof;
			for (i1 = 3; i1 < 6; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][3 + i1] + elementDOF[0].dof;
			for (i1 = 0; i1 < 6; i1++)
				B->JA[rowstart[1] + i1] = elementDOF[1].val[element[0]][3 + i1] + elementDOF[0].dof;
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		// b11
		count = 6;
		for (i = 0; i < 3; i++)
			patchnodes[i] = i;
		for (i = 3; i < 6; i++)
			patchnodes[i] = 3 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		  // b21
		for (i = 0; i < 6; i++)
			patchnodes[i] = 3 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1] + elementDOF[2].dof;
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[1] * phi2[1] + phi1[0] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				mini_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] -= 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b2->val[i + elementDOF[2].dof] -= 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesMINI2d_stressP1P0(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of MINI element for Stokes equation with stress discretized by piecewise P1P0 element
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesMINI2d_stressP1P0(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i, j, k, l;

	M->row = elementDOF[0].dof + elementDOF[1].dof;
	M->col = M->row;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial


									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(M->col * sizeof(int));
	for (i = 0; i<M->col; i++)
		index[i] = -1;
	// step 1M: Find first the structure IA of the stiffness matrix M
	M->IA[0] = 0;
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		M->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			if (k1 < 6)
				M->IA[i + 1] = 4;
			else if (k1 < 7)
				M->IA[i + 1] = 3;
			else
				M->IA[i + 1] = 6;
		}
	}

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];

	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)malloc(M->nnz * sizeof(int));
	// pressure p
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = M->IA[i]; j<M->IA[i + 1]; j++)
		{
			M->JA[j] = istart;
			istart = index[istart];
			index[M->JA[j]] = -1;
		}
	} // i
	free(index);
	// stress sigma
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = M->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 3)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				M->JA[rowstart[0] + 3] = elementDOF[1].val[k][7] + elementDOF[0].dof;
			}
			else if (i < 6)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
				M->JA[rowstart[0] + 3] = elementDOF[1].val[k][8] + elementDOF[0].dof;
			}
			else if (i < 7)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else if (i == 7)
			{
				for (j = 0; j < 3; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
				for (j = 3; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
			}
			else // (i == 8)
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][3 + j] + elementDOF[0].dof;
			}
		} // i
	} // k

	  // step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			// M11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * (2.0 / (2 * mu));
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M11\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			// M22
			if (k1 < 3)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				patchnodes[3] = 7;
			}
			else if (k1 < 6)
			{
				count = 4;
				for (j = 0; j < 3; j++)
					patchnodes[j] = 3 + j;
				patchnodes[3] = 8;
			}
			else if (k1 < 7)
			{
				count = 3;
				for (j = 0; j < 3; j++)
					patchnodes[j] = 6 + j;
			}
			else if (k1 == 7)
			{
				count = 6;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
				for (j = 3; j < 6; j++)
					patchnodes[j] = 3 + j;
			}
			else // k1==8 
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 3 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k1, phi1);
							miniSymTensorP1P0_basis(lambdas, s, eta, xi, k2, phi2);
							M->val[j1] += 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesMINIdBD2d_stressP1P0(dBDmat *M, dvector *diag, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of MINI element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesMINIdBD2d_stressP1P0(dBDmat *M, dvector *diag, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;
	int dof1 = elementDOF[0].dof;
	int dof2 = elementDOF[1].dof;
	//	int dop = elementDOF[1].dop;
	int col = 9;
	create_dbd_matrix(dof1 + dof2, dof1 + dof2, dof1 + elements->row, M);
	for (i = 0; i<dof1; i++)
		create_dden_matrix(1, 1, M->blk + i);
	for (i = dof1; i<M->nb; i++)
		create_dden_matrix(col, col, M->blk + i);

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;

	int num_qp21 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	for (i = 0; i < dof1; i++)
		M->blk[i].val[0][0] = diag->val[i];

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<col; k1++)
		{
			for (k2 = 0; k2<col; k2++)
			{
				for (i1 = 0; i1<num_qp21; i1++)
				{
					lambdas[0] = gauss21[i1][0];
					lambdas[1] = gauss21[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					miniSymTensorP1P0_basis(lambdas, s, eta, xi, k1, phi1);
					miniSymTensorP1P0_basis(lambdas, s, eta, xi, k2, phi2);
					M->blk[dof1 + k].val[k1][k2] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesMINI_psv_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_stokesMINI_psv_stressP2(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrix_stokesMINI_psv_stressP2(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);
	init_dvector(&uh[1], 1);///////////////////////////////////

							//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[2], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[2]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrix_stokesMINI_psv_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix of MINI element for Stokes equation with stress discretized by piecewise P2 element
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_stokesMINI_psv_stressP2(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	if (t < 0)
		t = 0;

	A->row = elementDOF[0].dof + elementDOF[1].dof;
	A->col = A->row;
	A->IA = (int*)malloc((A->row + 1) * sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[2].dof * 2;
	B->col = A->col;
	B->IA = (int*)malloc((B->row + 1) * sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp20 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss20[num_qp20][3];
	init_Gauss(num_qp20, 2, gauss20); // gauss intergation initial

	int num_qp21 = getNumQuadPoints(3, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(A->col * sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	A->IA[0] = 0;
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;
		A->IA[i + 1] += 12 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]);

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i
	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			if (k1 < 6 * 2)
				A->IA[i + 1] = 3 + 6;
			else
				A->IA[i + 1] = 6;

		}
	}

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)malloc(A->nnz * sizeof(int));
	// pressure p
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}

			//			j1 = (j - elementdofTran[0].IA[i]) * 12;
			//			rowstart[0] = A->IA[i + 1] - 12 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]) + j1;
			rowstart[0] = A->IA[i + 1] + 12 * (j - elementdofTran[0].IA[i + 1]);
			for (i1 = 0; i1 < 12; i1++)
				A->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1] + elementDOF[0].dof;
		}

		for (j = A->IA[i]; j<A->IA[i + 1] - 12 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]); j++)
		{
			A->JA[j] = istart;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);
	// stress sigma
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = A->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 6)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
				for (j = 0; j < 6; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
			else if (i < 12)
			{
				for (j = 0; j < elementDOF[0].col; j++)
					A->JA[rowstart[0] + j] = elementDOF[0].val[k][j];

				for (j = 0; j < 6; j++)
					A->JA[rowstart[0] + elementDOF[0].col + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else
			{
				for (j = 0; j < 6; j++)
					A->JA[rowstart[0] + j] = elementDOF[1].val[k][12 + j] + elementDOF[0].dof;
			}
		} // i
	} // k
	  //	free(index);

	  // step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			row21[0] = A->IA[i + 1] - A->IA[i] - 12 * (elementdofTran[0].IA[i + 1] - elementdofTran[0].IA[i]);
			// A11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp20; i1++)
						{
							lambdas[0] = gauss20[i1][0];
							lambdas[1] = gauss20[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] += 2 * s*gauss20[i1][2] * phi1[0] * phi2[0] * (1.0 / mu + t);
						}
						break;
					}
					if (j1 == A->IA[i] + row21[0])
					{
						printf("There is something wrong in constructing A11\n");
						exit(1);
					}
				}
			} // k2
			  // A12
			for (k2 = 0; k2 < 6 * 2; k2++)
			{
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + row21[0]; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//						crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi1[0] * (phi2[0] + phi2[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A12\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			row21[1] = 3;
			if (k1 >= 12)
				row21[1] = 0;

			// A21
			if (k1 < 12)
			{
				count = 3;
				for (j = 0; j < 3; j++)
					patchnodes[j] = j;
			}
			else
				count = 0;

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1 < A->IA[i] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							//				crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							A->val[j1] -= 2 * s*gauss21[i1][2] * phi2[0] * (phi1[0] + phi1[1]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i] + row21[1])
					{
						printf("There is something wrong in constructing A21\n");
						exit(1);
					}
				}
			} // k2

			  // A22
			if (k1 < 6)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = j;
			}
			else if (k1 < 12)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 6 + j;
			}
			else
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 12 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = A->IA[i] + row21[1]; j1 < A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//			crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							//			crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							A->val[j1] += 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == A->IA[i + 1])
					{
						printf("There is something wrong in constructing A22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	B->IA[0] = 0;
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran[1].IA[i + 1] - elementdofTran[1].IA[i]) * 12;
		B->IA[i + elementDOF[2].dof + 1] = (elementdofTran[1].IA[i + 1] - elementdofTran[1].IA[i]) * 12;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)malloc(B->nnz * sizeof(int));
	for (i = 0; i<elementDOF[2].dof; i++)
	{
		for (j = elementdofTran[1].IA[i]; j<elementdofTran[1].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[1].JA[j];
			j1 = (j - elementdofTran[1].IA[i]) * 12;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[2].dof] + j1;
			for (i1 = 0; i1 < 6; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][i1] + elementDOF[0].dof;
			for (i1 = 6; i1 < 12; i1++)
				B->JA[rowstart[0] + i1] = elementDOF[1].val[element[0]][6 + i1] + elementDOF[0].dof;
			for (i1 = 0; i1 < 12; i1++)
				B->JA[rowstart[1] + i1] = elementDOF[1].val[element[0]][6 + i1] + elementDOF[0].dof;
		} // j
	} // i

	  // step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		// b11
		count = 12;
		for (i = 0; i < 6; i++)
			patchnodes[i] = i;
		for (i = 6; i < 12; i++)
			patchnodes[i] = 6 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//			crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		  // b21
		for (i = 0; i < 12; i++)
			patchnodes[i] = 6 + i;

		for (k1 = 0; k1 < elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1] + elementDOF[2].dof;
			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];
				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							mini_basis1(lambdas, s, eta, xi, k1, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//		crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							B->val[j1] -= 2 * s*gauss22[i1][2] * (phi1[1] * phi2[1] + phi1[0] * phi2[2]);
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1

	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	/************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[2].col; k1++)
		{
			i = elementDOF[2].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				mini_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] -= 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b2->val[i + elementDOF[2].dof] -= 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesMINI2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of MINI element for Stokes equation with stress discretized by piecewise P2 element
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesMINI2d_stressP2(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i, j, k, l;

	M->row = elementDOF[0].dof + elementDOF[1].dof;
	M->col = M->row;
	M->IA = (int*)malloc((M->row + 1) * sizeof(int));
	M->JA = NULL;
	M->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count, patchnodes[15];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial


									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)malloc(M->col * sizeof(int));
	for (i = 0; i<M->col; i++)
		index[i] = -1;
	// step 1M: Find first the structure IA of the stiffness matrix M
	M->IA[0] = 0;
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		M->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (k = 0; k < elements->row; k++)
	{
		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;
			M->IA[i + 1] = 6;
		}
	}

	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];

	M->nnz = M->IA[M->row];

	// step 2M: Find the structure JA of the stiffness matrix M
	M->JA = (int*)malloc(M->nnz * sizeof(int));
	// pressure p
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran[0].IA[i]; j<elementdofTran[0].IA[i + 1]; j++)
		{
			element[0] = elementdofTran[0].JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = M->IA[i]; j<M->IA[i + 1]; j++)
		{
			M->JA[j] = istart;
			istart = index[istart];
			index[M->JA[j]] = -1;
		}
	} // i
	free(index);
	// stress sigma
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = M->IA[elementDOF[1].val[k][i] + elementDOF[0].dof];
			if (i < 6)
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][j] + elementDOF[0].dof;
			}
			else if (i < 12)
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][6 + j] + elementDOF[0].dof;
			}
			else
			{
				for (j = 0; j < 6; j++)
					M->JA[rowstart[0] + j] = elementDOF[1].val[k][12 + j] + elementDOF[0].dof;
			}
		} // i
	} // k

	  // step 3M: Loop element by element and compute the actual entries storing them in M
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			// M11
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, 1, phi1);
							lagrange_basis(lambdas, k2, 1, phi2);
							M->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] * (2.0 / (2 * mu));
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M11\n");
						exit(1);
					}
				}
			} // k2
		} // k1

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1] + elementDOF[0].dof;

			// M22
			if (k1 < 6)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = j;
			}
			else if (k1 < 12)
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 6 + j;
			}
			else
			{
				count = 6;
				for (j = 0; j < 6; j++)
					patchnodes[j] = 12 + j;
			}

			for (i2 = 0; i2 < count; i2++)
			{
				k2 = patchnodes[i2];

				j = elementDOF[1].val[k][k2] + elementDOF[0].dof;
				for (j1 = M->IA[i]; j1 < M->IA[i + 1]; j1++)
				{
					if (M->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrangeSymTensor_basis(lambdas, k1, 2, phi1);
							lagrangeSymTensor_basis(lambdas, k2, 2, phi2);
							//	crSymTensor_basis(lambdas, s, eta, xi, k1, phi1);
							//	crSymTensor_basis(lambdas, s, eta, xi, k2, phi2);
							M->val[j1] += 2 * s*gauss22[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) / (2 * mu);
						}
						break;
					}
					if (j1 == M->IA[i + 1])
					{
						printf("There is something wrong in constructing M22\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemblePressStressMassmatrixStokesMINIdBD2d_stressP2(dBDmat *M, dvector *diag, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble pressure-stress mass matrix of MINI element for Stokes equation
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assemblePressStressMassmatrixStokesMINIdBD2d_stressP2(dBDmat *M, dvector *diag, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k, l;
	int dof1 = elementDOF[0].dof;
	int dof2 = elementDOF[1].dof;
	int dop = elementDOF[1].dop;
	int col = (dop + 1)*(dop + 2) / 2;
	create_dbd_matrix(dof1 + dof2, dof1 + dof2, dof1 + elements->row * 3, M);
	for (i = 0; i<dof1; i++)
		create_dden_matrix(1, 1, M->blk + i);
	for (i = dof1; i<M->nb; i++)
		create_dden_matrix(col, col, M->blk + i);

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, i2, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;

	int num_qp21 = getNumQuadPoints(dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	for (i = 0; i < dof1; i++)
		M->blk[i].val[0][0] = diag->val[i];

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<col; k1++)
		{
			for (k2 = 0; k2<col; k2++)
			{
				for (i1 = 0; i1<num_qp21; i1++)
				{
					lambdas[0] = gauss21[i1][0];
					lambdas[1] = gauss21[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, k1, dop, &phi1[0]);
					lagrange_basis(lambdas, k2, dop, &phi2[0]);
					M->blk[dof1 + 3 * k].val[k1][k2] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[0] / (2 * mu);
				}
				M->blk[dof1 + 3 * k + 1].val[k1][k2] = M->blk[dof1 + 3 * k].val[k1][k2];
				M->blk[dof1 + 3 * k + 2].val[k1][k2] = M->blk[dof1 + 3 * k].val[k1][k2] * 2;
			} // k2
		} // k1
	} // k
}

/**
* \fn void assemble_stokesNcP1P0_elas(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix *C and righ hand side *ptr_b
* \param *ptr_A pointer to stiffness matrix
* \param *ptr_b pointer to right hand side
* \param *uh pointer to solution, i.e. stress and displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assemble_stokesNcP1P0_elas(dCSRmat *ptr_A, dvector *ptr_b, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	/**	A[0]=A, A[1]=B^T, A[2]=B
	Ax1 + B^Tx2 = 0
	Bx1 +       = b
	where x1: sigma_h, x2: u_h
	**/
	dCSRmat A, B;
	dvector b;

	assembleStiffmatrix_stokesNcP1P0_elas(&ptr_A[0], &B, &ptr_A[3], &ptr_b[0], &b, elements, elementEdge, edges, nodes, elementDOF, elementdofTran, mu, t);

	// initial solution
	create_dvector(ptr_b[0].row, &uh[0]);
	create_dvector(b.row, &uh[1]);
	init_dvector(&uh[1], 1);///////////////////////////////////

	//extract
	extractFreenodesVector2ElasDirichlet(&ptr_A[0], &B, &b, ptr_b, &elementDOF[1], &uh[1]);
	free_dvector(&b);
	extractFreenodesMatrix1r(&B, &ptr_A[2], &elementDOF[1]);
	free_csr_matrix(&B);
	getTransposeOfSparse(&ptr_A[2], &ptr_A[1]);
}

/**
* \fn void assembleStiffmatrix_stokesNcP1P0_elas(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
									the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \param t 1/t=lambda Lame constant
* \return void
*/
void assembleStiffmatrix_stokesNcP1P0_elas(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu, double t)
{
	int i, j, k, l;

	if (t < 0)
		t = 0;

	A->row = elementDOF[0].dof * 4;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof * 2;
	B->col = A->col;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart, colstart;
	int count;

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
//	int *index;
	int istart;
//	index = (int*)calloc(A->col, sizeof(int));
//	for (i = 0; i<A->col; i++)
//		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (k = 0; k < elements->row; k++)
	{
		A->IA[4 * k + 1] = 2;
		A->IA[4 * k + 2] = 2;
		A->IA[4 * k + 3] = 1;
		A->IA[4 * k + 4] = 1;
	}

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		j = A->IA[4 * k];
		A->JA[j] = 4 * k;
		A->JA[j+1] = 4 * k+1;

		j = A->IA[4 * k + 1];
		A->JA[j] = 4 * k;
		A->JA[j + 1] = 4 * k + 1;

		j = A->IA[4 * k + 2];
		A->JA[j] = 4 * k + 2;
		
		j = A->IA[4 * k + 3];
		A->JA[j] = 4 * k + 3;
	}
//	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		// end set parameters

		j = A->IA[4 * k];
		A->val[j] = (1 - 1 / (2 + 2 * mu*t))*s / (2 * mu);
		A->val[j + 1] = ( - 1 / (2 + 2 * mu*t))*s / (2 * mu);

		j = A->IA[4 * k + 1];
		A->val[j] = (-1 / (2 + 2 * mu*t))*s / (2 * mu);
		A->val[j + 1] = (1 - 1 / (2 + 2 * mu*t))*s / (2 * mu);

		j = A->IA[4 * k + 2];
		A->val[j] = s / (2 * mu);

		j = A->IA[4 * k + 3];
		A->val[j] = s / (2 * mu);
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		B->IA[i + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i])*elementDOF[0].col * 2;
		B->IA[i + elementDOF[1].dof + 1] = (elementdofTran->IA[i + 1] - elementdofTran->IA[i])*elementDOF[0].col * 2;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];
			j1 = (j - elementdofTran->IA[i])*elementDOF[0].col * 2;
			rowstart[0] = B->IA[i] + j1;
			rowstart[1] = B->IA[i + elementDOF[1].dof] + j1;
			colstart = element[0] * elementDOF[0].col * 4;
			for (i1 = 0; i1 < elementDOF[0].col; i1++)
			{
				B->JA[rowstart[0] + i1] = colstart + i1;
				B->JA[rowstart[0] + elementDOF[0].col + i1] = colstart + elementDOF[0].col * 2 + i1;
				B->JA[rowstart[1] + i1] = colstart + elementDOF[0].col + i1;
				B->JA[rowstart[1] + elementDOF[0].col + i1] = colstart + elementDOF[0].col * 3 + i1;
			}
		} // j
	} // i

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
//				j = elementDOF[0].val[k][k2];
				// b11
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == k * elementDOF[0].col * 4 + k2)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[0];
							B->val[j1 + elementDOF[0].col] -= 2 * s*gauss22[i1][2] * phi1[1];
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
				// b21
				for (j1 = B->IA[i + elementDOF[1].dof]; j1 < B->IA[i + elementDOF[1].dof + 1]; j1++)
				{
					if (B->JA[j1] == k * elementDOF[0].col * 4 + elementDOF[0].col + k2)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							ncp1_basis1(lambdas, s, eta, xi, k1, phi1);
							B->val[j1] -= 2 * s*gauss22[i1][2] * phi1[1];
							B->val[j1 + elementDOF[0].col] -= 2 * s*gauss22[i1][2] * phi1[0];
						}
						break;
					}
					if (j1 == B->IA[i + elementDOF[1].dof + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix C *****************************************************************/
	C->row = 0;
	C->col = 0;
	C->IA = NULL;
	C->JA = NULL;
	C->val = NULL;

	  /************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial
	
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				ncp1_basis(lambdas, k1, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] -= 2 * s*gauss0[i1][2] * f1(x, y, 0, 0)*phi;
				b2->val[i + elementDOF[1].dof] -= 2 * s*gauss0[i1][2] * f2(x, y, 0, 0)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStressMassmatrixStokesNcP1P02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
* \brief assemble mass matrix
* \param *M pointer to mass matrix
* \param *elements pointer to the structure of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assembleStressMassmatrixStokesNcP1P02d(dCSRmat *M, ELEMENT *elements, ELEMENT_DOF *elementDOF, double mu)
{
	int i, j, k;

	M->row = elementDOF->dof * 4;
	M->col = elementDOF->dof * 4;
	M->IA = (int*)calloc(M->row + 1, sizeof(int));
	M->JA = NULL;
	M->val = NULL;
	for (i = 0; i<M->row; i++)
		M->IA[i + 1] = 1;
	for (i = 0; i<M->row; i++)
		M->IA[i + 1] += M->IA[i];
	M->nnz = M->IA[M->row];
	M->JA = (int*)calloc(M->nnz, sizeof(int));
	for (i = 0; i<M->nnz; i++)
		M->JA[i] = i;
	M->val = (double*)calloc(M->nnz, sizeof(double));
	for (k = 0; k < elements->row; k++)
	{
		M->val[4 * k] = elements->vol[k] / (2 * mu);
		M->val[4 * k + 1] = elements->vol[k] / (2 * mu);
		M->val[4 * k + 2] = elements->vol[k] / (2 * mu);
		M->val[4 * k + 3] = elements->vol[k] / (2 * mu);
	}
}

/**
* \fn void assembleStiffmatrixHuZhang2d(dCSRmat *A, dCSRmat *B, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b1 pointer to right hand side
* \param *b2 pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleStiffmatrixHuZhang2d(dCSRmat *A, dCSRmat *B, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof * 2;
	B->col = elementDOF[0].dof;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(B->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;
	double **tensorBasis[3];

	//		create_dden_matrix(3, 3, &tensorBasis[i]);

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{

				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k2, elementDOF[0].dop, phi2);
							if (lambda>-0.5)
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
							else
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		B->IA[i + 1] += elementDOF[0].col;
		B->IA[i + elementDOF[1].dof + 1] += elementDOF[0].col;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = B->IA[elementDOF[1].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			rowstart[0] = B->IA[elementDOF[1].val[k][i] + elementDOF[1].dof];
			for (j = 0; j < elementDOF[0].col; j++)
				B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// b11
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], tensorBasis, k2, elementDOF[0].dop, phi1);
							B->val[j1] += 2 * s*gauss22[i1][2] * phi1[0] * phi;
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
				// b21
				for (j1 = B->IA[i + elementDOF[1].dof]; j1 < B->IA[i + elementDOF[1].dof + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], tensorBasis, k2, elementDOF[0].dop, phi1);
							B->val[j1] += 2 * s*gauss22[i1][2] * phi1[1] * phi;
						}
						break;
					}
					if (j1 == B->IA[i + elementDOF[1].dof + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k


	  /************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	int num_qp10 = 20; // the number of numerical intergation points
	double gauss10[num_qp10][2];
	init_Gauss1D(num_qp10, 1, gauss10); // gauss intergation initial

	int patchnodes[200];
	double elen;

	for (edge = 0; edge<edges->row; edge++)
	{
		if (edges->bdFlag[edge] <1 || edges->bdFlag[edge]>4) //////////////////////////////////////////
			continue;

		// only consider Dirichlet boundary
		k = edges->val[edge][2];

		// set parameters
		elen = edges->length[edge];
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		for (l = 0; l<3; l++)
		{
			if (elementEdge->val[k][l] == edge)
				break;
		}

		l1 = (l + 1) % 3;
		if (elements->val[k][l1] != edges->val[edge][0])
			l1 = (l + 2) % 3;
		l2 = 3 - l - l1;

		count = getlocaldofsEdge4HuZhang2d(patchnodes, l, elementDOF[0].dop);

		for (i = 0; i<count; i++)
		{
			k1 = patchnodes[i];
			for (i1 = 0; i1<num_qp10; i1++)
			{
				lambdas[l1] = gauss10[i1][0];
				lambdas[l2] = 1 - gauss10[i1][0];
				lambdas[l] = 0;
				huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b1->val[elementDOF[0].val[k][k1]] += elen*gauss10[i1][1] * (phi1[0] * elements->nvector[k][l][0] + phi1[2] * elements->nvector[k][l][1])*u1(x, y, lambda, mu);
				b1->val[elementDOF[0].val[k][k1]] += elen*gauss10[i1][1] * (phi1[2] * elements->nvector[k][l][0] + phi1[1] * elements->nvector[k][l][1])*u2(x, y, lambda, mu);
			}
		}
	} // edge

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
				b2->val[i + elementDOF[1].dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k
}

/**
* \fn void assembleStiffmatrixHuZhangA11_2d(dCSRmat *A, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *BT pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleStiffmatrixHuZhangA11_2d(dCSRmat *A, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;
	double **tensorBasis[3];

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

										/************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{

				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k2, elementDOF[0].dop, phi2);
							if (lambda>-0.5)
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
							else
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k
}

/**
* \fn void assembleStiffmatrixHuZhangIP2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *B pointer to stiffness matrix
* \param *C pointer to stiffness matrix
* \param *b1 pointer to right hand side
* \param *b2 pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleStiffmatrixHuZhangIP2d(dCSRmat *A, dCSRmat *B, dCSRmat *C, dvector *b1, dvector *b2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double lambda, double mu)
{
	int i, j, k, l;

	A->row = elementDOF[0].dof;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	B->row = elementDOF[1].dof * 2;
	B->col = elementDOF[0].dof;
	B->IA = (int*)calloc(B->row + 1, sizeof(int));
	B->JA = NULL;
	B->val = NULL;

	C->row = elementDOF[1].dof * 2;
	C->col = C->row;
	C->IA = (int*)calloc(C->row + 1, sizeof(int));
	C->JA = NULL;
	C->val = NULL;

	create_dvector(A->row, b1);
	create_dvector(C->row, b2);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[3], phi2[3];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], taustart;
	int count;
	double **tensorBasis[3];

	int num_qp21 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp22 = getNumQuadPoints(elementDOF[0].dop + elementDOF[1].dop - 1, 2); // the number of numerical intergation points
	double gauss22[num_qp22][3];
	init_Gauss(num_qp22, 2, gauss22); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF[1].dop * 2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial

										/************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{
				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF[0].col; k++)
			{

				node = elementDOF[0].val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		//	xi=elements->xi[k];
		//	eta=elements->eta[k];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		for (k1 = 0; k1<elementDOF[0].col; k1++)
		{
			i = elementDOF[0].val[k][k1];
			for (k2 = 0; k2<elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i + 1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
							huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k2, elementDOF[0].dop, phi2);
							if (lambda>-0.5)
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - lambda / (2 * lambda + 2 * mu)*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
							else
								A->val[j1] += 2 * s*gauss21[i1][2] * ((phi1[0] * phi2[0] + phi1[1] * phi2[1] + 2 * phi1[2] * phi2[2]) - 1.0 / 2.0*(phi1[0] + phi1[1])*(phi2[0] + phi2[1])) / (2 * mu);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	  /************************************************** stiffness matrix B *****************************************************************/
	  // step 1B: Find first the structure IA of the stiffness matrix B
	for (i = 0; i<elementDOF[1].dof; i++)
	{
		B->IA[i + 1] += elementDOF[0].col;
		B->IA[i + elementDOF[1].dof + 1] += elementDOF[0].col;
	} // i

	for (i = 0; i<B->row; i++)
		B->IA[i + 1] += B->IA[i];

	B->nnz = B->IA[B->row];

	// step 2B: Find the structure JA of the stiffness matrix B
	B->JA = (int*)calloc(B->nnz, sizeof(int));
	for (k = 0; k < elements->row; k++)
	{
		for (i = 0; i < elementDOF[1].col; i++)
		{
			rowstart[0] = B->IA[elementDOF[1].val[k][i]];
			for (j = 0; j < elementDOF[0].col; j++)
				B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
			rowstart[0] = B->IA[elementDOF[1].val[k][i] + elementDOF[1].dof];
			for (j = 0; j < elementDOF[0].col; j++)
				B->JA[rowstart[0] + j] = elementDOF[0].val[k][j];
		}
	}

	// step 3B: Loop element by element and compute the actual entries storing them in B
	B->val = (double*)calloc(B->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		for (k1 = 0; k1 < elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (k2 = 0; k2 < elementDOF[0].col; k2++)
			{
				j = elementDOF[0].val[k][k2];
				// b11
				for (j1 = B->IA[i]; j1 < B->IA[i + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], tensorBasis, k2, elementDOF[0].dop, phi1);
							B->val[j1] += 2 * s*gauss22[i1][2] * phi1[0] * phi;
						}
						break;
					}
					if (j1 == B->IA[i + 1])
					{
						printf("There is something wrong in constructing b11 of B\n");
						exit(1);
					}
				}
				// b21
				for (j1 = B->IA[i + elementDOF[1].dof]; j1 < B->IA[i + elementDOF[1].dof + 1]; j1++)
				{
					if (B->JA[j1] == j)
					{
						for (i1 = 0; i1 < num_qp22; i1++)
						{
							lambdas[0] = gauss22[i1][0];
							lambdas[1] = gauss22[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
							huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], tensorBasis, k2, elementDOF[0].dop, phi1);
							B->val[j1] += 2 * s*gauss22[i1][2] * phi1[1] * phi;
						}
						break;
					}
					if (j1 == B->IA[i + elementDOF[1].dof + 1])
					{
						printf("There is something wrong in constructing b21 of B\n");
						exit(1);
					}
				}
			} // k2
		} // k1
	} // k


	  /************************************************** stiffness matrix C *****************************************************************/
	int curnode[2];
	// step 1C: Find first the structure IA of the stiffness matrix C
	for (k = 0; k<elements->row; k++)
	{
		if (elementDOF[1].dop == 0)
		{
			curnode[0] = elementDOF[1].val[k][0];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			C->IA[curnode[0] + 1] += 2;
			C->IA[curnode[1] + 1] += 2;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2;
					C->IA[curnode[1] + 1] += 2;
				}
			}
			continue;
		}

		for (i = 0; i<3 * elementDOF[1].dop; i++) //  for each node
		{
			curnode[0] = elementDOF[1].val[k][i];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			if (i<3)
			{
				C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop * 2 + 1);
				C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop * 2 + 1);
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
					C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				}
				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
					C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				}
			}
			else
			{
				C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
				C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				edge = elementEdge->val[k][(i - 3) / (elementDOF[1].dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					C->IA[curnode[0] + 1] += 2 * (elementDOF[1].dop + 1);
					C->IA[curnode[1] + 1] += 2 * (elementDOF[1].dop + 1);
				}
			}
		} // i
	} // k

	for (i = 0; i<C->row; i++)
		C->IA[i + 1] += C->IA[i];

	C->nnz = C->IA[C->row];

	// step 2C: Find the structure JA of the stiffness matrix C
	C->JA = (int*)calloc(C->nnz, sizeof(int));
	for (k = 0; k<elements->row; k++)
	{
		if (elementDOF[1].dop == 0)
		{
			curnode[0] = elementDOF[1].val[k][0];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = C->IA[curnode[j]];
				row21[j] = (C->IA[curnode[j] + 1] - C->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (j = 0; j<2; j++)
			{
				C->JA[rowstart[j] + count] = curnode[0];
				C->JA[rowstart[j] + count + row21[j]] = curnode[1];
			}
			count++;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = elementDOF[1].val[element[0]][0];
						C->JA[rowstart[j] + count + row21[j]] = elementDOF[1].val[element[0]][0] + elementDOF[1].dof;
					}
					count++;
				}
			}
			continue;
		}

		// elementDOF[1].dop > 0
		for (i = 0; i<3 * elementDOF[1].dop; i++) //  for each node
		{
			curnode[0] = elementDOF[1].val[k][i];
			curnode[1] = curnode[0] + elementDOF[1].dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = C->IA[curnode[j]];
				row21[j] = (C->IA[curnode[j] + 1] - C->IA[curnode[j]]) / 2;
			}
			count = 0;
			if (i<3)
			{
				for (i1 = 0; i1<3; i1++)
				{
					node = elementDOF[1].val[k][i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}
				l = (i + 1) % 3;
				for (i1 = 0; i1<elementDOF[1].dop - 1; i1++)
				{
					node = elementDOF[1].val[k][3 + l*(elementDOF[1].dop - 1) + i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}
				l = (i + 2) % 3;
				for (i1 = 0; i1<elementDOF[1].dop - 1; i1++)
				{
					node = elementDOF[1].val[k][3 + l*(elementDOF[1].dop - 1) + i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}

				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(C, count, element[0], edge, elementEdge, &elementDOF[1], rowstart, row21);

				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(C, count, element[0], edge, elementEdge, &elementDOF[1], rowstart, row21);
			} // if(i<3)
			else
			{
				l = (i - 3) / (elementDOF[1].dop - 1);
				node = elementDOF[1].val[k][(l + 1) % 3];
				for (j = 0; j<2; j++)
				{
					C->JA[rowstart[j] + count] = node;
					C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
				}
				count++;
				node = elementDOF[1].val[k][(l + 2) % 3];
				for (j = 0; j<2; j++)
				{
					C->JA[rowstart[j] + count] = node;
					C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
				}
				count++;
				for (i1 = 0; i1<elementDOF[1].dop - 1; i1++)
				{
					node = elementDOF[1].val[k][3 + l*(elementDOF[1].dop - 1) + i1];
					for (j = 0; j<2; j++)
					{
						C->JA[rowstart[j] + count] = node;
						C->JA[rowstart[j] + count + row21[j]] = node + elementDOF[1].dof;
					}
					count++;
				}

				edge = elementEdge->val[k][l];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(C, count, element[0], edge, elementEdge, &elementDOF[1], rowstart, row21);
			}  // if(i<3) else
		} // i
	} // k

	  // step 3C: Loop edge by edge and compute the actual entries storing them in C
	double elen, C11;
	int patchnodes[200];
	C->val = (double*)calloc(C->nnz, sizeof(double));
	for (edge = 0; edge<edges->row; edge++)
	{
		//		edgeNode[0]=edges->val[edge][0];
		//		edgeNode[1]=edges->val[edge][1];
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

		//		C11 = elen;

		count = 0;
		if (elementDOF[1].dop == 0)
		{
			patchnodes[count] = elementDOF[1].val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF[1].val[element[1]][0];
				count++;
			}
		} // if(elementDOF[1].dop==0)
		else
		{
			for (i = 0; i<3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 1) % 3];
			count++;
			patchnodes[count] = elementDOF[1].val[element[0]][(i + 2) % 3];
			count++;
			for (j = 0; j<elementDOF[1].dop - 1; j++)
			{
				patchnodes[count] = elementDOF[1].val[element[0]][3 + i*(elementDOF[1].dop - 1) + j];
				count++;
			}

			if (element[1] != -1)
			{
				for (i = 0; i<3; i++)
				{
					if (elementEdge->val[element[1]][i] == edge)
						break;
				}
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 1) % 3];
				count++;
				patchnodes[count] = elementDOF[1].val[element[1]][(i + 2) % 3];
				count++;
				for (j = 0; j<elementDOF[1].dop - 1; j++)
				{
					patchnodes[count] = elementDOF[1].val[element[1]][3 + i*(elementDOF[1].dop - 1) + j];
					count++;
				}
			}
		} // if(elementDOF[1].dop==0) else

		for (k1 = 0; k1<count; k1++)
		{
			curnode[0] = patchnodes[k1];
			curnode[1] = curnode[0] + elementDOF[1].dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = C->IA[curnode[j]];
				row21[j] = (C->IA[curnode[j] + 1] - C->IA[curnode[j]]) / 2;
			}
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				// c11
				for (j1 = C->IA[curnode[0]]; j1<C->IA[curnode[0]] + row21[0]; j1++)
				{
					if (C->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], curnode[0], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], j, phi2);
							C->val[j1] += gauss11[i1][1] * phi1[0] * phi2[0];
						}
						break;
					}
				}
				// c22
				for (j1 = C->IA[curnode[1]] + row21[1]; j1<C->IA[curnode[1] + 1]; j1++)
				{
					if (C->JA[j1] == (j + elementDOF[1].dof))
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], curnode[1], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], j + elementDOF[1].dof, phi2);
							C->val[j1] += gauss11[i1][1] * phi1[1] * phi2[1];
						}
						break;
					}
				}
			} // k2
		} // k1
	} // edge



	  /************************************************** right hand side b *****************************************************************/
	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	int num_qp10 = 20; // the number of numerical intergation points
	double gauss10[num_qp10][2];
	init_Gauss1D(num_qp10, 1, gauss10); // gauss intergation initial

	for (edge = 0; edge<edges->row; edge++)
	{
		if (edges->bdFlag[edge] <1 || edges->bdFlag[edge]>4) //////////////////////////////////////////
			continue;

		// only consider Dirichlet boundary
		k = edges->val[edge][2];

		// set parameters
		elen = edges->length[edge];
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		for (l = 0; l<3; l++)
		{
			if (elementEdge->val[k][l] == edge)
				break;
		}

		l1 = (l + 1) % 3;
		if (elements->val[k][l1] != edges->val[edge][0])
			l1 = (l + 2) % 3;
		l2 = 3 - l - l1;

		count = getlocaldofsEdge4HuZhang2d(patchnodes, l, elementDOF[0].dop);

		for (i = 0; i<count; i++)
		{
			k1 = patchnodes[i];
			for (i1 = 0; i1<num_qp10; i1++)
			{
				lambdas[l1] = gauss10[i1][0];
				lambdas[l2] = 1 - gauss10[i1][0];
				lambdas[l] = 0;
				huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b1->val[elementDOF[0].val[k][k1]] += elen*gauss10[i1][1] * (phi1[0] * elements->nvector[k][l][0] + phi1[2] * elements->nvector[k][l][1])*u1(x, y, lambda, mu);
				b1->val[elementDOF[0].val[k][k1]] += elen*gauss10[i1][1] * (phi1[2] * elements->nvector[k][l][0] + phi1[1] * elements->nvector[k][l][1])*u2(x, y, lambda, mu);
			}
		}
	} // edge

	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF[1].col; k1++)
		{
			i = elementDOF[1].val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b2->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
				b2->val[i + elementDOF[1].dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k
}

/**
 * \fn void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
 * \brief assemble stiffness matrix 
 * \param *A pointer to stiffness matrix
 * \param *b pointer to right hand side
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *elementdofTran pointer to transpose of elementDOF
 * \param mu Lame constant or Poisson ratio of plate
 * \return void
 */
void assembleStiffmatrixElasLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i,j,k,l;

	A->row = elementDOF->dof * 2;
	A->col = A->row;
	A->IA=(int*)calloc(A->row+1, sizeof(int));
	A->JA=NULL;
	A->val=NULL;
	
//	create_dvector(A->row, b);

	int nvertices=nodes->row;
	int nedges=edges->row;
	int element[2], edge, node;
	
	double phi, phi1[2], phi2[2];
	int k1,k2,i1,j1,l1,l2,ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;
	
	int num_qp21=getNumQuadPoints(elementDOF->dop*2-2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial
	
	/************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index=(int*)calloc(A->col, sizeof(int));
	for(i=0;i<A->col;i++)
		index[i]=-1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for(i=0;i<elementDOF->dof;i++)
	{
		count=0;
		istart=-2;
		for(j=elementdofTran->IA[i];j<elementdofTran->IA[i+1];j++)
		{
			element[0]=elementdofTran->JA[j];

			for(k=0;k<elementDOF->col;k++)
			{
				node=elementDOF->val[element[0]][k];
				if(index[node]==-1)
				{
					index[node]=istart;
					istart=node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[elementDOF->dof + i + 1] = count * 2;

		for(j=0;j<count;j++)
		{
			l=istart;
			istart=index[l];
			index[l]=-1;
		}		
	} // i
	
	for(i=0;i<A->row;i++)
		A->IA[i+1]+=A->IA[i];
	
	A->nnz=A->IA[A->row];
	
	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA=(int*)calloc(A->nnz,sizeof(int));
	for(i=0;i<elementDOF->dof;i++)
	{
		istart=-2;
		for(j=elementdofTran->IA[i];j<elementdofTran->IA[i+1];j++)
		{
			element[0]=elementdofTran->JA[j];

			for(k=0;k<elementDOF->col;k++)
			{

				node=elementDOF->val[element[0]][k];
				if(index[node]==-1)
				{
					index[node]=istart;
					istart=node;
				}
			}
		}

		row21[0] = (A->IA[i + 1] - A->IA[i]) / 2;
		rowstart[1] = A->IA[i + elementDOF->dof];
		for (j = A->IA[i]; j<A->IA[i] + row21[0]; j++)
		{
			A->JA[j] = istart;
			A->JA[j + row21[0]] = istart + elementDOF->dof;
			A->JA[j - A->IA[i] + rowstart[1]] = A->JA[j];
			A->JA[j - A->IA[i] + rowstart[1] + row21[0]] = A->JA[j + row21[0]];
			istart=index[istart];
			index[A->JA[j]]=-1;
		}
	} // i
	free(index);
	
	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val=(double*)calloc(A->nnz, sizeof(double));
	for(k=0;k<elements->row;k++)
	{
		// set parameters
		s=elements->vol[k];
		xi=elements->xi[k];
		eta=elements->eta[k];
		// end set parameters

		for(k1=0;k1<elementDOF->col;k1++)
		{
			curnode[0] = elementDOF->val[k][k1];
			curnode[1] = curnode[0] + elementDOF->dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for(k2=0;k2<elementDOF->col;k2++)
			{
				j=elementDOF->val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if(A->JA[j1]==j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2) * 2 * mu;
						}
						break;
					}
				} // j1
				// a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF->dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[1] * phi2[0]/ 2 * 2 * mu;
						}
						break;
					}
				} // j1
				// a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] / 2 * 2 * mu;
						}
						break;
					}
				} // j1
				// a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	
	/************************************************** right hand side b *****************************************************************/
/*	int num_qp0=49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for(k=0;k<elements->row;k++)
	{
		for(i=0;i<3;i++)
		{
			j=elements->val[k][i];
			xs[i]=nodes->val[j][0];
			ys[i]=nodes->val[j][1];
		}
		s=elements->vol[k];
		
		for(k1=0;k1<elementDOF->col;k1++)
		{
			i=elementDOF->val[k][k1];
			for(i1=0;i1<num_qp0;i1++)
			{
				lambdas[0]=gauss0[i1][0];
				lambdas[1]=gauss0[i1][1];
				lambdas[2]=1-lambdas[0]-lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
				x=xs[0]*lambdas[0]+xs[1]*lambdas[1]+xs[2]*lambdas[2];
				y=ys[0]*lambdas[0]+ys[1]*lambdas[1]+ys[2]*lambdas[2];
				b->val[i]+=2*s*gauss0[i1][2]*f1(x,y, lambda, mu)*phi;
				b->val[i+elementDOF->dof]+=2*s*gauss0[i1][2]*f2(x,y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k  */
}

/**
* \fn void assembleStiffmatrixVecLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleStiffmatrixVecLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i, j, k, l;

	A->row = elementDOF->dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF->dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count * 2;
		A->IA[elementDOF->dof + i + 1] = count * 2;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF->dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{

				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		row21[0] = (A->IA[i + 1] - A->IA[i]) / 2;
		rowstart[1] = A->IA[i + elementDOF->dof];
		for (j = A->IA[i]; j<A->IA[i] + row21[0]; j++)
		{
			A->JA[j] = istart;
			A->JA[j + row21[0]] = istart + elementDOF->dof;
			A->JA[j - A->IA[i] + rowstart[1]] = A->JA[j];
			A->JA[j - A->IA[i] + rowstart[1] + row21[0]] = A->JA[j + row21[0]];
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			curnode[0] = elementDOF->val[k][k1];
			curnode[1] = curnode[0] + elementDOF->dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0]  + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k


	  /************************************************** right hand side b *****************************************************************/
	  /*	int num_qp0=49; // the number of numerical intergation points
	  double gauss0[num_qp0][3];
	  init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	  for(k=0;k<elements->row;k++)
	  {
	  for(i=0;i<3;i++)
	  {
	  j=elements->val[k][i];
	  xs[i]=nodes->val[j][0];
	  ys[i]=nodes->val[j][1];
	  }
	  s=elements->vol[k];

	  for(k1=0;k1<elementDOF->col;k1++)
	  {
	  i=elementDOF->val[k][k1];
	  for(i1=0;i1<num_qp0;i1++)
	  {
	  lambdas[0]=gauss0[i1][0];
	  lambdas[1]=gauss0[i1][1];
	  lambdas[2]=1-lambdas[0]-lambdas[1];
	  lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
	  x=xs[0]*lambdas[0]+xs[1]*lambdas[1]+xs[2]*lambdas[2];
	  y=ys[0]*lambdas[0]+ys[1]*lambdas[1]+ys[2]*lambdas[2];
	  b->val[i]+=2*s*gauss0[i1][2]*f1(x,y, lambda, mu)*phi;
	  b->val[i+elementDOF->dof]+=2*s*gauss0[i1][2]*f2(x,y, lambda, mu)*phi;
	  } // i1
	  } // k1
	  } // k  */
}

/**
* \fn void assembleStiffmatrixLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param *elementdofTran pointer to transpose of elementDOF
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void assembleStiffmatrixLagrange(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, iCSRmat *elementdofTran, double mu)
{
	int i, j, k, l;

	A->row = elementDOF->dof;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

									  /************************************************** stiffness matrix A *****************************************************************/
	int *index;
	int istart;
	index = (int*)calloc(A->col, sizeof(int));
	for (i = 0; i<A->col; i++)
		index[i] = -1;
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (i = 0; i<elementDOF->dof; i++)
	{
		count = 0;
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{
				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
					count++;
				}
			}
		}
		A->IA[i + 1] = count;

		for (j = 0; j<count; j++)
		{
			l = istart;
			istart = index[l];
			index[l] = -1;
		}
	} // i

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	for (i = 0; i<elementDOF->dof; i++)
	{
		istart = -2;
		for (j = elementdofTran->IA[i]; j<elementdofTran->IA[i + 1]; j++)
		{
			element[0] = elementdofTran->JA[j];

			for (k = 0; k<elementDOF->col; k++)
			{

				node = elementDOF->val[element[0]][k];
				if (index[node] == -1)
				{
					index[node] = istart;
					istart = node;
				}
			}
		}

		for (j = A->IA[i]; j<A->IA[i + 1]; j++)
		{
			A->JA[j] = istart;
			istart = index[istart];
			index[A->JA[j]] = -1;
		}
	} // i
	free(index);

	// step 3A: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				for (j1 = A->IA[i]; j1<A->IA[i+1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]) * 2 * mu;
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k


	  /************************************************** right hand side b *****************************************************************/
	  /*	int num_qp0=49; // the number of numerical intergation points
	  double gauss0[num_qp0][3];
	  init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	  for(k=0;k<elements->row;k++)
	  {
	  for(i=0;i<3;i++)
	  {
	  j=elements->val[k][i];
	  xs[i]=nodes->val[j][0];
	  ys[i]=nodes->val[j][1];
	  }
	  s=elements->vol[k];

	  for(k1=0;k1<elementDOF->col;k1++)
	  {
	  i=elementDOF->val[k][k1];
	  for(i1=0;i1<num_qp0;i1++)
	  {
	  lambdas[0]=gauss0[i1][0];
	  lambdas[1]=gauss0[i1][1];
	  lambdas[2]=1-lambdas[0]-lambdas[1];
	  lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
	  x=xs[0]*lambdas[0]+xs[1]*lambdas[1]+xs[2]*lambdas[2];
	  y=ys[0]*lambdas[0]+ys[1]*lambdas[1]+ys[2]*lambdas[2];
	  b->val[i]+=2*s*gauss0[i1][2]*f(x,y)*phi;
	  } // i1
	  } // k1
	  } // k  */
}


/**
* \fn void assembleStiffmatrixElasIPDG0(dCSRmat *A, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assembleStiffmatrixElasIPDG0(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int dop = elementDOF->dop;
	int dof = elementDOF->dof;
	int i, j, k, l;

	A->row = elementDOF->dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF->dop * 2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial


	/************************************************** stiffness matrix A *****************************************************************/
	
	// step 1A: Find first the structure IA of the stiffness matrix A
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2;
			A->IA[curnode[1] + 1] += 2;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2;
					A->IA[curnode[1] + 1] += 2;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2 * elementDOF->col;
			A->IA[curnode[1] + 1] += 2 * elementDOF->col;
			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
		} // i
	} // k

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
//	count = 0;
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (j = 0; j<2; j++)
			{
				A->JA[rowstart[j] + count] = curnode[0];
				A->JA[rowstart[j] + count + row21[j]] = curnode[1];
			}
			count++;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					for (j = 0; j<2; j++)
					{
						A->JA[rowstart[j] + count] = elementDOF->val[element[0]][0];
						A->JA[rowstart[j] + count + row21[j]] = elementDOF->val[element[0]][0] + dof;
					}
					count++;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (i1 = 0; i1<elementDOF->col; i1++)
			{
				node = elementDOF->val[k][i1];
				for (j = 0; j<2; j++)
				{
					A->JA[rowstart[j] + count] = node;
					A->JA[rowstart[j] + count + row21[j]] = node + dof;
				}
				count++;
			}

			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);

				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
		} // i
	} // k
	
	// step 3A1: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
			continue;

		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			curnode[0] = elementDOF->val[k][k1];
			curnode[1] = curnode[0] + elementDOF->dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2);
						}
						break;
					}
				} // j1
				 // a12
				for (j1 = A->IA[curnode[0]] + row21[0]; j1<A->IA[curnode[0] + 1]; j1++)
				{
					if (A->JA[j1] == j + elementDOF->dof)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[1] * phi2[0] / 2;
						}
						break;
					}
				} // j1
				  // a21
				for (j1 = A->IA[curnode[1]]; j1<A->IA[curnode[1]] + row21[1]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * phi1[0] * phi2[1] / 2;
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	 // step 3A2: Loop edge by edge and compute the actual entries storing them in A
	double elen;
	int patchnodes[100];
	for (edge = 0; edge<edges->row; edge++)
	{
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

		count = 0;
		if (dop == 0)
		{
			patchnodes[count] = elementDOF->val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF->val[element[1]][0];
				count++;
			}
		} // if(dop==0)
		else
		{
			for (i = 0; i<3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			patchnodes[count] = elementDOF->val[element[0]][(i + 1) % 3];
			count++;
			patchnodes[count] = elementDOF->val[element[0]][(i + 2) % 3];
			count++;
			for (j = 0; j<dop - 1; j++)
			{
				patchnodes[count] = elementDOF->val[element[0]][3 + i*(dop - 1) + j];
				count++;
			}

			if (element[1] != -1)
			{
				for (i = 0; i<3; i++)
				{
					if (elementEdge->val[element[1]][i] == edge)
						break;
				}
				patchnodes[count] = elementDOF->val[element[1]][(i + 1) % 3];
				count++;
				patchnodes[count] = elementDOF->val[element[1]][(i + 2) % 3];
				count++;
				for (j = 0; j<dop - 1; j++)
				{
					patchnodes[count] = elementDOF->val[element[1]][3 + i*(dop - 1) + j];
					count++;
				}
			}
		} // if(dop==0) else

		for (k1 = 0; k1<count; k1++)
		{
			curnode[0] = patchnodes[k1];
			curnode[1] = curnode[0] + elementDOF->dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[0], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j, phi2);
							A->val[j1] += gauss11[i1][1]*phi1[0] * phi2[0];
						}
						break;
					}
				}
				// a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[1], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j + elementDOF->dof, phi2);
							A->val[j1] += gauss11[i1][1]*phi1[1] * phi2[1];
						}
						break;
					}
				}
			} // k2
		} // k1
	} // edge
	
	  /************************************************** right hand side b *****************************************************************/
/*	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k<elements->row; k++)
	{
		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}
		s = elements->vol[k];

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			i = elementDOF->val[k][k1];
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
				b->val[i + elementDOF->dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
			} // i1
		} // k1 
	} // k  */
}

/**
* \fn void assembleStiffmatrixVecIPDG0(dCSRmat *A, dvector *b, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief assemble stiffness matrix
* \param *A pointer to stiffness matrix
* \param *b pointer to right hand side
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void assembleStiffmatrixVecIPDG0(dCSRmat *A, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int dop = elementDOF->dop;
	int dof = elementDOF->dof;
	int i, j, k, l;

	A->row = elementDOF->dof * 2;
	A->col = A->row;
	A->IA = (int*)calloc(A->row + 1, sizeof(int));
	A->JA = NULL;
	A->val = NULL;

	//	create_dvector(A->row, b);

	int nvertices = nodes->row;
	int nedges = edges->row;
	int element[2], edge, node;

	double phi, phi1[2], phi2[2];
	int k1, k2, i1, j1, l1, l2, ej;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi;
	int rowstart[2], row21[2], curnode[2];
	int count;

	int num_qp21 = getNumQuadPoints(elementDOF->dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss21[num_qp21][3];
	init_Gauss(num_qp21, 2, gauss21); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF->dop * 2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial


										/************************************************** stiffness matrix A *****************************************************************/

										// step 1A: Find first the structure IA of the stiffness matrix A
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2;
			A->IA[curnode[1] + 1] += 2;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2;
					A->IA[curnode[1] + 1] += 2;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			A->IA[curnode[0] + 1] += 2 * elementDOF->col;
			A->IA[curnode[1] + 1] += 2 * elementDOF->col;
			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					A->IA[curnode[0] + 1] += 2 * (dop + 1);
					A->IA[curnode[1] + 1] += 2 * (dop + 1);
				}
			}
		} // i
	} // k

	for (i = 0; i<A->row; i++)
		A->IA[i + 1] += A->IA[i];

	A->nnz = A->IA[A->row];

	// step 2A: Find the structure JA of the stiffness matrix A
	A->JA = (int*)calloc(A->nnz, sizeof(int));
	//	count = 0;
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
		{
			curnode[0] = elementDOF->val[k][0];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (j = 0; j<2; j++)
			{
				A->JA[rowstart[j] + count] = curnode[0];
				A->JA[rowstart[j] + count + row21[j]] = curnode[1];
			}
			count++;
			for (i = 0; i<3; i++)
			{
				edge = elementEdge->val[k][i];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
				{
					for (j = 0; j<2; j++)
					{
						A->JA[rowstart[j] + count] = elementDOF->val[element[0]][0];
						A->JA[rowstart[j] + count + row21[j]] = elementDOF->val[element[0]][0] + dof;
					}
					count++;
				}
			}
			continue;
		}

		for (i = 0; i<elementDOF->col; i++) //  for each node
		{
			curnode[0] = elementDOF->val[k][i];
			curnode[1] = curnode[0] + dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			count = 0;
			for (i1 = 0; i1<elementDOF->col; i1++)
			{
				node = elementDOF->val[k][i1];
				for (j = 0; j<2; j++)
				{
					A->JA[rowstart[j] + count] = node;
					A->JA[rowstart[j] + count + row21[j]] = node + dof;
				}
				count++;
			}

			if (i<3)
			{
				edge = elementEdge->val[k][(i + 1) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);

				edge = elementEdge->val[k][(i + 2) % 3];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
			else if (i<3 * dop)
			{
				edge = elementEdge->val[k][(i - 3) / (dop - 1)];
				element[0] = edges->val[edge][2] + edges->val[edge][3] - k;
				if (element[0] != -1)
					count = getEdgeDOFsVector(A, count, element[0], edge, elementEdge, elementDOF, rowstart, row21);
			}
		} // i
	} // k

	  // step 3A1: Loop element by element and compute the actual entries storing them in A
	A->val = (double*)calloc(A->nnz, sizeof(double));
	for (k = 0; k<elements->row; k++)
	{
		if (dop == 0)
			continue;

		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		// end set parameters

		for (k1 = 0; k1<elementDOF->col; k1++)
		{
			curnode[0] = elementDOF->val[k][k1];
			curnode[1] = curnode[0] + elementDOF->dof;

			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}

			for (k2 = 0; k2<elementDOF->col; k2++)
			{
				j = elementDOF->val[k][k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
				  // a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp21; i1++)
						{
							lambdas[0] = gauss21[i1][0];
							lambdas[1] = gauss21[i1][1];
							lambdas[2] = 1 - lambdas[0] - lambdas[1];
							lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF->dop, phi1);
							lagrange_basis1(lambdas, s, eta, xi, k2, elementDOF->dop, phi2);
							A->val[j1] += 2 * s*gauss21[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1]);
						}
						break;
					}
				} // j1
			} // k2
		} // k1
	} // k

	  // step 3A2: Loop edge by edge and compute the actual entries storing them in A
	double elen;
	int patchnodes[100];
	for (edge = 0; edge<edges->row; edge++)
	{
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];

		count = 0;
		if (dop == 0)
		{
			patchnodes[count] = elementDOF->val[element[0]][0];
			count++;
			if (element[1] != -1)
			{
				patchnodes[count] = elementDOF->val[element[1]][0];
				count++;
			}
		} // if(dop==0)
		else
		{
			for (i = 0; i<3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			patchnodes[count] = elementDOF->val[element[0]][(i + 1) % 3];
			count++;
			patchnodes[count] = elementDOF->val[element[0]][(i + 2) % 3];
			count++;
			for (j = 0; j<dop - 1; j++)
			{
				patchnodes[count] = elementDOF->val[element[0]][3 + i*(dop - 1) + j];
				count++;
			}

			if (element[1] != -1)
			{
				for (i = 0; i<3; i++)
				{
					if (elementEdge->val[element[1]][i] == edge)
						break;
				}
				patchnodes[count] = elementDOF->val[element[1]][(i + 1) % 3];
				count++;
				patchnodes[count] = elementDOF->val[element[1]][(i + 2) % 3];
				count++;
				for (j = 0; j<dop - 1; j++)
				{
					patchnodes[count] = elementDOF->val[element[1]][3 + i*(dop - 1) + j];
					count++;
				}
			}
		} // if(dop==0) else

		for (k1 = 0; k1<count; k1++)
		{
			curnode[0] = patchnodes[k1];
			curnode[1] = curnode[0] + elementDOF->dof;
			for (j = 0; j<2; j++)
			{
				rowstart[j] = A->IA[curnode[j]];
				row21[j] = (A->IA[curnode[j] + 1] - A->IA[curnode[j]]) / 2;
			}
			for (k2 = 0; k2<count; k2++)
			{
				j = patchnodes[k2];
				// a11
				for (j1 = A->IA[curnode[0]]; j1<A->IA[curnode[0]] + row21[0]; j1++)
				{
					if (A->JA[j1] == j)
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[0], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j, phi2);
							A->val[j1] += gauss11[i1][1] * phi1[0] * phi2[0];
						}
						break;
					}
				}
				// a22
				for (j1 = A->IA[curnode[1]] + row21[1]; j1<A->IA[curnode[1] + 1]; j1++)
				{
					if (A->JA[j1] == (j + elementDOF->dof))
					{
						for (i1 = 0; i1<num_qp11; i1++)
						{
							lambdas[0] = gauss11[i1][0];
							lambdas[1] = 1 - lambdas[0];
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, curnode[1], phi1);
							jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, elementDOF, j + elementDOF->dof, phi2);
							A->val[j1] += gauss11[i1][1] * phi1[1] * phi2[1];
						}
						break;
					}
				}
			} // k2
		} // k1
	} // edge

	  /************************************************** right hand side b *****************************************************************/
	  /*	int num_qp0 = 49; // the number of numerical intergation points
	  double gauss0[num_qp0][3];
	  init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	  for (k = 0; k<elements->row; k++)
	  {
	  for (i = 0; i<3; i++)
	  {
	  j = elements->val[k][i];
	  xs[i] = nodes->val[j][0];
	  ys[i] = nodes->val[j][1];
	  }
	  s = elements->vol[k];

	  for (k1 = 0; k1<elementDOF->col; k1++)
	  {
	  i = elementDOF->val[k][k1];
	  for (i1 = 0; i1<num_qp0; i1++)
	  {
	  lambdas[0] = gauss0[i1][0];
	  lambdas[1] = gauss0[i1][1];
	  lambdas[2] = 1 - lambdas[0] - lambdas[1];
	  lagrange_basis(lambdas, k1, elementDOF->dop, &phi);
	  x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
	  y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
	  b->val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi;
	  b->val[i + elementDOF->dof] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi;
	  } // i1
	  } // k1
	  } // k  */
}

/**
 * \fnvoid sumNormalDerivativeEta(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, ddenmat *etas, double *sum)
 * \brief the sum operator for normal derivative of basis function with penalty parameter
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *etas pointer to penalty parameter 
 * \param *sum pointer to the result of sum operator
 * \return void
 */
void sumNormalDerivativeEta(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, ddenmat *etas, double *sum)
{
	int i,j,l;
	double lambdas[3];
	double s, *eta, *xi, phi[2], nve[2];
	int nodeindex1,nodeindex2;
	int element[2], edgeNode[2];
	int dop=elementDOF->dop;

	*sum=0;
	if(dop<=0)
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	element[0]=edges->val[edge][2];
	element[1]=edges->val[edge][3];
	nve[0]=edges->nvector[edge][0];
	nve[1]=edges->nvector[edge][1];
	
	for(nodeindex1=0;nodeindex1<elementDOF->col;nodeindex1++)
    {
		if(elementDOF->val[element[0]][nodeindex1]==node)
			break;
	}

	if(nodeindex1<elementDOF->col)
	{
		for(i=0;i<3;i++)
		{
			if(elements->val[element[0]][i]==edgeNode[0])
				break;
		}
	
		for(j=0;j<3;j++)
		{
			if(elements->val[element[0]][j]==edgeNode[1])
				break;
		}
	
		l=3-i-j;
		lambdas[i]=lambda1;
		lambdas[j]=lambda2;
		lambdas[l]=0;
		s=elements->vol[element[0]];
		xi=elements->xi[element[0]];
		eta=elements->eta[element[0]];
		lagrange_basis1(lambdas, s, eta, xi, nodeindex1, dop, phi);
		*sum+=(phi[0]*nve[0]+phi[1]*nve[1])*etas->val[element[0]][l];
	}

	if(element[1]==-1)
			return;

	for(nodeindex2=0;nodeindex2<elementDOF->col;nodeindex2++)
    {
		if(elementDOF->val[element[1]][nodeindex2]==node)
			break;
	}

	if(nodeindex2<elementDOF->col)
	{
		for(i=0;i<3;i++)
		{
			if(elements->val[element[1]][i]==edgeNode[0])
				break;
		}
	
		for(j=0;j<3;j++)
		{
			if(elements->val[element[1]][j]==edgeNode[1])
				break;
		}
	
		l=3-i-j;
		lambdas[i]=lambda1;
		lambdas[j]=lambda2;
		lambdas[l]=0;
		s=elements->vol[element[1]];
		xi=elements->xi[element[1]];
		eta=elements->eta[element[1]];
		lagrange_basis1(lambdas, s, eta, xi, nodeindex2, dop, phi);
		*sum+=(phi[0]*nve[0]+phi[1]*nve[1])*etas->val[element[1]][l];
	}
}

/**
 * \fn void jumpOperatorVectorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, DOF *dofs, int node, double *jump)
 * \brief the jump operator for vector otimes normal derivative
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorVectorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex, ni2;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2];
	double phi;
	double elen=edges->length[edge];

	jump[0]=0;
	jump[1]=0;
	jump[2]=0;

	if(node<0 && node>=elementDOF->dof*2)
		return;

	element=node/(elementDOF->col*2);
	ni2=node%(elementDOF->col*2);

	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;

	if(elementDOF->dop==0)
	{
		if(ni2<elementDOF->col)
		{
			jump[0]=nv[0];
			jump[1]=0;
			jump[2]=nv[1]/2;
		}
		else
		{
			jump[0]=0;
			jump[1]=nv[1];
			jump[2]=nv[0]/2;
		}
		return;
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);


	if(ni2<elementDOF->col)
	{
		jump[0]=phi*nv[0];
		jump[1]=0;
		jump[2]=phi*nv[1]/2;
	}
	else
	{
		jump[0]=0;
		jump[1]=phi*nv[1];
		jump[2]=phi*nv[0]/2;
	}
}

/**
* \fn void jumpOperatorVector(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, DOF *dofs, int node, double *jump)
* \brief the jump operator for vector
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorVector(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i, j, l;
	int nodeindex, edgeindex, ni2;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], nve[2];
	double phi;
	double elen = edges->length[edge];

	jump[0] = 0;
	jump[1] = 0;

	if (node<0 && node >= elementDOF->dof * 2)
		return;

	element = (node % elementDOF->dof) / elementDOF->col;
	int li = node / elementDOF->dof;

	for (edgeindex = 0; edgeindex<3; edgeindex++)
	{
		if (elementEdge->val[element][edgeindex] == edge)
			break;
	}

	if (edgeindex == 3)
		return;

	xi = elements->xi[element];
	eta = elements->eta[element];
	nv[0] = -eta[edgeindex] / elen;
	nv[1] = xi[edgeindex] / elen;
	nve[0] = edges->nvector[edge][0];
	nve[1] = edges->nvector[edge][1];
	double sgn = nve[0] * nv[0] + nve[1] * nv[1];

	if (elementDOF->dop == 0)
	{
		jump[li] = sgn;
		return;
	}

	nodeindex = node%elementDOF->col;

	if (nodeindex<3 && nodeindex >= 0)
	{
		if (nodeindex == edgeindex)
			return;
	}
	else if (nodeindex<3 * elementDOF->dop && nodeindex >= 3)
	{
		if ((nodeindex - 3) / (elementDOF->dop - 1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
	for (i = 0; i<3; i++)
	{
		if (elements->val[element][i] == edgeNode[0])
			break;
	}
	for (j = 0; j<3; j++)
	{
		if (elements->val[element][j] == edgeNode[1])
			break;
	}
	l = 3 - i - j;
	lambdas[i] = lambda1;
	lambdas[j] = lambda2;
	lambdas[l] = 0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	jump[li] = phi*sgn;
}


/**
 * \fn void jumpOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for normal derivative of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i;
	int nodeindex1,nodeindex2, edgeindex;
	int element[2], edgeNode[2];
	double jumpminus;
	int isInedge=0;
	int nedges=edges->row;
	int dop=elementDOF->dop;

	*jump=0;
	if(dop<=0)
		return;

	double elen=edges->length[edge];
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	element[0]=edges->val[edge][2];
	element[1]=edges->val[edge][3];
	
	for(nodeindex1=0;nodeindex1<elementDOF->col;nodeindex1++)
    {
		if(elementDOF->val[element[0]][nodeindex1]==node)
			break;
	}
	
	if(element[1]==-1)
	{
		if(nodeindex1==elementDOF->col)
			return;
	}
	else
	{
		for(nodeindex2=0;nodeindex2<elementDOF->col;nodeindex2++)
		{
			if(elementDOF->val[element[1]][nodeindex2]==node)
				break;
		}
		
		if(nodeindex1==elementDOF->col && nodeindex2==elementDOF->col)
			return;
	}
	
	for(edgeindex=0;edgeindex<elementEdge->col;edgeindex++)
	{
		if(elementEdge->val[element[0]][edgeindex]==edge)
			break;
	}
	
	for(i=0;i<2;i++)
	{
		if(edgeNode[i]==node)
		{
			isInedge=1;
			break;
		}
	}
	
	if(isInedge==0)
	{
		for(i=0;i<dop-1;i++)
		{
			if(elementDOF->val[element[0]][3+edgeindex*(dop-1)+i]==node)
			{
				isInedge=1;
				break;
			}
		}
	}
	
	if(isInedge==1)
	{
		getinfo4jump(lambda1, lambda2, edgeNode, elen, element[0], elements, dop, nodeindex1, jump);
		if(element[1]!=-1)
		{			
			getinfo4jump(lambda1, lambda2, edgeNode, elen, element[1], elements, dop, nodeindex2, &jumpminus);
				*jump += jumpminus;
		}
	}
	else
	{
		if(nodeindex1<elementDOF->col)
			getinfo4jump(lambda1, lambda2, edgeNode, elen, element[0], elements, dop, nodeindex1, jump);
		else{
				if(element[1]!=-1)
				{
					if(nodeindex2<elementDOF->col)
						getinfo4jump(lambda1, lambda2, edgeNode, elen, element[1], elements, dop, nodeindex2, jump);
				}
			}
	}
}

/**
 * \fn void getinfo4jump(double lambda1, double lambda2, int *edgeNode, double elen, int element, ELEMENT *elements, int dop, int index, double *jump)
 * \brief the normal derivative of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param *edgeNode pointer to the two indices of current edge
 * \param elen the length of cuurent edge
 * \param element the index of current element
 * \param *elements pointer to the structure of the triangulation
 * \param dop number of degrees of polynomial
 * \param index index of current node variable
 * \param *jump pointer to the normal tensor of gradient of basis function
 * \return void
 */
void getinfo4jump(double lambda1, double lambda2, int *edgeNode, double elen, int element, ELEMENT *elements, int dop, int index, double *jump)
{
	int i,j,l;
	double lambdas[3];
	double s, *eta, *xi, phi[2], nv[2];
		
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;
	
	s=elements->vol[element];
	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[l]/elen;
	nv[1]=xi[l]/elen;
	
	lagrange_basis1(lambdas, s, eta, xi, index, dop, phi);
		
	*jump=phi[0]*nv[0]+phi[1]*nv[1];
}

/**
 * \fn void TangentDerivative4Edge(double lambda1, double lambda2, int edge, int element, ELEMENT *elements, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *val)
 * \brief the tangential derivative of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *edgeNode pointer to the two indices of current edge
 * \param element the index of current element
 * \param *elements pointer to the structure of the triangulation
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *val pointer to the tangential derivative of basis function on edge
 * \return void
 */
void TangentDerivative4Edge(double lambda1, double lambda2, int edge, int element, ELEMENT *elements, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *val)
{
	int i,j,l;
	double lambdas[3];
	double s, *eta, *xi, phi[2], tve[2];
	int nodeindex, edgeNode[2];
	int dop=elementDOF->dop;

	*val=0;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];

	for(nodeindex=0;nodeindex<elementDOF->col;nodeindex++)
    {
		if(elementDOF->val[element][nodeindex]==node)
			break;
	}

	if(nodeindex==elementDOF->col)
			return;
		
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;
	
	s=elements->vol[element];
	xi=elements->xi[element];
	eta=elements->eta[element];
	tve[0]=edges->tvector[edge][0];
	tve[1]=edges->tvector[edge][1];
	
	lagrange_basis1(lambdas, s, eta, xi, nodeindex, dop, phi);
		
	*val=phi[0]*tve[0]+phi[1]*tve[1];
}

/**
 * \fn void averageOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *average)
 * \brief the average operator for gradient of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void averageOperator(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *average)
{
	int i;
	int nodeindex1,nodeindex2, edgeindex;
	int element[2], edgeNode[2];
	double averageplus[2];
	int isInedge=0;
	int nedges=edges->row;
	int dop=elementDOF->dop;

	average[0]=0;
	average[1]=0;
	if(dop<=0)
		return;

	double elen=edges->length[edge];
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	element[0]=edges->val[edge][2];
	element[1]=edges->val[edge][3];
	
	for(nodeindex1=0;nodeindex1<elementDOF->col;nodeindex1++)
    {
		if(elementDOF->val[element[0]][nodeindex1]==node)
			break;
	}
	
	if(element[1]==-1)
	{
		if(nodeindex1==elementDOF->col)
			return;
	}
	else
	{
		for(nodeindex2=0;nodeindex2<elementDOF->col;nodeindex2++)
		{
			if(elementDOF->val[element[1]][nodeindex2]==node)
				break;
		}
		
		if(nodeindex1==elementDOF->col && nodeindex2==elementDOF->col)
			return;
	}
	
	for(edgeindex=0;edgeindex<elementEdge->col;edgeindex++)
	{
		if(elementEdge->val[element[0]][edgeindex]==edge)
			break;
	}
	
	for(i=0;i<2;i++)
	{
		if(edgeNode[i]==node)
		{
			isInedge=1;
			break;
		}
	}
	
	if(isInedge==0)
	{
		for(i=0;i<dop-1;i++)
		{
			if(elementDOF->val[element[0]][3+edgeindex*(dop-1)+i]==node)
			{
				isInedge=1;
				break;
			}
		}
	}
	
	if(isInedge==1)
	{
		getinfo4average(lambda1, lambda2, edgeNode, element[0], elements, dop, nodeindex1, average);
		if(element[1]!=-1)
		{			
			getinfo4average(lambda1, lambda2, edgeNode, element[1], elements, dop, nodeindex2, averageplus);
			average[0] += averageplus[0];
			average[1] += averageplus[1];
		}
	}
	else
	{
		if(nodeindex1<elementDOF->col)
			getinfo4average(lambda1, lambda2, edgeNode, element[0], elements, dop, nodeindex1, average);
		else{
				if(element[1]!=-1)
				{
					if(nodeindex2<elementDOF->col)
						getinfo4average(lambda1, lambda2, edgeNode, element[1], elements, dop, nodeindex2, average);
				}
			}
	}

	if(element[1]!=-1)
	{
		average[0]/=2.0;
		average[1]/=2.0;
	}
}

/**
 * \fn void getinfo4average(double lambda1, double lambda2, int *edgeNode, int element, ELEMENT *elements, int dop, int index, double *average)
 * \brief the gradient of basis function
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param *edgeNode pointer to the two indices of current edge
 * \param element the index of current element
 * \param *elements pointer to the structure of the triangulation
 * \param dop number of degrees of polynomial
 * \param index index of current node variable
 * \param *average pointer to gradient of basis function
 * \return void
 */
void getinfo4average(double lambda1, double lambda2, int *edgeNode, int element, ELEMENT *elements, int dop, int index, double *average)
{
	int i,j,l;
	double lambdas[3];
	double s, *eta, *xi;
		
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;
	
	s=elements->vol[element];
	xi=elements->xi[element];
	eta=elements->eta[element];
	
	lagrange_basis1(lambdas, s, eta, xi, index, dop, average);
}

/**
 * \fn void jumpOperatorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for tensor
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensor(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2];
	double phi;
	double elen=edges->length[edge];

	jump[0]=0;
	jump[1]=0;

	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;

	if(elementDOF->dop==0)
	{
		if(node%3<1)
		{
			jump[0]=nv[0];
			jump[1]=0;
		}
		else if(node%3<2)
		{
			jump[0]=0;
			jump[1]=nv[1];
		}
		else
		{
			jump[0]=nv[1];
			jump[1]=nv[0];
		}
		return;
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	int lnode=node%(3*elementDOF->col);
	if(lnode<elementDOF->col)
	{
		jump[0]=phi*nv[0];
		jump[1]=0;
	}
	else if(lnode<elementDOF->col*2)
	{
		jump[0]=0;
		jump[1]=phi*nv[1];
	}
	else
	{
		jump[0]=phi*nv[1];
		jump[1]=phi*nv[0];
	}
}

/**
 * \fn void jumpOperatorTensorNormal(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for M_{nn}(\tau)
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorNormal(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], nve[2];
	double phi;
	double elen=edges->length[edge];

	*jump=0;

	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	nve[0]=edges->nvector[edge][0];
	nve[1]=edges->nvector[edge][1];
	
	int lnode=node%(3*elementDOF->col);
	if(elementDOF->dop==0)
	{
		if(lnode<elementDOF->col)
			*jump=nve[0]*nv[0];
		else if(lnode<elementDOF->col*2)
			*jump=nve[1]*nv[1];
		else
			*jump=nve[0]*nv[1]+nve[1]*nv[0];
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	if(lnode<elementDOF->col)
		*jump=phi*nve[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump=phi*nve[1]*nv[1];
	else
		*jump=phi*nve[0]*nv[1]+phi*nve[1]*nv[0];
}

/**
 * \fn void jumpOperatorTensoTangent (double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for M_{nt}(\tau)
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorTangent(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], tve[2];
	double phi;
	double elen=edges->length[edge];

	*jump=0;

	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	tve[0]=edges->tvector[edge][0];
	tve[1]=edges->tvector[edge][1];

	int lnode=node%(3*elementDOF->col);
	if(elementDOF->dop==0)
	{
		if(lnode<elementDOF->col)
			*jump=tve[0]*nv[0];
		else if(lnode<elementDOF->col*2)
			*jump=tve[1]*nv[1];
		else
			*jump=tve[0]*nv[1]+tve[1]*nv[0];
		return;
	}

	nodeindex=node%elementDOF->col;

	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex==edgeindex)
			return;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) != edgeindex)
			return;
	}
	else
		return;

	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis(lambdas, nodeindex, elementDOF->dop, &phi);

	if(lnode<elementDOF->col)
		*jump=phi*tve[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump=phi*tve[1]*nv[1];
	else
		*jump=phi*tve[0]*nv[1]+phi*tve[1]*nv[0];
}

/**
 * \fn void jumpOperatorTensorDIVPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for divergence of tensor and tangential derivative of M_{nt}(\tau)
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorDIVPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double s, *eta, *xi, nv[2], tv[2];
	double phi[3];
	double elen=edges->length[edge];

	*jump=0;

	if(elementDOF->dop<1)
		return;
	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	nodeindex=node%elementDOF->col;

	s=elements->vol[element];
	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	tv[0]=-nv[1];
	tv[1]=nv[0];
	
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis1(lambdas, s, eta, xi, nodeindex, elementDOF->dop, phi);

	int lnode=node%(3*elementDOF->col);
	if(lnode<elementDOF->col)
		*jump+=((1+tv[0]*tv[0])*phi[0] + tv[0]*tv[1]*phi[1])*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump+=(tv[0]*tv[1]*phi[0] + (1+tv[1]*tv[1])*phi[1])*nv[1];
	else
		*jump+=((tv[0]*tv[1]*phi[0] + (1+tv[1]*tv[1])*phi[1])*nv[0] + ((1+tv[0]*tv[0])*phi[0] + tv[0]*tv[1]*phi[1])*nv[1]);

	/*if(lnode<elementDOF->col)
		*jump+=phi[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump+=phi[1]*nv[1];
	else
		*jump+=phi[0]*nv[1]+phi[1]*nv[0];

	int Pt=0;
	if(nodeindex<3 && nodeindex>=0)
	{
		if(nodeindex!=edgeindex)
			Pt=1;
	}
	else if(nodeindex<3*elementDOF->dop && nodeindex>=3)
	{
		if((nodeindex-3)/(elementDOF->dop-1) == edgeindex)
			Pt=1;
	}

	if(Pt==1)
	{
		phi[2]=phi[0]*tv[0]+phi[1]*tv[1];
		if(lnode<elementDOF->col)
			*jump+=phi[2]*tv[0]*nv[0];
		else if(lnode<elementDOF->col*2)
			*jump+=phi[2]*tv[1]*nv[1];
		else
			*jump+=phi[2]*(tv[0]*nv[1]+tv[1]*nv[0]);
	}*/
}

 /**
 * \fn void jumpOperatorTensorDIV(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
 * \brief the jump operator for divergence of tensor
 * \param lambda1 the first length coordiante
 * \param lambda2 the second length coordiante
 * \param edge the index of current edge
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param node index of current node variable
 * \param *jump pointer to the result of jump operator
 * \return void
 */
void jumpOperatorTensorDIV(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, ELEMENT_DOF *elementDOF, int node, double *jump)
{
	int i,j,l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double s, *eta, *xi, nv[2], tv[2];
	double phi[3];
	double elen=edges->length[edge];

	*jump=0;

	if(elementDOF->dop<1)
		return;
	if(node<0 && node>=elementDOF->dof*3)
		return;

	element=node/(3*elementDOF->col);
	
	for(edgeindex=0;edgeindex<3;edgeindex++)
	{
		if(elementEdge->val[element][edgeindex]==edge)
			break;
	}

	if(edgeindex==3)
		return;

	nodeindex=node%elementDOF->col;

	s=elements->vol[element];
	xi=elements->xi[element];
	eta=elements->eta[element];
	nv[0]=-eta[edgeindex]/elen;
	nv[1]=xi[edgeindex]/elen;
	
	edgeNode[0]=edges->val[edge][0];
	edgeNode[1]=edges->val[edge][1];
	for(i=0;i<3;i++)
	{
		if(elements->val[element][i]==edgeNode[0])
			break;
	}
	for(j=0;j<3;j++)
	{
		if(elements->val[element][j]==edgeNode[1])
			break;
	}
	l=3-i-j;
	lambdas[i]=lambda1;
	lambdas[j]=lambda2;
	lambdas[l]=0;

	lagrange_basis1(lambdas, s, eta, xi, nodeindex, elementDOF->dop, phi);

	int lnode=node%(3*elementDOF->col);
	if(lnode<elementDOF->col)
		*jump+=phi[0]*nv[0];
	else if(lnode<elementDOF->col*2)
		*jump+=phi[1]*nv[1];
	else
		*jump+=phi[0]*nv[1]+phi[1]*nv[0];
}

/**
* \fn void jumpOperatorATensorTangent2(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
* \brief the jump of M_{tt}(Asigmah) for Hu-Zhang element
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorATensorTangent2(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
{
	int i, j, k, l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], nve[2], tv[2];
	double phi[3], sgn, value[2];
	double **tensorBasis[3];
	
	*jump = 0;
	
	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
	nve[0] = edges->nvector[edge][0];
	nve[1] = edges->nvector[edge][1];

	for (k = 0; k < 2; k++)
	{
		element = edges->val[edge][2 + k];
		if (element == -1)
			continue;

		for (nodeindex = 0; nodeindex<elementDOF->col; nodeindex++)
		{
			if (elementDOF->val[element][nodeindex] == node)
				break;
		}
		if (nodeindex == elementDOF->col)
			continue;

		for (edgeindex = 0; edgeindex<3; edgeindex++)
		{
			if (elementEdge->val[element][edgeindex] == edge)
				break;
		}

		for (i = 0; i<3; i++)
		{
			if (elements->val[element][i] == edgeNode[0])
				break;
		}
		for (j = 0; j<3; j++)
		{
			if (elements->val[element][j] == edgeNode[1])
				break;
		}
		l = 3 - i - j;
		lambdas[i] = lambda1;
		lambdas[j] = lambda2;
		lambdas[l] = 0;

		xi = elements->xi[element];
		eta = elements->eta[element];
		nv[0] = -eta[edgeindex] / elen;
		nv[1] = xi[edgeindex] / elen;
		tv[0] = -nv[1];
		tv[1] = nv[0];
		sgn = nve[0] * nv[0] + nve[1] * nv[1];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[element][i]];

		huzhang_basis(lambdas, elements->nvector[element], elements->tvector[element], tensorBasis, nodeindex, elementDOF->dop, phi);

		value[0] = phi[0] * tv[0] * tv[0] + phi[1] * tv[1] * tv[1] + 2 * phi[2] * tv[0] * tv[1];
		value[1] = phi[0] + phi[1];

		*jump += (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu)*sgn;
	} // k	
}

/**
* \fn void jumpOperatorATensorTangent2Dirichlet(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
* \brief the jump of M_{tt}(Asigmah) for Hu-Zhang element
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorATensorTangent2Dirichlet(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
{
	int i, j, k, l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double *eta, *xi, nv[2], tv[2];
	double phi[3], sgn, value[2];
	double **tensorBasis[3];

	*jump = 0;

	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
//	nve[0] = edges->nvector[edge][0];
//	nve[1] = edges->nvector[edge][1];

	for (k = 0; k < 2; k++)
	{
		element = edges->val[edge][2 + k];
		if (element == -1)
			continue;

		for (nodeindex = 0; nodeindex<elementDOF->col; nodeindex++)
		{
			if (elementDOF->val[element][nodeindex] == node)
				break;
		}
		if (nodeindex == elementDOF->col)
			continue;

		for (edgeindex = 0; edgeindex<3; edgeindex++)
		{
			if (elementEdge->val[element][edgeindex] == edge)
				break;
		}

		for (i = 0; i<3; i++)
		{
			if (elements->val[element][i] == edgeNode[0])
				break;
		}
		for (j = 0; j<3; j++)
		{
			if (elements->val[element][j] == edgeNode[1])
				break;
		}
		l = 3 - i - j;
		lambdas[i] = lambda1;
		lambdas[j] = lambda2;
		lambdas[l] = 0;

		xi = elements->xi[element];
		eta = elements->eta[element];
		nv[0] = -eta[edgeindex] / elen;
		nv[1] = xi[edgeindex] / elen;
		tv[0] = -nv[1];
		tv[1] = nv[0];
//		sgn = nve[0] * nv[0] + nve[1] * nv[1];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[element][i]];

		huzhang_basis(lambdas, elements->nvector[element], elements->tvector[element], tensorBasis, nodeindex, elementDOF->dop, phi);

		value[0] = phi[0] * tv[0] * tv[0] + phi[1] * tv[1] * tv[1] + 2 * phi[2] * tv[0] * tv[1];
		value[1] = phi[0] + phi[1];

		*jump += (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
	} // k	
}

/**
* \fn void jumpOperatorRotATensorTangentPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
* \brief the jump of rot(Asigmah) t - \partial_t(M_{nt}(Asigmah)) for Hu-Zhang element
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \param *jump pointer to the result of jump operator
* \return void
*/
void jumpOperatorRotATensorTangentPt(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
{
	int i, j, k, l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double lambdas[3];
	double s, *eta, *xi, nv[2], tv[2];
	double phi1[2], phi2[2], phi3[3][2], sgn, value[3];
	double **tensorBasis[3];

	*jump = 0;

	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];
//	nve[0] = edges->nvector[edge][0];
//	nve[1] = edges->nvector[edge][1];

	for (k = 0; k < 2; k++)
	{
		element = edges->val[edge][2 + k];
		if (element == -1)
			continue;

		for (nodeindex = 0; nodeindex<elementDOF->col; nodeindex++)
		{
			if (elementDOF->val[element][nodeindex] == node)
				break;
		}
		if (nodeindex == elementDOF->col)
			continue;

		for (edgeindex = 0; edgeindex<3; edgeindex++)
		{
			if (elementEdge->val[element][edgeindex] == edge)
				break;
		}

		for (i = 0; i<3; i++)
		{
			if (elements->val[element][i] == edgeNode[0])
				break;
		}
		for (j = 0; j<3; j++)
		{
			if (elements->val[element][j] == edgeNode[1])
				break;
		}
		l = 3 - i - j;
		lambdas[i] = lambda1;
		lambdas[j] = lambda2;
		lambdas[l] = 0;

		s = elements->vol[element];
		xi = elements->xi[element];
		eta = elements->eta[element];
		nv[0] = -eta[edgeindex] / elen;
		nv[1] = xi[edgeindex] / elen;
		tv[0] = -nv[1];
		tv[1] = nv[0];
//		sgn = nve[0] * nv[0] + nve[1] * nv[1];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[element][i]];

		huzhang_basisROT(lambdas, s, eta, xi, elements->nvector[element], elements->tvector[element], tensorBasis, nodeindex, elementDOF->dop, phi1);
		huzhang_basisCurlTrace(lambdas, s, eta, xi, elements->nvector[element], elements->tvector[element], tensorBasis, nodeindex, elementDOF->dop, phi2);

		value[0] = phi1[0] * tv[0] + phi1[1] * tv[1];
		value[1] = phi2[0] * tv[0] + phi2[1] * tv[1];

		*jump += (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);

		huzhang_basis1(lambdas, s, eta, xi, elements->nvector[element], elements->tvector[element], tensorBasis, nodeindex, elementDOF->dop, phi3);
		value[0] = phi3[0][0] * tv[0] + phi3[0][1] * tv[1];
		value[1] = phi3[1][0] * tv[0] + phi3[1][1] * tv[1];
		value[2] = phi3[2][0] * tv[0] + phi3[2][1] * tv[1];

		*jump -= (value[0] * nv[0] * tv[0] + value[1] * nv[1] * tv[1] + value[2] * (nv[0] * tv[1] + nv[1] * tv[0])) / (2 * mu);
	} // k	
}

/**
* \fn void NeumannStress(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
* \brief (sigmah)n on Neumann boundary for Hu-Zhang element
* \param lambda1 the first length coordiante
* \param lambda2 the second length coordiante
* \param edge the index of current edge
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \param node index of current node variable
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \param *jump pointer to the result of jump operator
* \return void
*/
void NeumannStress(double lambda1, double lambda2, int edge, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, int node, double lambda, double mu, double *jump)
{
	int i, j, l;
	int nodeindex, edgeindex;
	int element, edgeNode[2];
	double x, y, lambdas[3], xs[3], ys[3];
	double *eta, *xi, nv[2], tv[2];
	double phi[3], value[2];
	double **tensorBasis[3];

	jump[0] = 0;
	jump[1] = 0;

	if (edges->bdFlag[edge] != 1)
		return;

	if (node<0 && node >= elementDOF->dof)
		return;

	double elen = edges->length[edge];
	edgeNode[0] = edges->val[edge][0];
	edgeNode[1] = edges->val[edge][1];

	element = edges->val[edge][2];

	for (nodeindex = 0; nodeindex < elementDOF->col; nodeindex++)
	{
		if (elementDOF->val[element][nodeindex] == node)
			break;
	}
	if (nodeindex == elementDOF->col)
		return;

	for (edgeindex = 0; edgeindex < 3; edgeindex++)
	{
		if (elementEdge->val[element][edgeindex] == edge)
			break;
	}

	for (i = 0; i < 3; i++)
	{
		if (elements->val[element][i] == edgeNode[0])
			break;
	}
	for (j = 0; j < 3; j++)
	{
		if (elements->val[element][j] == edgeNode[1])
			break;
	}
	l = 3 - i - j;
	lambdas[i] = lambda1;
	lambdas[j] = lambda2;
	lambdas[l] = 0;

	xi = elements->xi[element];
	eta = elements->eta[element];
	nv[0] = -eta[edgeindex] / elen;
	nv[1] = xi[edgeindex] / elen;
	tv[0] = -nv[1];
	tv[1] = nv[0];
	for (i = 0; i<3; i++)
	{
		j = elements->val[element][i];
		xs[i] = nodes->val[j][0];
		ys[i] = nodes->val[j][1];
	}
	x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
	y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
	for (i = 0; i < 3; i++)
		tensorBasis[i] = nodes->tensorBasis[elements->val[element][i]];


	huzhang_basis(lambdas, elements->nvector[element], elements->tvector[element], tensorBasis, nodeindex, elementDOF->dop, phi);

	jump[0] = phi[0] * nv[0] + phi[2] * nv[1];
	jump[1] = phi[2] * nv[0] + phi[1] * nv[1];
}


/**
 * \fn void getElementEdgeGeoInfo(ELEMENT *elements, EDGE *edges, dennode *nodes)
 * \brief compute the geometric information of elements and edges
 * \param *elements pointer to triangulation: the first 3 columns store the indexes of vertices
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
 * \return void
 */
void getElementEdgeGeoInfo(ELEMENT *elements, EDGE *edges, dennode *nodes)
{
	int i,j;
	int vertex1, vertex2, vertex3;
	int element1, element2;
/*	if(!elements->h)
		free(elements->h);
	if(!edges->h)
		free(edges->h);
	elements->h=(double*)calloc(elements->row, sizeof(double));
	edges->h=(double*)calloc(edges->row, sizeof(double)); */
	
	for(i=0;i<elements->row;i++)
	{
		vertex1 = elements->val[i][0];
		vertex2 = elements->val[i][1];
		vertex3 = elements->val[i][2];
		elements->vol[i] = area(nodes->val[vertex1][0],nodes->val[vertex2][0],nodes->val[vertex3][0],nodes->val[vertex1][1],nodes->val[vertex2][1],nodes->val[vertex3][1]);
		elements->xi[i][0] = nodes->val[vertex2][0]-nodes->val[vertex3][0];
		elements->xi[i][1] = nodes->val[vertex3][0]-nodes->val[vertex1][0];
		elements->xi[i][2] = nodes->val[vertex1][0]-nodes->val[vertex2][0];
		elements->eta[i][0] = nodes->val[vertex2][1]-nodes->val[vertex3][1];
		elements->eta[i][1] = nodes->val[vertex3][1]-nodes->val[vertex1][1];
		elements->eta[i][2] = nodes->val[vertex1][1]-nodes->val[vertex2][1];
		for(j=0;j<3;j++)
		{
			elements->edgeslength[i][j]=sqrt(elements->xi[i][j]*elements->xi[i][j]+elements->eta[i][j]*elements->eta[i][j]);
			elements->nvector[i][j][0]=-elements->eta[i][j]/elements->edgeslength[i][j];
		    elements->nvector[i][j][1]=elements->xi[i][j]/elements->edgeslength[i][j];
			elements->tvector[i][j][0]=-elements->nvector[i][j][1];
			elements->tvector[i][j][1]=elements->nvector[i][j][0];
		}
	}
	
	for(i=0;i<edges->row;i++)
	{
		vertex1=edges->val[i][0];
		vertex2=edges->val[i][1];
		edges->xi[i]=nodes->val[vertex1][0]-nodes->val[vertex2][0];
		edges->eta[i]=nodes->val[vertex1][1]-nodes->val[vertex2][1];
		edges->length[i]=sqrt(edges->xi[i]*edges->xi[i]+edges->eta[i]*edges->eta[i]);
		edges->nvector[i][0]=-edges->eta[i]/edges->length[i];
	    edges->nvector[i][1]=edges->xi[i]/edges->length[i];
		edges->tvector[i][0]=-edges->nvector[i][1];
	    edges->tvector[i][1]=edges->nvector[i][0];

		element1=edges->val[i][2];
		element2=edges->val[i][3];
		if(element2>-1)
			edges->h[i]=sqrt((elements->vol[element1]+elements->vol[element2])/2.0);
		else
			edges->h[i]=sqrt(elements->vol[element1]);
	}
}

/**
* \fn void getFreenodesInfoHuZhang(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Hu-Zhang element
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                  the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoHuZhang(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;
	
	int nn = nodes->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	for (i = 0; i < nodes->row; i++)
	{
		if (nodes->bdFlag[i] == 5 || nodes->bdFlag[i] == 6 || nodes->bdFlag[i] == 45 || nodes->bdFlag[i] == 61) // Neumann boundary vertex or Dirichlet-Neumann boundary vertex
		{
			nfFlag->val[i] = 1;
			nfFlag->val[i + nn] = 1;
			nnf += 2;
		}
		else if (nodes->bdFlag[i] == 56) // Neumann-Neumann boundary vertex
		{
			nfFlag->val[i] = 1;
			nfFlag->val[i + nn] = 1;
			nfFlag->val[i + nn * 2] = 1;
			nnf += 3;
		}
	}

	for (j = 0; j<edges->row; j++)
	{
		if (edges->bdFlag[j] == 5 || edges->bdFlag[j] == 6) // Neumann boundary
		{
			for (i = 0; i < dop - 1; i++)
			{
				nfFlag->val[nn * 3 + j*(dop - 1) + i] = 1;
				nfFlag->val[nn * 3 + (ne + j)*(dop - 1) + i] = 1;
			}
			nnf += 2 * (dop - 1);
		}
	}

	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoLagrange(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Lagrange element
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
*                                 the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoLagrange(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	nnf = 0; // number of non-free nodes
	for (i = 0; i < nodes->row; i++)
	{
		j = nodes->bdFlag[i] / 10;
		k = nodes->bdFlag[i] % 10;
		if (j == 1 || j == 2 || j == 3 || j == 4 || k == 1 || k == 2 || k == 3 || k == 4)
		{
			nfFlag->val[i] = 1;
			nnf += 1;
		}
	}

	for (j = 0; j<edges->row; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			for (i = 0; i < dop - 1; i++)
				nfFlag->val[nn + j*(dop - 1) + i] = 1;
			nnf += (dop - 1);
		}
	}

	create_ivector(nnf, nfreenodes);
	create_ivector(dof - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<dof; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoNcP1(EDGE *edges, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of nonconforming P1 element
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoNcP1(EDGE *edges, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int ne = edges->row;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(ne, nfFlag);
	create_ivector(ne, index);

	for (j = 0; j<ne; j++)
	{
		if (edges->val[j][3] == -1)
		{
			nfFlag->val[j] = 1;
		}
	}

	nnf = 0; // number of non-free nodes
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1)
			nnf++;
	}
	create_ivector(nnf, nfreenodes);
	create_ivector(nfFlag->row - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoMorley(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Morley element
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
*                                 the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoMorley(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	//	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof, nfFlag);
	create_ivector(dof, index);

	for (i = 0; i < nodes->row; i++)
	{
		j = nodes->bdFlag[i] / 10;
		k = nodes->bdFlag[i] % 10;
		if (j == 1 || j == 2 || j == 3 || j == 4 || k == 1 || k == 2 || k == 3 || k == 4)
		{
			nfFlag->val[i] = 1;
		}
	}

	for (j = 0; j<edges->row; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			nfFlag->val[nn + j] = 1;
		}
	}

	nnf = 0; // number of non-free nodes
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1)
			nnf++;
	}
	create_ivector(nnf, nfreenodes);
	create_ivector(nfFlag->row - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoNcP1Vector2(EDGE *edges, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of nonconforming P1 element in 2d-vector version
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoNcP1Vector2(EDGE *edges, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int ne = edges->row;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(ne * 2, nfFlag);
	create_ivector(ne * 2, index);

	for (j = 0; j<ne; j++)
	{
		if (edges->val[j][3] == -1)
		{
			nfFlag->val[j] = 1;
			nfFlag->val[j + ne] = 1;
		}
	}

	nnf = 0; // number of non-free nodes
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1)
			nnf++;
	}
	create_ivector(nnf, nfreenodes);
	create_ivector(nfFlag->row - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoMINIVector2(dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of MINI element for Stokes equation
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoMINIVector2(dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int dof = elementDOF->dof;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof * 2, nfFlag);
	create_ivector(dof * 2, index);

	for (i = 0; i < nodes->row; i++)
	{
		j = nodes->bdFlag[i] / 10;
		k = nodes->bdFlag[i] % 10;
		if (j == 1 || j == 2 || j == 3 || j == 4 || k == 1 || k == 2 || k == 3 || k == 4)
		{
			nfFlag->val[i] = 1;
			nfFlag->val[i + dof] = 1;
		}
	}

	nnf = 0; // number of non-free nodes
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1)
			nnf++;
	}
	create_ivector(nnf, nfreenodes);
	create_ivector(nfFlag->row - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
* \fn void getFreenodesInfoCRVector2(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief get freenodes information of Crouzeix每Raviart element for Stokes equation
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
*                                 the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \return void
*/
void getFreenodesInfoCRVector2(EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int i, j, k, nnf;

	int nn = nodes->row;
	int ne = edges->row;
	int dof = elementDOF->dof;
	//	int dop = elementDOF->dop;

	ivector *nfFlag = &elementDOF->nfFlag;
	ivector *freenodes = &elementDOF->freenodes;
	ivector *nfreenodes = &elementDOF->nfreenodes;
	ivector *index = &elementDOF->index;

	create_ivector(dof * 2, nfFlag);
	create_ivector(dof * 2, index);

	for (i = 0; i < nodes->row; i++)
	{
		j = nodes->bdFlag[i] / 10;
		k = nodes->bdFlag[i] % 10;
		if (j == 1 || j == 2 || j == 3 || j == 4 || k == 1 || k == 2 || k == 3 || k == 4)
		{
			nfFlag->val[i] = 1;
			nfFlag->val[i + dof] = 1;
		}
	}

	for (j = 0; j<edges->row; j++)
	{
		if (edges->bdFlag[j] == 1 || edges->bdFlag[j] == 2 || edges->bdFlag[j] == 3 || edges->bdFlag[j] == 4) // Dirichlet boundary
		{
			nfFlag->val[nn + j] = 1;
			nfFlag->val[nn + j + dof] = 1;
		}
	}

	nnf = 0; // number of non-free nodes
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1)
			nnf++;
	}
	create_ivector(nnf, nfreenodes);
	create_ivector(nfFlag->row - nnf, freenodes);

	j = 0; k = 0;
	for (i = 0; i<nfFlag->row; i++)
	{
		if (nfFlag->val[i] == 1) //  non-free node
		{
			nfreenodes->val[k] = i;
			index->val[i] = k;
			k++;
		}
		else // free variable
		{
			freenodes->val[j] = i;
			index->val[i] = j;
			j++;
		}
	}
}

/**
 * \fn void getBoundaryInfoNormalTrace(EDGE *edges, ELEMENT_DOF *elementDOF, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
 * \brief get information of triangulation on the boundary
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *dirichlet pointer to the index of dirichlet nodes
 * \param *nondirichlet pointer to the index of nondirichlet nodes
 * \param *index pointer to the transpose of dirichlet and nondirichlet nodes
 * \return void
 */
void getBoundaryInfoNormalTrace(EDGE *edges, ELEMENT_DOF *elementDOF, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
{
	int i,j,k,dr;

	int dof=elementDOF->dof;
	create_ivector(dof, isInNode);
	create_ivector(dof, index);
	
	for(j=0;j<edges->row;j++)
	{
		if(edges->val[j][3]==-1)
		{
			for(i=0;i<elementDOF->col;i++)
				isInNode->val[elementDOF->val[j][i]]=-1;
		}
	}
	
	dr=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1)
			dr++;
	}
	create_ivector(dr, dirichlet);
	create_ivector(isInNode->row-dr, nondirichlet);

	j=0;k=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1) //  Dirichlet boundary node
		{
			dirichlet->val[k]=i;
			index->val[i]=k;
			k++;
		}
		else // free variable
		{
			nondirichlet->val[j]=i;
			index->val[i]=j;
			j++;
		}
	}
}

/**
 * \fn void getBoundaryInfoEdge(EDGE *edges, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
 * \brief get information of triangulation on the boundary
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *isInNode pointer to boundary information of nodes: if the node is interior node, it will be 0; if the node is on the boundary, it will be -1
 * \param *dirichlet pointer to the index of dirichlet nodes
 * \param *nondirichlet pointer to the index of nondirichlet nodes
 * \param *index pointer to the transpose of dirichlet and nondirichlet nodes
 * \return void
 */
void getBoundaryInfoEdge(EDGE *edges, ivector *isInNode, ivector *dirichlet, ivector *nondirichlet, ivector *index)
{
	int i,j,k,dr;

	int ne=edges->row;
	create_ivector(ne, isInNode);
	create_ivector(ne, index);
	
	for(j=0;j<ne;j++)
	{
		if(edges->val[j][3]==-1)
			isInNode->val[j]=-1;
	}
	
	dr=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1)
			dr++;
	}
	create_ivector(dr, dirichlet);
	create_ivector(isInNode->row-dr, nondirichlet);

	j=0;k=0;
	for(i=0;i<isInNode->row;i++)
	{
		if(isInNode->val[i]==-1) //  Dirichlet boundary node
		{
			dirichlet->val[k]=i;
			index->val[i]=k;
			k++;
		}
		else // free variable
		{
			nondirichlet->val[j]=i;
			index->val[i]=j;
			j++;
		}
	}
}

/**
* \fn void getTensorBasis4HuZhang2d(EDGE *edges, iCSRmat *edgesTran, dennode *nodes)
* \brief get tensor basis of nodes for Hu-Zhang element in 2d
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                  the fourth column stores -1 if the edge is on boundary
* \param *edgesTran the relation between nodes and edges. JA stores edge index, A stores another vertex
* \param *nodes pointer to the nodes location of the triangulation
* \return void
*/
void getTensorBasis4HuZhang2d(EDGE *edges, iCSRmat *edgesTran, dennode *nodes)
{

	int i, j;
	int edge;
	double *nv, *tv, **tensorBasis;

	for (i = 0; i < nodes->row; i++)
	{
		tensorBasis = nodes->tensorBasis[i];
		tensorBasis[0][0] = 1; tensorBasis[0][1] = 0; tensorBasis[0][2] = 0;
		tensorBasis[1][0] = 0; tensorBasis[1][1] = 1; tensorBasis[1][2] = 0;
		tensorBasis[2][0] = 0; tensorBasis[2][1] = 0; tensorBasis[2][2] = 1.0 / sqrt(2);

		if (nodes->bdFlag[i] == 5 || nodes->bdFlag[i] == 6 || nodes->bdFlag[i] == 45 || nodes->bdFlag[i] == 61) // Neumann or Dirichlet-Neumann
		{
			for (j = edgesTran->IA[i]; j < edgesTran->IA[i + 1]; j++)
			{
				edge = edgesTran->JA[j];
				if (edges->bdFlag[edge] == 5 || edges->bdFlag[edge] == 6)
				{
					nv = edges->nvector[edge];
					tv = edges->tvector[edge];
					tensorBasis[0][0] = nv[0] * nv[0]; tensorBasis[0][1] = nv[1] * nv[1]; tensorBasis[0][2] = nv[0] * nv[1];
					tensorBasis[1][0] = nv[0] * tv[0] * sqrt(2); tensorBasis[1][1] = nv[1] * tv[1] * sqrt(2); tensorBasis[1][2] = (nv[0] * tv[1] + nv[1] * tv[0]) / sqrt(2);
					tensorBasis[2][0] = tv[0] * tv[0]; tensorBasis[2][1] = tv[1] * tv[1]; tensorBasis[2][2] = tv[0] * tv[1];
					break;
				}
			}
		}
	}
}


/**
* \fn int getlocaldofsEdge4HuZhang2d(int *patchnodes, int edgeindex, int dop)
* \brief get local index of dofs on edge for Hu-Zhang element in 2d
* \param *patchnodes pointer to local index of dofs on edge
* \param edgeindex local index of current edge
* \param dop dop of Hu-Zhang element in 2d
* \return void
*/
int getlocaldofsEdge4HuZhang2d(int *patchnodes, int edgeindex, int dop)
{
	int i, j;
	int l = edgeindex;
	int count = 0;
	for (i = 0; i < 3; i++)
	{
		j = (l + 1) % 3 + i * 3;
		patchnodes[count] = j;
		count++;

		j = (l + 2) % 3 + i * 3;
		patchnodes[count] = j;
		count++;
	}

	for (i = 0; i < dop - 1; i++)
	{
		j = 9 + l*(dop - 1) + i;
		patchnodes[count] = j;
		count++;

		j = 9 + (l + 3) * (dop - 1) + i;
		patchnodes[count] = j;
		count++;
	}

	return count;
}

/**
* \fn void NeumannbdValHuZhang(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief get Neumann boundary value of Hu-Zhang element
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to nodes: the first column stores the x coordinate of points, the second column stores the y coordinate of points
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void NeumannbdValHuZhang(dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int i, j, k, l;

	int nnd = nodes->row;
	int ne = edges->row;
	int dop = elementDOF->dop;

	double x, y, xs[2], ys[2];
	double **tensorBasis, *nv, *tv, nn[3], nt[3];
	int edgenode[2];

	for (i = 0; i < nodes->row; i++)
	{
		if (nodes->bdFlag[i] == 5 || nodes->bdFlag[i] == 6 || nodes->bdFlag[i] == 45 || nodes->bdFlag[i] == 61) // Neumann boundary vertex or Dirichlet-Neumann boundary vertex
		{
			x = nodes->val[i][0];
			y = nodes->val[i][1];
			tensorBasis = nodes->tensorBasis[i];
			uh->val[i] = sigma11(x, y, lambda, mu)*tensorBasis[0][0] + sigma22(x, y, lambda, mu)*tensorBasis[0][1] + 2 * sigma12(x, y, lambda, mu)*tensorBasis[0][2];
			uh->val[i + nnd] = sigma11(x, y, lambda, mu)*tensorBasis[1][0] + sigma22(x, y, lambda, mu)*tensorBasis[1][1] + 2 * sigma12(x, y, lambda, mu)*tensorBasis[1][2];
		}
		else if (nodes->bdFlag[i] == 56) // Neumann-Neumann boundary vertex
		{
			x = nodes->val[i][0];
			y = nodes->val[i][1];
			tensorBasis = nodes->tensorBasis[i];
/*			uh->val[i] = 0;
			uh->val[i + nnd] = 0;
			uh->val[i + nnd * 2] = 0;*/
			uh->val[i] = sigma11(x, y, lambda, mu)*tensorBasis[0][0] + sigma22(x, y, lambda, mu)*tensorBasis[0][1] + 2 * sigma12(x, y, lambda, mu)*tensorBasis[0][2];
			uh->val[i + nnd] = sigma11(x, y, lambda, mu)*tensorBasis[1][0] + sigma22(x, y, lambda, mu)*tensorBasis[1][1] + 2 * sigma12(x, y, lambda, mu)*tensorBasis[1][2];
			uh->val[i + nnd * 2] = sigma11(x, y, lambda, mu)*tensorBasis[2][0] + sigma22(x, y, lambda, mu)*tensorBasis[2][1] + 2 * sigma12(x, y, lambda, mu)*tensorBasis[2][2];
		}
	}

	for (j = 0; j<edges->row; j++)
	{
		if (edges->bdFlag[j] == 5 || edges->bdFlag[j] == 6) // Neumann boundary
		{
			k = edges->val[j][2];
			for (l = 0; l<3; l++)
			{
				if (elementEdge->val[k][l] == j)
					break;
			}

			edgenode[0] = elements->val[k][(l + 1) % 3];
			edgenode[1] = elements->val[k][(l + 2) % 3];

			xs[0] = nodes->val[edgenode[0]][0];
			ys[0] = nodes->val[edgenode[0]][1];
			xs[1] = nodes->val[edgenode[1]][0];
			ys[1] = nodes->val[edgenode[1]][1];

			nv = edges->nvector[j];
			tv = edges->tvector[j];
			nn[0] = nv[0] * nv[0]; nn[1] = nv[1] * nv[1]; nn[2] = nv[0] * nv[1];
			nt[0] = nv[0] * tv[0] * sqrt(2); nt[1] = nv[1] * tv[1] * sqrt(2); nt[2] = (nv[0] * tv[1] + nv[1] * tv[0]) / sqrt(2);

			for (i = 0; i < dop - 1; i++)
			{
				x = xs[0] * (double)(dop - i - 1) / (double)dop + xs[1] * (double)(i + 1) / (double)dop;
				y = ys[0] * (double)(dop - i - 1) / (double)dop + ys[1] * (double)(i + 1) / (double)dop;
				uh->val[nnd * 3 + j*(dop - 1) + i] = sigma11(x, y, lambda, mu)*nn[0] + sigma22(x, y, lambda, mu)*nn[1] + 2 * sigma12(x, y, lambda, mu)*nn[2];
				uh->val[nnd * 3 + (ne + j)*(dop - 1) + i] = sigma11(x, y, lambda, mu)*nt[0] + sigma22(x, y, lambda, mu)*nt[1] + 2 * sigma12(x, y, lambda, mu)*nt[2];
			}
		}
	}
}
