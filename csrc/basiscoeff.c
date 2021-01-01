/*
 *  basiscoeff.c
 *
 *  Created by Xuehai Huang on 7/9/2012.
 *  Copyright 2012 WZU. All rights reserved.
 *
 */

/*! \file basiscoeff.c
 *  \brief coefficients of Basis functions
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "matvec.h"


/**
 * \fn void generateBasisCoeffs(ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes)
 * \brief generate the coefficients of basis functions
 * \param *basisCoeffs pointer to the coefficients of basis functions
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \return void
 */
void generateBasisCoeffs(ddenmat3 *basisCoeffs, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes)
{
	int i, j, k, count1, count2;
	dCSRmat A;
	ddenmat B, C;
	int iter=0;
	
	A.row=basisCoeffs->col;
	A.col=A.row;
	A.IA=(int*)calloc(A.row+1, sizeof(int));
	A.IA[1]=5;A.IA[2]=5;A.IA[3]=5;
	A.IA[4]=5;A.IA[5]=5;A.IA[6]=5;
	A.IA[7]=5;A.IA[8]=5;A.IA[9]=5;
	A.IA[10]=11;A.IA[11]=11;A.IA[12]=9;A.IA[13]=9;
	A.IA[14]=11;A.IA[15]=11;A.IA[16]=9;A.IA[17]=9;
	A.IA[18]=11;A.IA[19]=11;A.IA[20]=9;A.IA[21]=9;
	A.IA[22]=10;A.IA[23]=10;A.IA[24]=10;
	for(i=0;i<A.row;i++)
		A.IA[i+1]+=A.IA[i];
	A.nnz=A.IA[A.row];
	A.JA=(int*)calloc(A.nnz,sizeof(int));
	A.val=(double*)calloc(A.nnz,sizeof(double));

	
	double x[3], y[3], nv[3][2], nve[3][2];
	double lambdas[3], xx, yy;

	int num_qp=9; // the number of numerical intergation points
	double gauss[num_qp][3];
	init_Gauss(num_qp, 2, gauss); // gauss intergation initial

	int num_qp1=4; // the number of numerical intergation points
	double gauss1[num_qp1][2];
	init_Gauss1D(num_qp1, 1, gauss1); // gauss intergation initial

	for(k=0;k<elements->row;k++)
	{
		for(i=0;i<3;i++)
		{
			j=elements->val[k][i];
			x[i]=nodes->val[j][0];
			y[i]=nodes->val[j][1];
		}
		for(i=0;i<3;i++)
		{
			j=elementEdge->val[k][i];
			nv[i][0]=-elements->eta[k][i]/edges->length[j];
			nv[i][1]=elements->xi[k][i]/edges->length[j];
			nve[i][0]=edges->nvector[j][0];
			nve[i][1]=edges->nvector[j][1];
		}

		/********************************generate A************************************/
		// row 1
		A.JA[0] = 0; A.val[0] = 1;
		A.JA[1] = 18; A.val[1] = y[0]*y[0]*y[0];
		A.JA[2] = 21; A.val[2] = -3*x[0]*y[0]*y[0];
		A.JA[3] = 22; A.val[3] = -x[0]*x[0]*x[0]/3.0;
		A.JA[4] = 23; A.val[4] = -x[0]*x[0]*y[0];
		// row 2
		A.JA[5] = 1; A.val[5] = 1;
		A.JA[6] = 19; A.val[6] = x[0]*x[0]*x[0];
		A.JA[7] = 20; A.val[7] = -3*x[0]*x[0]*y[0];
		A.JA[8] = 22; A.val[8] = -x[0]*y[0]*y[0];
		A.JA[9] = 23; A.val[9] = -y[0]*y[0]*y[0]/3.0;
		// row 3
		A.JA[10] = 2; A.val[10] = 1;
		A.JA[11] = 20; A.val[11] = x[0]*x[0]*x[0];
		A.JA[12] = 21; A.val[12] = y[0]*y[0]*y[0];
		A.JA[13] = 22; A.val[13] = x[0]*x[0]*y[0];
		A.JA[14] = 23; A.val[14] = x[0]*y[0]*y[0];
		// row 4
		A.JA[15] = 3; A.val[15] = 1;
		A.JA[16] = 18; A.val[16] = y[1]*y[1]*y[1];
		A.JA[17] = 21; A.val[17] = -3*x[1]*y[1]*y[1];
		A.JA[18] = 22; A.val[18] = -x[1]*x[1]*x[1]/3.0;
		A.JA[19] = 23; A.val[19] = -x[1]*x[1]*y[1];
		// row 5
		A.JA[20] = 4; A.val[20] = 1;
		A.JA[21] = 19; A.val[21] = x[1]*x[1]*x[1];
		A.JA[22] = 20; A.val[22] = -3*x[1]*x[1]*y[1];
		A.JA[23] = 22; A.val[23] = -x[1]*y[1]*y[1];
		A.JA[24] = 23; A.val[24] = -y[1]*y[1]*y[1]/3.0;
		// row 6
		A.JA[25] = 5; A.val[25] = 1;
		A.JA[26] = 20; A.val[26] = x[1]*x[1]*x[1];
		A.JA[27] = 21; A.val[27] = y[1]*y[1]*y[1];
		A.JA[28] = 22; A.val[28] = x[1]*x[1]*y[1];
		A.JA[29] = 23; A.val[29] = x[1]*y[1]*y[1];
		// row 7
		A.JA[30] = 6; A.val[30] = 1;
		A.JA[31] = 18; A.val[31] = y[2]*y[2]*y[2];
		A.JA[32] = 21; A.val[32] = -3*x[2]*y[2]*y[2];
		A.JA[33] = 22; A.val[33] = -x[2]*x[2]*x[2]/3.0;
		A.JA[34] = 23; A.val[34] = -x[2]*x[2]*y[2];
		// row 8
		A.JA[35] = 7; A.val[35] = 1;
		A.JA[36] = 19; A.val[36] = x[2]*x[2]*x[2];
		A.JA[37] = 20; A.val[37] = -3*x[2]*x[2]*y[2];
		A.JA[38] = 22; A.val[38] = -x[2]*y[2]*y[2];
		A.JA[39] = 23; A.val[39] = -y[2]*y[2]*y[2]/3.0;
		// row 9
		A.JA[40] = 8; A.val[40] = 1;
		A.JA[41] = 20; A.val[41] = x[2]*x[2]*x[2];
		A.JA[42] = 21; A.val[42] = y[2]*y[2]*y[2];
		A.JA[43] = 22; A.val[43] = x[2]*x[2]*y[2];
		A.JA[44] = 23; A.val[44] = x[2]*y[2]*y[2];
		// row 10
		A.JA[45] = 3; A.val[45] = nve[0][0]/2.0;
		A.JA[46] = 6; A.val[46] = nve[0][0]/2.0;
		A.JA[47] = 9; A.val[47] = nve[0][0]/6.0;
		A.JA[48] = 5; A.val[48] = nve[0][1]/2.0;
		A.JA[49] = 8; A.val[49] = nve[0][1]/2.0;
		A.JA[50] = 11; A.val[50] = nve[0][1]/6.0;
		A.JA[51] = 18; A.val[51] = 0;
		A.JA[52] = 20; A.val[52] = 0;
		A.JA[53] = 21; A.val[53] = 0;
		A.JA[54] = 22; A.val[54] = 0;
		A.JA[55] = 23; A.val[55] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[2]=gauss1[i][0];
			lambdas[1]=1-lambdas[2];
			xx=x[1]*lambdas[1]+x[2]*lambdas[2];
			yy=y[1]*lambdas[1]+y[2]*lambdas[2];
			A.val[51] += gauss1[i][1]*pow(yy,3)*nve[0][0];
			A.val[52] += gauss1[i][1]*pow(xx,3)*nve[0][1];
			A.val[53] += gauss1[i][1]*(-3*xx*yy*yy*nve[0][0] + pow(yy,3)*nve[0][1]);
			A.val[54] += gauss1[i][1]*(-pow(xx,3)*nve[0][0]/3.0 + xx*xx*yy*nve[0][1]);
			A.val[55] += gauss1[i][1]*(-xx*xx*yy*nve[0][0] + xx*yy*yy*nve[0][1]);
		}
		// row 11
		A.JA[56] = 4; A.val[56] = nve[0][1]/2.0;
		A.JA[57] = 7; A.val[57] = nve[0][1]/2.0;
		A.JA[58] = 10; A.val[58] = nve[0][1]/6.0;
		A.JA[59] = 5; A.val[59] = nve[0][0]/2.0;
		A.JA[60] = 8; A.val[60] = nve[0][0]/2.0;
		A.JA[61] = 11; A.val[61] = nve[0][0]/6.0;
		A.JA[62] = 19; A.val[62] = 0;
		A.JA[63] = 20; A.val[63] = 0;
		A.JA[64] = 21; A.val[64] = 0;
		A.JA[65] = 22; A.val[65] = 0;
		A.JA[66] = 23; A.val[66] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[2]=gauss1[i][0];
			lambdas[1]=1-lambdas[2];
			xx=x[1]*lambdas[1]+x[2]*lambdas[2];
			yy=y[1]*lambdas[1]+y[2]*lambdas[2];
			A.val[62] += gauss1[i][1]*pow(xx,3)*nve[0][1];
			A.val[63] += gauss1[i][1]*(pow(xx,3)*nve[0][0] - 3*xx*xx*yy*nve[0][1]);
			A.val[64] += gauss1[i][1]*pow(yy,3)*nve[0][0];
			A.val[65] += gauss1[i][1]*(xx*xx*yy*nve[0][0] - xx*yy*yy*nve[0][1]);
			A.val[66] += gauss1[i][1]*(xx*yy*yy*nve[0][0] - pow(yy,3)*nve[0][1]/3.0);
		}
		// row 12
		A.JA[67] = 3; A.val[67] = -nv[0][0]/12.0;
		A.JA[68] = 6; A.val[68] = nv[0][0]/12.0;
		A.JA[69] = 5; A.val[69] = -nv[0][1]/12.0;
		A.JA[70] = 8; A.val[70] = nv[0][1]/12.0;
		A.JA[71] = 18; A.val[71] = 0;
		A.JA[72] = 20; A.val[72] = 0;
		A.JA[73] = 21; A.val[73] = 0;
		A.JA[74] = 22; A.val[74] = 0;
		A.JA[75] = 23; A.val[75] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[2]=gauss1[i][0];
			lambdas[1]=1-lambdas[2];
			xx=x[1]*lambdas[1]+x[2]*lambdas[2];
			yy=y[1]*lambdas[1]+y[2]*lambdas[2];
			A.val[71] += gauss1[i][1]*(lambdas[2]-0.5)*pow(yy,3)*nv[0][0];
			A.val[72] += gauss1[i][1]*(lambdas[2]-0.5)*pow(xx,3)*nv[0][1];
			A.val[73] += gauss1[i][1]*(lambdas[2]-0.5)*(-3*xx*yy*yy*nv[0][0] + pow(yy,3)*nv[0][1]);
			A.val[74] += gauss1[i][1]*(lambdas[2]-0.5)*(-pow(xx,3)*nv[0][0]/3.0 + xx*xx*yy*nv[0][1]);
			A.val[75] += gauss1[i][1]*(lambdas[2]-0.5)*(-xx*xx*yy*nv[0][0] + xx*yy*yy*nv[0][1]);
		}
		// row 13
		A.JA[76] = 4; A.val[76] = -nv[0][1]/12.0;
		A.JA[77] = 7; A.val[77] = nv[0][1]/12.0;
		A.JA[78] = 5; A.val[78] = -nv[0][0]/12.0;
		A.JA[79] = 8; A.val[79] = nv[0][0]/12.0;
		A.JA[80] = 19; A.val[80] = 0;
		A.JA[81] = 20; A.val[81] = 0;
		A.JA[82] = 21; A.val[82] = 0;
		A.JA[83] = 22; A.val[83] = 0;
		A.JA[84] = 23; A.val[84] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[2]=gauss1[i][0];
			lambdas[1]=1-lambdas[2];
			xx=x[1]*lambdas[1]+x[2]*lambdas[2];
			yy=y[1]*lambdas[1]+y[2]*lambdas[2];
			A.val[80] += gauss1[i][1]*(lambdas[2]-0.5)*pow(xx,3)*nv[0][1];
			A.val[81] += gauss1[i][1]*(lambdas[2]-0.5)*(pow(xx,3)*nv[0][0] - 3*xx*xx*yy*nv[0][1]);
			A.val[82] += gauss1[i][1]*(lambdas[2]-0.5)*pow(yy,3)*nv[0][0];
			A.val[83] += gauss1[i][1]*(lambdas[2]-0.5)*(xx*xx*yy*nv[0][0] - xx*yy*yy*nv[0][1]);
			A.val[84] += gauss1[i][1]*(lambdas[2]-0.5)*(xx*yy*yy*nv[0][0] - pow(yy,3)*nv[0][1]/3.0);
		}
		// row 14
		A.JA[85] = 0; A.val[85] = nve[1][0]/2.0;
		A.JA[86] = 6; A.val[86] = nve[1][0]/2.0;
		A.JA[87] = 12; A.val[87] = nve[1][0]/6.0;
		A.JA[88] = 2; A.val[88] = nve[1][1]/2.0;
		A.JA[89] = 8; A.val[89] = nve[1][1]/2.0;
		A.JA[90] = 14; A.val[90] = nve[1][1]/6.0;
		A.JA[91] = 18; A.val[91] = 0;
		A.JA[92] = 20; A.val[92] = 0;
		A.JA[93] = 21; A.val[93] = 0;
		A.JA[94] = 22; A.val[94] = 0;
		A.JA[95] = 23; A.val[95] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[0]=gauss1[i][0];
			lambdas[2]=1-lambdas[0];
			xx=x[0]*lambdas[0]+x[2]*lambdas[2];
			yy=y[0]*lambdas[0]+y[2]*lambdas[2];
			A.val[91] += gauss1[i][1]*pow(yy,3)*nve[1][0];
			A.val[92] += gauss1[i][1]*pow(xx,3)*nve[1][1];
			A.val[93] += gauss1[i][1]*(-3*xx*yy*yy*nve[1][0] + pow(yy,3)*nve[1][1]);
			A.val[94] += gauss1[i][1]*(-pow(xx,3)*nve[1][0]/3.0 + xx*xx*yy*nve[1][1]);
			A.val[95] += gauss1[i][1]*(-xx*xx*yy*nve[1][0] + xx*yy*yy*nve[1][1]);
		}
		// row 15
		A.JA[96] = 1; A.val[96] = nve[1][1]/2.0;
		A.JA[97] = 7; A.val[97] = nve[1][1]/2.0;
		A.JA[98] = 13; A.val[98] = nve[1][1]/6.0;
		A.JA[99] = 2; A.val[99] = nve[1][0]/2.0;
		A.JA[100] = 8; A.val[100] = nve[1][0]/2.0;
		A.JA[101] = 14; A.val[101] = nve[1][0]/6.0;
		A.JA[102] = 19; A.val[102] = 0;
		A.JA[103] = 20; A.val[103] = 0;
		A.JA[104] = 21; A.val[104] = 0;
		A.JA[105] = 22; A.val[105] = 0;
		A.JA[106] = 23; A.val[106] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[0]=gauss1[i][0];
			lambdas[2]=1-lambdas[0];
			xx=x[0]*lambdas[0]+x[2]*lambdas[2];
			yy=y[0]*lambdas[0]+y[2]*lambdas[2];
			A.val[102] += gauss1[i][1]*pow(xx,3)*nve[1][1];
			A.val[103] += gauss1[i][1]*(pow(xx,3)*nve[1][0] - 3*xx*xx*yy*nve[1][1]);
			A.val[104] += gauss1[i][1]*pow(yy,3)*nve[1][0];
			A.val[105] += gauss1[i][1]*(xx*xx*yy*nve[1][0] - xx*yy*yy*nve[1][1]);
			A.val[106] += gauss1[i][1]*(xx*yy*yy*nve[1][0] - pow(yy,3)*nve[1][1]/3.0);
		}
		// row 16
		A.JA[107] = 0; A.val[107] = nv[1][0]/12.0;
		A.JA[108] = 6; A.val[108] = -nv[1][0]/12.0;
		A.JA[109] = 2; A.val[109] = nv[1][1]/12.0;
		A.JA[110] = 8; A.val[110] = -nv[1][1]/12.0;
		A.JA[111] = 18; A.val[111] = 0;
		A.JA[112] = 20; A.val[112] = 0;
		A.JA[113] = 21; A.val[113] = 0;
		A.JA[114] = 22; A.val[114] = 0;
		A.JA[115] = 23; A.val[115] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[0]=gauss1[i][0];
			lambdas[2]=1-lambdas[0];
			xx=x[0]*lambdas[0]+x[2]*lambdas[2];
			yy=y[0]*lambdas[0]+y[2]*lambdas[2];
			A.val[111] += gauss1[i][1]*(lambdas[0]-0.5)*pow(yy,3)*nv[1][0];
			A.val[112] += gauss1[i][1]*(lambdas[0]-0.5)*pow(xx,3)*nv[1][1];
			A.val[113] += gauss1[i][1]*(lambdas[0]-0.5)*(-3*xx*yy*yy*nv[1][0] + pow(yy,3)*nv[1][1]);
			A.val[114] += gauss1[i][1]*(lambdas[0]-0.5)*(-pow(xx,3)*nv[1][0]/3.0 + xx*xx*yy*nv[1][1]);
			A.val[115] += gauss1[i][1]*(lambdas[0]-0.5)*(-xx*xx*yy*nv[1][0] + xx*yy*yy*nv[1][1]);
		}
		// row 17
		A.JA[116] = 1; A.val[116] = nv[1][1]/12.0;
		A.JA[117] = 7; A.val[117] = -nv[1][1]/12.0;
		A.JA[118] = 2; A.val[118] = nv[1][0]/12.0;
		A.JA[119] = 8; A.val[119] = -nv[1][0]/12.0;
		A.JA[120] = 19; A.val[120] = 0;
		A.JA[121] = 20; A.val[121] = 0;
		A.JA[122] = 21; A.val[122] = 0;
		A.JA[123] = 22; A.val[123] = 0;
		A.JA[124] = 23; A.val[124] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[0]=gauss1[i][0];
			lambdas[2]=1-lambdas[0];
			xx=x[0]*lambdas[0]+x[2]*lambdas[2];
			yy=y[0]*lambdas[0]+y[2]*lambdas[2];
			A.val[120] += gauss1[i][1]*(lambdas[0]-0.5)*pow(xx,3)*nv[1][1];
			A.val[121] += gauss1[i][1]*(lambdas[0]-0.5)*(pow(xx,3)*nv[1][0] - 3*xx*xx*yy*nv[1][1]);
			A.val[122] += gauss1[i][1]*(lambdas[0]-0.5)*pow(yy,3)*nv[1][0];
			A.val[123] += gauss1[i][1]*(lambdas[0]-0.5)*(xx*xx*yy*nv[1][0] - xx*yy*yy*nv[1][1]);
			A.val[124] += gauss1[i][1]*(lambdas[0]-0.5)*(xx*yy*yy*nv[1][0] - pow(yy,3)*nv[1][1]/3.0);
		}
		// row 18
		A.JA[125] = 0; A.val[125] = nve[2][0]/2.0;
		A.JA[126] = 3; A.val[126] = nve[2][0]/2.0;
		A.JA[127] = 15; A.val[127] = nve[2][0]/6.0;
		A.JA[128] = 2; A.val[128] = nve[2][1]/2.0;
		A.JA[129] = 5; A.val[129] = nve[2][1]/2.0;
		A.JA[130] = 17; A.val[130] = nve[2][1]/6.0;
		A.JA[131] = 18; A.val[131] = 0;
		A.JA[132] = 20; A.val[132] = 0;
		A.JA[133] = 21; A.val[133] = 0;
		A.JA[134] = 22; A.val[134] = 0;
		A.JA[135] = 23; A.val[135] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[1]=gauss1[i][0];
			lambdas[0]=1-lambdas[1];
			xx=x[0]*lambdas[0]+x[1]*lambdas[1];
			yy=y[0]*lambdas[0]+y[1]*lambdas[1];
			A.val[131] += gauss1[i][1]*pow(yy,3)*nve[2][0];
			A.val[132] += gauss1[i][1]*pow(xx,3)*nve[2][1];
			A.val[133] += gauss1[i][1]*(-3*xx*yy*yy*nve[2][0] + pow(yy,3)*nve[2][1]);
			A.val[134] += gauss1[i][1]*(-pow(xx,3)*nve[2][0]/3.0 + xx*xx*yy*nve[2][1]);
			A.val[135] += gauss1[i][1]*(-xx*xx*yy*nve[2][0] + xx*yy*yy*nve[2][1]);
		}
		// row 19
		A.JA[136] = 1; A.val[136] = nve[2][1]/2.0;
		A.JA[137] = 4; A.val[137] = nve[2][1]/2.0;
		A.JA[138] = 16; A.val[138] = nve[2][1]/6.0;
		A.JA[139] = 2; A.val[139] = nve[2][0]/2.0;
		A.JA[140] = 5; A.val[140] = nve[2][0]/2.0;
		A.JA[141] = 17; A.val[141] = nve[2][0]/6.0;
		A.JA[142] = 19; A.val[142] = 0;
		A.JA[143] = 20; A.val[143] = 0;
		A.JA[144] = 21; A.val[144] = 0;
		A.JA[145] = 22; A.val[145] = 0;
		A.JA[146] = 23; A.val[146] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[1]=gauss1[i][0];
			lambdas[0]=1-lambdas[1];
			xx=x[0]*lambdas[0]+x[1]*lambdas[1];
			yy=y[0]*lambdas[0]+y[1]*lambdas[1];
			A.val[142] += gauss1[i][1]*pow(xx,3)*nve[2][1];
			A.val[143] += gauss1[i][1]*(pow(xx,3)*nve[2][0] - 3*xx*xx*yy*nve[2][1]);
			A.val[144] += gauss1[i][1]*pow(yy,3)*nve[2][0];
			A.val[145] += gauss1[i][1]*(xx*xx*yy*nve[2][0] - xx*yy*yy*nve[2][1]);
			A.val[146] += gauss1[i][1]*(xx*yy*yy*nve[2][0] - pow(yy,3)*nve[2][1]/3.0);
		}
		// row 20
		A.JA[147] = 0; A.val[147] = -nv[2][0]/12.0;
		A.JA[148] = 3; A.val[148] = nv[2][0]/12.0;
		A.JA[149] = 2; A.val[149] = -nv[2][1]/12.0;
		A.JA[150] = 5; A.val[150] = nv[2][1]/12.0;
		A.JA[151] = 18; A.val[151] = 0;
		A.JA[152] = 20; A.val[152] = 0;
		A.JA[153] = 21; A.val[153] = 0;
		A.JA[154] = 22; A.val[154] = 0;
		A.JA[155] = 23; A.val[155] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[1]=gauss1[i][0];
			lambdas[0]=1-lambdas[1];
			xx=x[0]*lambdas[0]+x[1]*lambdas[1];
			yy=y[0]*lambdas[0]+y[1]*lambdas[1];
			A.val[151] += gauss1[i][1]*(lambdas[1]-0.5)*pow(yy,3)*nv[2][0];
			A.val[152] += gauss1[i][1]*(lambdas[1]-0.5)*pow(xx,3)*nv[2][1];
			A.val[153] += gauss1[i][1]*(lambdas[1]-0.5)*(-3*xx*yy*yy*nv[2][0] + pow(yy,3)*nv[2][1]);
			A.val[154] += gauss1[i][1]*(lambdas[1]-0.5)*(-pow(xx,3)*nv[2][0]/3.0 + xx*xx*yy*nv[2][1]);
			A.val[155] += gauss1[i][1]*(lambdas[1]-0.5)*(-xx*xx*yy*nv[2][0] + xx*yy*yy*nv[2][1]);
		}
		// row 21
		A.JA[156] = 1; A.val[156] = -nv[2][1]/12.0;
		A.JA[157] = 4; A.val[157] = nv[2][1]/12.0;
		A.JA[158] = 2; A.val[158] = -nv[2][0]/12.0;
		A.JA[159] = 5; A.val[159] = nv[2][0]/12.0;
		A.JA[160] = 19; A.val[160] = 0;
		A.JA[161] = 20; A.val[161] = 0;
		A.JA[162] = 21; A.val[162] = 0;
		A.JA[163] = 22; A.val[163] = 0;
		A.JA[164] = 23; A.val[164] = 0;
		for(i=0;i<num_qp1;i++)
		{
			lambdas[1]=gauss1[i][0];
			lambdas[0]=1-lambdas[1];
			xx=x[0]*lambdas[0]+x[1]*lambdas[1];
			yy=y[0]*lambdas[0]+y[1]*lambdas[1];
			A.val[160] += gauss1[i][1]*(lambdas[1]-0.5)*pow(xx,3)*nv[2][1];
			A.val[161] += gauss1[i][1]*(lambdas[1]-0.5)*(pow(xx,3)*nv[2][0] - 3*xx*xx*yy*nv[2][1]);
			A.val[162] += gauss1[i][1]*(lambdas[1]-0.5)*pow(yy,3)*nv[2][0];
			A.val[163] += gauss1[i][1]*(lambdas[1]-0.5)*(xx*xx*yy*nv[2][0] - xx*yy*yy*nv[2][1]);
			A.val[164] += gauss1[i][1]*(lambdas[1]-0.5)*(xx*yy*yy*nv[2][0] - pow(yy,3)*nv[2][1]/3.0);
		}
		// row 22
		A.JA[165] = 0; A.val[165] = 1.0/3.0;
		A.JA[166] = 3; A.val[166] = 1.0/3.0;
		A.JA[167] = 6; A.val[167] = 1.0/3.0;
		A.JA[168] = 9; A.val[168] = 1.0/12.0;
		A.JA[169] = 12; A.val[169] = 1.0/12.0;
		A.JA[170] = 15; A.val[170] = 1.0/12.0;
		A.JA[171] = 18; A.val[171] = 0;
		A.JA[172] = 21; A.val[172] = 0;
		A.JA[173] = 22; A.val[173] = 0;
		A.JA[174] = 23; A.val[174] = 0;
		for(i=0;i<num_qp;i++)
		{
			lambdas[0]=gauss[i][0];
			lambdas[1]=gauss[i][1];
			lambdas[2]=1-lambdas[0]-lambdas[1];
			xx=x[0]*lambdas[0]+x[1]*lambdas[1]+x[2]*lambdas[2];
			yy=y[0]*lambdas[0]+y[1]*lambdas[1]+y[2]*lambdas[2];
			A.val[171] += 2*gauss[i][2]*pow(yy,3);
			A.val[172] += 2*gauss[i][2]*(-3*xx*yy*yy);
			A.val[173] += 2*gauss[i][2]*(-pow(xx,3)/3.0);
			A.val[174] += 2*gauss[i][2]*(-xx*xx*yy);
		}
		// row 23
		A.JA[175] = 1; A.val[175] = 1.0/3.0;
		A.JA[176] = 4; A.val[176] = 1.0/3.0;
		A.JA[177] = 7; A.val[177] = 1.0/3.0;
		A.JA[178] = 10; A.val[178] = 1.0/12.0;
		A.JA[179] = 13; A.val[179] = 1.0/12.0;
		A.JA[180] = 16; A.val[180] = 1.0/12.0;
		A.JA[181] = 19; A.val[181] = 0;
		A.JA[182] = 20; A.val[182] = 0;
		A.JA[183] = 22; A.val[183] = 0;
		A.JA[184] = 23; A.val[184] = 0;
		for(i=0;i<num_qp;i++)
		{
			lambdas[0]=gauss[i][0];
			lambdas[1]=gauss[i][1];
			lambdas[2]=1-lambdas[0]-lambdas[1];
			xx=x[0]*lambdas[0]+x[1]*lambdas[1]+x[2]*lambdas[2];
			yy=y[0]*lambdas[0]+y[1]*lambdas[1]+y[2]*lambdas[2];
			A.val[181] += 2*gauss[i][2]*pow(xx,3);
			A.val[182] += 2*gauss[i][2]*(-3*xx*xx*yy);
			A.val[183] += 2*gauss[i][2]*(-xx*yy*yy);
			A.val[184] += 2*gauss[i][2]*(-pow(yy,3)/3.0);
		}
		// row 24
		A.JA[185] = 2; A.val[185] = 1.0/3.0;
		A.JA[186] = 5; A.val[186] = 1.0/3.0;
		A.JA[187] = 8; A.val[187] = 1.0/3.0;
		A.JA[188] = 11; A.val[188] = 1.0/12.0;
		A.JA[189] = 14; A.val[189] = 1.0/12.0;
		A.JA[190] = 17; A.val[190] = 1.0/12.0;
		A.JA[191] = 20; A.val[191] = 0;
		A.JA[192] = 21; A.val[192] = 0;
		A.JA[193] = 22; A.val[193] = 0;
		A.JA[194] = 23; A.val[194] = 0;
		for(i=0;i<num_qp;i++)
		{
			lambdas[0]=gauss[i][0];
			lambdas[1]=gauss[i][1];
			lambdas[2]=1-lambdas[0]-lambdas[1];
			xx=x[0]*lambdas[0]+x[1]*lambdas[1]+x[2]*lambdas[2];
			yy=y[0]*lambdas[0]+y[1]*lambdas[1]+y[2]*lambdas[2];
			A.val[191] += 2*gauss[i][2]*pow(xx,3);
			A.val[192] += 2*gauss[i][2]*pow(yy,3);
			A.val[193] += 2*gauss[i][2]*(xx*xx*yy);
			A.val[194] += 2*gauss[i][2]*(xx*yy*yy);
		}
		/********************************form A End************************************/

		/********************************  full A  ************************************/
		//   B=full(A)
		create_dden_matrix(A.row, A.col, &B);
		for(i=0;i<A.row;i++)
		{
			for(j=A.IA[i];j<A.IA[i+1];j++)
				B.val[i][A.JA[j]]=A.val[j];
		}
		/********************************full A End************************************/

		create_dden_matrix(A.row, basisCoeffs->row, &C);
		for(i=0;i<C.col && i< C.row;i++)
			C.val[i][i]=1;
	
		AxBrref(&B, &C);
		for(i=0;i<basisCoeffs->row;i++)
		{
			for(j=0;j<basisCoeffs->col;j++)
				basisCoeffs->val[k][i][j]=C.val[j][i];
		}		

		/*********************************************************************************************
		FILE *outputFile;
		
		if(k==1)
		{
			outputFile=fopen("output/basisCoeffA.dat", "w");
			for(i=0;i<A.row;i++)
			{
				for(j=A.IA[i];j<A.IA[i+1];j++)
					fprintf(outputFile,"%d\t%d\t%e\n",i+1,A.JA[j]+1,A.val[j]);
			}
			fclose(outputFile);
			
			outputFile=fopen("output/basisCoeffs.dat", "w");
			for(j=0;j<basisCoeffs->col;j++)
			{
				for(i=0;i<basisCoeffs->row-1;i++)
					fprintf(outputFile,"%e\t",basisCoeffs->val[k][i][j]);
				fprintf(outputFile,"%e\n",basisCoeffs->val[k][i][j]);
			}
			fclose(outputFile);
		}
		********************************************************************************************/
		
		free_dden_matrix(&B);
		free_dden_matrix(&C);
	} // k
	free_csr_matrix(&A);
}

/**
 * \fn void generateBasisCoeffsEdgeProj(dBDmat *B, dBDmat *bC, EDGE *edges, int dopl, int dopk)
 * \brief generate the coefficients of basis functions caused by projection from space with dop of k onto space with dop of l on edges
 * \param *B pointer to the right hand side matrix
 * \param *bC pointer to the coefficients of basis functions
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
                                   the fourth column stores -1 if the edge is on boundary
 * \param dopl degree of polynomial of  "noto" space 
 * \param dopk degree of polynomial of  "from" space
 * \return void
 */
void generateBasisCoeffsEdgeProj(dBDmat *B, dBDmat *bC, EDGE *edges, int dopl, int dopk)
{
	dBDmat M;
	int i,k1,k2,i1;
	int ne=edges->row;
	int edge;
	
	double phi1, phi2, lambda, h;	

	create_dbd_matrix((dopl+1)*ne, (dopl+1)*ne, ne, &M);
	for(i=0;i<M.nb;i++)
		create_dden_matrix(dopl+1, dopl+1, M.blk+i);

	create_dbd_matrix((dopl+1)*ne, (dopk+1)*ne, ne, B);
	for(i=0;i<B->nb;i++)
		create_dden_matrix(dopl+1, dopk+1, B->blk+i);

	create_dbd_matrix((dopl+1)*ne, (dopk+1)*ne, ne, bC);
	for(i=0;i<bC->nb;i++)
		create_dden_matrix(dopl+1, dopk+1, bC->blk+i);

	int num_qp11=getNumQuadPoints(dopl*2, 1); // the number of numerical intergation points
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial

	int num_qp12=getNumQuadPoints(dopl+dopk, 1); // the number of numerical intergation points
	double gauss12[num_qp12][2];
	init_Gauss1D(num_qp12, 1, gauss12); // gauss intergation initial

	for(edge=0;edge<ne;edge++)
	{
		h=edges->length[edge];
	    // form M
		for(k1=0;k1<dopl+1;k1++)
		{
			for(k2=0;k2<dopl+1;k2++)
			{
				for(i1=0;i1<num_qp11;i1++)
				{
					lambda=gauss11[i1][0];
					lagrange1D_basis(lambda, k1, dopl, &phi1);
					lagrange1D_basis(lambda, k2, dopl, &phi2);
					M.blk[edge].val[k1][k2]+=h*gauss11[i1][1]*phi1*phi2; 
				}
			} // k2
		} // k1

		// form B
		for(k1=0;k1<dopl+1;k1++)
		{
			for(k2=0;k2<dopk+1;k2++)
			{
				for(i1=0;i1<num_qp12;i1++)
				{
					lambda=gauss12[i1][0];
					lagrange1D_basis(lambda, k1, dopl, &phi1);
					lagrange1D_basis(lambda, k2, dopk, &phi2);
					B->blk[edge].val[k1][k2]+=h*gauss12[i1][1]*phi1*phi2; 
				}
				bC->blk[edge].val[k1][k2]=B->blk[edge].val[k1][k2];
			} // k2
		} // k1

		AxBrref(M.blk+edge, bC->blk+edge);
	} // k

	free_dbd_matrix(&M);
}