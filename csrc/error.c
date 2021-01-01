/*
 *  error.c
 *  FEM
 *
 *  Created by Xuehai Huang on 06/10/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file error.c
 *  \brief compute error between numerical solution and exact solution
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "header.h"
#include "matvec.h"

 /**
 * \fn void geterrors(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
 * \brief compute error between numerical solution and exact solution: L2 norm, Energy norm
 * \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
 * \param *sigmah pointer to numerical solution
 * \param *uh pointer to numerical solution
 * \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
 * \param *uhstar pointer to postprocessed displacement
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
 the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 */
void geterrors(double *errors, dvector *sigmah, dvector *uh, dvector *Qhu, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int Nt = elements->row;

	int i, j, k, l;

	for (i = 0; i < 10; i++)
		errors[i] = 0;

	double phi0, phi1[3], phi2[3];
	int k1, i1, j1, l1, l2;
	double x, y, xs[3], ys[3], lambdas[3], s, *eta, *xi;
	double value[4];
	double **tensorBasis[3];

	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	int num_qp1 = 49; // the number of numerical intergation points
	double gauss1[num_qp1][3];
	init_Gauss(num_qp1, 2, gauss1); // gauss intergation initial

	int num_qp2 = 49; // the number of numerical intergation points
	double gauss2[num_qp2][3];
	init_Gauss(num_qp2, 2, gauss2); // gauss intergation initial

	int num_qp3 = getNumQuadPoints(elementDOF[1].dop * 2 - 2, 2); // the number of numerical intergation points
	double gauss3[num_qp3][3];
	init_Gauss(num_qp3, 2, gauss3); // gauss intergation initial

	int num_qp4 = getNumQuadPoints(elementDOF[1].dop * 2, 2); // the number of numerical intergation points
	double gauss4[num_qp4][3];
	init_Gauss(num_qp4, 2, gauss4); // gauss intergation initial

	int num_qp11 = getNumQuadPoints(elementDOF[1].dop * 2, 1); // the number of numerical intergation points
	if (num_qp11>5)
		num_qp11 = 5;
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial

	for (k = 0; k<Nt; k++)
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
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		// Energy norm and L2 norm of sigma-sigma_h
		for (i1 = 0; i1<num_qp0; i1++)
		{
			lambdas[0] = gauss0[i1][0];
			lambdas[1] = gauss0[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
				value[2] += phi1[2] * sigmah->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= sigma11(x, y, lambda, mu);
			value[1] -= sigma22(x, y, lambda, mu);
			value[2] -= sigma12(x, y, lambda, mu);
			value[3] = value[0] + value[1];
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			value[3] = value[3] * value[3];
			if (lambda>-0.5)
				errors[0] += 2 * s*gauss0[i1][2] * (value[0] + value[1] + 2 * value[2] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			else
				errors[0] += 2 * s*gauss0[i1][2] * (value[0] + value[1] + 2 * value[2] - value[3] * 1.0 / 2.0) / (2 * mu);

			errors[9] += 2 * s*gauss0[i1][2] * (value[0] + value[1] + 2 * value[2]);
		}

		// L2 norm of divergence of sigma-sigma_h
		for (i1 = 0; i1<num_qp1; i1++)
		{
			lambdas[0] = gauss1[i1][0];
			lambdas[1] = gauss1[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				huzhang_basisDIV(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] += f1(x, y, lambda, mu);
			value[1] += f2(x, y, lambda, mu);
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			errors[1] += 2 * s*gauss1[i1][2] * (value[0] + value[1]);
		}

		// L2 norm of u-u_h
		for (i1 = 0; i1<num_qp2; i1++)
		{
			lambdas[0] = gauss2[i1][0];
			lambdas[1] = gauss2[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*uh->val[j1];
				j1 += elementDOF[1].dof;
				value[1] += phi0*uh->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= u1(x, y, lambda, mu);
			value[1] -= u2(x, y, lambda, mu);
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			errors[3] += 2 * s*gauss2[i1][2] * (value[0] + value[1]);
		}

		// H1 semi-norm of Qhu-u_h
		for (i1 = 0; i1<num_qp3; i1++)
		{
			lambdas[0] = gauss3[i1][0];
			lambdas[1] = gauss3[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			value[3] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF[1].dop, phi1);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi1[0] * (Qhu->val[j1] - uh->val[j1]);
				value[1] += phi1[1] * (Qhu->val[j1] - uh->val[j1]);
				j1 += elementDOF[1].dof;
				value[2] += phi1[0] * (Qhu->val[j1] - uh->val[j1]);
				value[3] += phi1[1] * (Qhu->val[j1] - uh->val[j1]);
			}
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			value[3] = value[3] * value[3];
			errors[5] += 2 * s*gauss3[i1][2] * (value[0] + value[1] + value[2] + value[3]);
		}

		// L2 norm of Qhu-u_h
		for (i1 = 0; i1<num_qp4; i1++)
		{
			lambdas[0] = gauss4[i1][0];
			lambdas[1] = gauss4[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*(Qhu->val[j1]- uh->val[j1]);
				j1 += elementDOF[1].dof;
				value[1] += phi0*(Qhu->val[j1] - uh->val[j1]);
			}
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			errors[6] += 2 * s*gauss4[i1][2] * (value[0] + value[1]);
		}

		// H1 semi-norm of u-uhstar
		for (i1 = 0; i1<num_qp2; i1++)
		{
			lambdas[0] = gauss2[i1][0];
			lambdas[1] = gauss2[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			value[3] = 0;
			for (k1 = 0; k1<elementDOF[2].col; k1++)
			{
				lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF[2].dop, phi1);
				j1 = elementDOF[2].val[k][k1];
				value[0] += phi1[0] * uhstar->val[j1];
				value[1] += phi1[1] * uhstar->val[j1];
				j1 += elementDOF[2].dof;
				value[2] += phi1[0] * uhstar->val[j1];
				value[3] += phi1[1] * uhstar->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= u1_x(x, y, lambda, mu);
			value[1] -= u1_y(x, y, lambda, mu);
			value[2] -= u2_x(x, y, lambda, mu);
			value[3] -= u2_y(x, y, lambda, mu);
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			value[3] = value[3] * value[3];
			errors[7] += 2 * s*gauss2[i1][2] * (value[0] + value[1] + value[2] + value[3]);
		}

		// L2 norm of u-uhstar
		for (i1 = 0; i1<num_qp2; i1++)
		{
			lambdas[0] = gauss2[i1][0];
			lambdas[1] = gauss2[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[2].col; k1++)
			{
				lagrange_basis(lambdas, k1, elementDOF[2].dop, &phi0);
				j1 = elementDOF[2].val[k][k1];
				value[0] += phi0*uhstar->val[j1];
				j1 += elementDOF[2].dof;
				value[1] += phi0*uhstar->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= u1(x, y, lambda, mu);
			value[1] -= u2(x, y, lambda, mu);
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			errors[8] += 2 * s*gauss2[i1][2] * (value[0] + value[1]);
		}
	}

	double elen, C11;
	int element[2], edgeNode[2];
	int count;
	int patchnodes[100];
	int edge;
	// error of lambdah for energy norm
	for (edge = 0; edge<edges->row; edge++)
	{
		edgeNode[0] = edges->val[edge][0];
		edgeNode[1] = edges->val[edge][1];
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];
		C11 = elen;

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

		for (i1 = 0; i1<num_qp11; i1++)
		{
			for (i = 0; i<3; i++)
				phi2[i] = 0;

			lambdas[0] = gauss11[i1][0];
			lambdas[1] = 1 - lambdas[0];
			for (k1 = 0; k1<count; k1++)
			{
				i = patchnodes[k1];
				jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], i, phi1);

				for (j = 0; j<2; j++)
					phi2[j] += uh->val[i] * phi1[j];

				jumpOperatorVector(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, &elementDOF[1], i + elementDOF[1].dof, phi1);
				for (j = 0; j<2; j++)
					phi2[j] += uh->val[i + elementDOF[1].dof] * phi1[j];
			}

			errors[2] += elen*gauss11[i1][1] * C11*(phi2[0] * phi2[0] + phi2[1] * phi2[1]);
		}
	} // e

	errors[4] = errors[0] + errors[1] + errors[2];

	for (i = 0; i<10; i++)
		errors[i] = sqrt(errors[i]);
}


/**
* \fn void getposteriorierrors(double *errors, dvector *sigmah, dvector *uh, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief compute error between numerical solution and exact solution: L2 norm, Energy norm
* \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
* \param *sigmah pointer to numerical solution
* \param *uh pointer to numerical solution
* \param *uhstar pointer to postprocessed displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
*/
void getposteriorierrors(double *errors, dvector *sigmah, dvector *uh, dvector *uhstar, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	dvector Qhf;
	create_dvector(uh->row, &Qhf);
	projPiecewiseLagrangeRHS(&Qhf, elements, nodes, &elementDOF[1], lambda, mu);

	int Nt = elements->row;

	int i, j, k;

	for (i = 0; i<8; i++)
		errors[i] = 0;

	double phi0, phi1[3], phi2[3];
	int k1, i1, j1, l1, l2;
	double x, y, xs[3], ys[3], lambdas[3], s, *eta, *xi, nv[2], tv[2];
	double value[4];
	double **tensorBasis[3];

	int num_qp0 = getNumQuadPoints((elementDOF[0].dop - 2) * 2, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	int num_qp1 = getNumQuadPoints(elementDOF[0].dop * 2, 2); // the number of numerical intergation points
	double gauss1[num_qp1][3];
	init_Gauss(num_qp1, 2, gauss1); // gauss intergation initial

	int num_qp = 49; // the number of numerical intergation points
	double gauss[num_qp][3];
	init_Gauss(num_qp, 2, gauss); // gauss intergation initial

//	int num_qp11 = getNumQuadPoints(elementDOF[0].dop * 2, 1); // the number of numerical intergation points
//	if (num_qp11>5)
//		num_qp11 = 5;
	int num_qp11 = 20;
	double gauss11[num_qp11][2];
	init_Gauss1D(num_qp11, 1, gauss11); // gauss intergation initial

//	int num_qp12 = getNumQuadPoints((elementDOF[0].dop - 1) * 2, 1); // the number of numerical intergation points
//	if (num_qp12>5)
//		num_qp12 = 5;
	int num_qp12 = 20;
	double gauss12[num_qp12][2];
	init_Gauss1D(num_qp12, 1, gauss12); // gauss intergation initial

	for (k = 0; k<Nt; k++)
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
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		// rotrot(Asigmah)
		for (i1 = 0; i1<num_qp0; i1++)
		{
			lambdas[0] = gauss0[i1][0];
			lambdas[1] = gauss0[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			if (elementDOF[0].dop>1)
			{
				for (k1 = 0; k1<elementDOF[0].col; k1++)
				{
					j1 = elementDOF[0].val[k][k1];
					huzhang_basisROTROT(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, &phi0);
					value[0] += phi0 * sigmah->val[j1];
					huzhang_basisLaplaceTrace(lambdas, s, eta, xi, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, &phi0);
					value[1] += phi0 * sigmah->val[j1];
				}
			}
			value[2] = (value[0] - value[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			value[2] = value[2] * value[2];
			errors[0] += 2 * s*gauss0[i1][2] * value[2] * s*s;
		}

		// oscillation: f-Qhf
		for (i1 = 0; i1<num_qp; i1++)
		{
			lambdas[0] = gauss[i1][0];
			lambdas[1] = gauss[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1<elementDOF[1].col; k1++)
			{
				lagrange_basis(lambdas, k1, elementDOF[1].dop, &phi0);
				j1 = elementDOF[1].val[k][k1];
				value[0] += phi0*Qhf.val[j1];
				j1 += elementDOF[1].dof;
				value[1] += phi0*Qhf.val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= f1(x, y, lambda, mu);
			value[1] -= f2(x, y, lambda, mu);
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			errors[2] += 2 * s*gauss[i1][2] * (value[0] + value[1])*s;
		}

		continue;////////////////////////////////////////////////////////////////////////
				 
		// Asigmah - \varepsilon_h(uhstar)
		for (i1 = 0; i1<num_qp1; i1++)
		{
			lambdas[0] = gauss1[i1][0];
			lambdas[1] = gauss1[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			// sigmah
			for (k1 = 0; k1<elementDOF[0].col; k1++)
			{
				huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi1[0] * sigmah->val[j1];
				value[1] += phi1[1] * sigmah->val[j1];
				value[2] += phi1[2] * sigmah->val[j1];
			}
			// Asigmah
			value[3] = value[0] + value[1];
			value[0] = (value[0] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			value[1] = (value[1] - value[3] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			value[2] = value[2] / (2 * mu);
			// Asigmah - \varepsilon_h(uhstar)
			for (k1 = 0; k1<elementDOF[2].col; k1++)
			{
				lagrange_basis1(lambdas, s, eta, xi, k1, elementDOF[2].dop, phi1);
				j1 = elementDOF[2].val[k][k1];
				value[0] -= phi1[0] * uhstar->val[j1];
				value[2] -= phi1[1] * uhstar->val[j1] / 2;
				j1 += elementDOF[2].dof;
				value[1] -= phi1[1] * uhstar->val[j1];
				value[2] -= phi1[0] * uhstar->val[j1] / 2;
			}
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			errors[1] += 2 * s*gauss1[i1][2] * (value[0] + value[1] + 2 * value[2]);
		}
	} // k

	int element[2];
	int count, istart;
	int *index;
	index = (int*)calloc(elementDOF[0].dof, sizeof(int));
	for (i = 0; i<elementDOF[0].dof; i++)
		index[i] = -1;

	int patchnodes[200];
	int edge;
	double elen;

	for (edge = 0; edge<edges->row; edge++)
	{
		element[0] = edges->val[edge][2];
		element[1] = edges->val[edge][3];
		elen = edges->length[edge];
		for (i = 0; i<2; i++)
		{
			j = edges->val[edge][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}

		istart = -2;
		count = 0;
		for (i = 0; i<elementDOF[0].col; i++)
		{
			j = elementDOF[0].val[element[0]][i];
			patchnodes[count] = j;
			count++;
			index[j] = istart;
			istart = j;
		}

		if (element[1] != -1)
		{
			for (i = 0; i<elementDOF[0].col; i++)
			{
				j = elementDOF[0].val[element[1]][i];
				if (index[j] == -1)
				{
					patchnodes[count] = j;
					count++;
					index[j] = istart;
					istart = j;

				}
			}
		}

		for (j = 0; j<count; j++)
		{
			j1 = istart;
			istart = index[j1];
			index[j1] = -1;
		}

		// M_{tt}(Asigmah)
		if (edges->bdFlag[edge] == 0) // interior edge
		{
			for (i1 = 0; i1<num_qp11; i1++)
			{
				value[0] = 0;

				lambdas[0] = gauss11[i1][0];
				lambdas[1] = 1 - lambdas[0];
				for (k1 = 0; k1<count; k1++)
				{
					i = patchnodes[k1];
					jumpOperatorATensorTangent2(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, nodes, &elementDOF[0], i, lambda, mu, &phi0);
					value[0] += sigmah->val[i] * phi0;
				}

				errors[3] += elen*gauss11[i1][1] * value[0] * value[0] * elen;
			}
		}
		else if (edges->bdFlag[edge] >0 && edges->bdFlag[edge]<5) // Dirichlet boundary edge
		{
			xi = elements->xi[element[0]];
			eta = elements->eta[element[0]];
			for (i = 0; i < 3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			nv[0] = -eta[i] / elen;
			nv[1] = xi[i] / elen;
			tv[0] = -nv[1];
			tv[1] = nv[0];

			for (i1 = 0; i1<num_qp11; i1++)
			{
				value[0] = 0;

				lambdas[0] = gauss11[i1][0];
				lambdas[1] = 1 - lambdas[0];
				for (k1 = 0; k1<count; k1++)
				{
					i = patchnodes[k1];
					jumpOperatorATensorTangent2Dirichlet(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, nodes, &elementDOF[0], i, lambda, mu, &phi0);
					value[0] += sigmah->val[i] * phi0;
				}

				x = xs[0] * lambdas[0] + xs[1] * lambdas[1];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1];
				value[0] -= ut_t(x, y, lambda, mu, edges->bdFlag[edge]);
//				value[0] -= (u1_x(x, y, lambda, mu)*tv[0] + u1_y(x, y, lambda, mu)*tv[1]) * tv[0];
//				value[0] -= (u2_x(x, y, lambda, mu)*tv[0] + u2_y(x, y, lambda, mu)*tv[1]) * tv[1];

//				if(edges->bdFlag[edge] ==4)
				errors[3] += elen*gauss11[i1][1] * value[0] * value[0] * elen;
			}
		}

		// rot(Asigmah) t - \partial_t(M_{nt}(Asigmah))
		if (edges->bdFlag[edge] == 0) // interior edge
		{
			for (i1 = 0; i1<num_qp12; i1++)
			{
				value[0] = 0;

				lambdas[0] = gauss12[i1][0];
				lambdas[1] = 1 - lambdas[0];
				for (k1 = 0; k1<count; k1++)
				{
					i = patchnodes[k1];
					jumpOperatorRotATensorTangentPt(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, nodes, &elementDOF[0], i, lambda, mu, &phi0);
					value[0] += sigmah->val[i] * phi0;
				}

				errors[4] += elen*gauss12[i1][1] * value[0] * value[0] * elen * elen * elen;
			}
		}
		else if (edges->bdFlag[edge] >0 && edges->bdFlag[edge]<5) // Dirichlet boundary edge
		{
			xi = elements->xi[element[0]];
			eta = elements->eta[element[0]];
			for (i = 0; i < 3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			nv[0] = -eta[i] / elen;
			nv[1] = xi[i] / elen;
			tv[0] = -nv[1];
			tv[1] = nv[0];

			for (i1 = 0; i1<num_qp12; i1++)
			{
				value[0] = 0;

				lambdas[0] = gauss12[i1][0];
				lambdas[1] = 1 - lambdas[0];
				for (k1 = 0; k1<count; k1++)
				{
					i = patchnodes[k1];
					jumpOperatorRotATensorTangentPt(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, nodes, &elementDOF[0], i, lambda, mu, &phi0);
					value[0] += sigmah->val[i] * phi0;
				}

				x = xs[0] * lambdas[0] + xs[1] * lambdas[1];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1];
				value[0] += un_tt(x, y, lambda, mu, edges->bdFlag[edge]);
//				value[0] += (u1_xx(x, y, lambda, mu)*tv[0] * tv[0] + u1_yy(x, y, lambda, mu)*tv[1] * tv[1] + (u1_xy(x, y, lambda, mu) + u1_yx(x, y, lambda, mu))*tv[0] * tv[1]) * nv[0];
//				value[0] += (u2_xx(x, y, lambda, mu)*tv[0] * tv[0] + u2_yy(x, y, lambda, mu)*tv[1] * tv[1] + (u2_xy(x, y, lambda, mu) + u2_yx(x, y, lambda, mu))*tv[0] * tv[1]) * nv[1];

//				if(edges->bdFlag[edge] ==4)
				errors[4] += elen*gauss12[i1][1] * value[0] * value[0] * elen * elen * elen;
			}
		}

		continue;////////////////////////////////////////////////////////////////////////
		// Neumann Oscillation
		if (edges->bdFlag[edge] == 5 || edges->bdFlag[edge] == 6) // Neumann boundary edge
		{
			xi = elements->xi[element[0]];
			eta = elements->eta[element[0]];
			for (i = 0; i < 3; i++)
			{
				if (elementEdge->val[element[0]][i] == edge)
					break;
			}
			nv[0] = -eta[i] / elen;
			nv[1] = xi[i] / elen;
			tv[0] = -nv[1];
			tv[1] = nv[0];

			for (i1 = 0; i1<num_qp11; i1++)
			{
				value[0] = 0;
				value[1] = 0;

				lambdas[0] = gauss11[i1][0];
				lambdas[1] = 1 - lambdas[0];
				for (k1 = 0; k1<count; k1++)
				{
					i = patchnodes[k1];
					NeumannStress(lambdas[0], lambdas[1], edge, elements, elementEdge, edges, nodes, &elementDOF[0], i, lambda, mu, phi1);
					value[0] += sigmah->val[i] * phi1[0];  // sigmah is wrong, it should be IhgN
					value[1] += sigmah->val[i] * phi1[1];  // sigmah is wrong, it should be IhgN
				}

				x = xs[0] * lambdas[0] + xs[1] * lambdas[1];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1];
				value[0] -= (sigma11(x, y, lambda, mu)*nv[0] + sigma12(x, y, lambda, mu)*nv[1]);
				value[1] -= (sigma21(x, y, lambda, mu)*nv[0] + sigma22(x, y, lambda, mu)*nv[1]);

				errors[5] += elen*gauss11[i1][1] * value[0] * value[0] * elen;
			}
		}
	} // edge
	free(index);

	free_dvector(&Qhf);

	errors[6] = errors[0] + errors[1] + errors[3] + errors[4];
	errors[7] = errors[6] + errors[2] + errors[5];
//	errors[5] = errors[0] + errors[1] + errors[3] + errors[4];
//	errors[6] = errors[5] + errors[2];

	for (i = 0; i<8; i++)
		errors[i] = sqrt(errors[i]);
}

/**
* \fn void geterrors_kirchhoffpalte_morley(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief compute error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
* \param *errors pointer to error between numerical solution and exact solution: L2 norm, H1 norm, Energy norm
* \param *sigmah pointer to numerical solution
* \param *uh pointer to numerical solution
* \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
* \param *uhstar pointer to postprocessed displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
*/
void geterrors_kirchhoffpalte_morley(double *errors, dvector *uh, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int Nt = elements->row;

	int i, j, k, l;

	for (i = 0; i < 3; i++)
		errors[i] = 0;

	double phi[3];
	int k1, i1, j1, l1, l2;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	double value[4];
	double **tensorBasis[3];

	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	int num_qp1 = 49; // the number of numerical intergation points
	double gauss1[num_qp1][3];
	init_Gauss(num_qp1, 2, gauss1); // gauss intergation initial

	int num_qp2 = 49; // the number of numerical intergation points
	double gauss2[num_qp2][3];
	init_Gauss(num_qp2, 2, gauss2); // gauss intergation initial

	for (k = 0; k < Nt; k++)
	{
		for (i = 0; i < 3; i++)
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

		// L2 norm of u-u_h
		for (i1 = 0; i1 < num_qp0; i1++)
		{
			lambdas[0] = gauss0[i1][0];
			lambdas[1] = gauss0[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			for (k1 = 0; k1 < elementDOF[0].col; k1++)
			{
				morley_basis(lambdas, s, elen, eta, xi, nv, nve, k1, phi);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi[0] * uh->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= u(x, y);
			value[0] = value[0] * value[0];
			errors[0] += 2 * s*gauss0[i1][2] * value[0];
		}

		// H1 norm of u-u_h
		for (i1 = 0; i1 < num_qp1; i1++)
		{
			lambdas[0] = gauss1[i1][0];
			lambdas[1] = gauss1[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1 < elementDOF[0].col; k1++)
			{
				morley_basis1(lambdas, s, elen, eta, xi, nv, nve, k1, phi);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi[0] * uh->val[j1];
				value[1] += phi[1] * uh->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= u_x(x, y);
			value[1] -= u_y(x, y);
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			errors[1] += 2 * s*gauss1[i1][2] * (value[0] + value[1]);
		}

		// H2 norm of u-u_h
		for (i1 = 0; i1 < num_qp2; i1++)
		{
			lambdas[0] = gauss2[i1][0];
			lambdas[1] = gauss2[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			for (k1 = 0; k1 < elementDOF[0].col; k1++)
			{
				morley_basis2(lambdas, s, elen, eta, xi, nv, nve, k1, phi);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi[0] * uh->val[j1];
				value[1] += phi[1] * uh->val[j1];
				value[2] += phi[2] * uh->val[j1];
			}
			x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
			y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
			value[0] -= u_xx(x, y);
			value[1] -= u_yy(x, y);
			value[2] -= u_xy(x, y);
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			errors[2] += 2 * s*gauss1[i1][2] * (value[0] + value[1] + 2 * value[2]);
		}
	}
		errors[1] += errors[0];
		errors[2] += errors[1];

	for (i = 0; i<3; i++)
		errors[i] = sqrt(errors[i]);
}

/**
* \fn void geterrors2discr_kirchhoffpalte_morley(double *errors, dvector *uh1, dvector *uh2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
* \brief compute error between two numerical solutions: L2 norm, H1 norm, Energy norm
* \param *errors pointer to error between two numerical solutions: L2 norm, H1 norm, Energy norm
* \param *sigmah pointer to numerical solution
* \param *uh pointer to numerical solution
* \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
* \param *uhstar pointer to postprocessed displacement
* \param *elements pointer to the structure of the triangulation
* \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
* \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
the fourth column stores -1 if the edge is on boundary
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
*/
void geterrors2discr_kirchhoffpalte_morley(double *errors, dvector *uh1, dvector *uh2, ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF)
{
	int Nt = elements->row;

	int i, j, k, l;

	for (i = 0; i < 3; i++)
		errors[i] = 0;

	double phi[3];
	int k1, i1, j1, l1, l2;
	double x, y, lambdas[3], xs[3], ys[3], s, *eta, *xi, *elen, **nv, *nve[3];
	double value[4];
	double **tensorBasis[3];

	int num_qp0 = getNumQuadPoints(4, 2); // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	int num_qp1 = getNumQuadPoints(2, 2); // the number of numerical intergation points
	double gauss1[num_qp1][3];
	init_Gauss(num_qp1, 2, gauss1); // gauss intergation initial

	int num_qp2 = getNumQuadPoints(0, 2); // the number of numerical intergation points
	double gauss2[num_qp2][3];
	init_Gauss(num_qp2, 2, gauss2); // gauss intergation initial

	for (k = 0; k < Nt; k++)
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

		// L2 norm of u-u_h
		for (i1 = 0; i1 < num_qp0; i1++)
		{
			lambdas[0] = gauss0[i1][0];
			lambdas[1] = gauss0[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			for (k1 = 0; k1 < elementDOF[0].col; k1++)
			{
				morley_basis(lambdas, s, elen, eta, xi, nv, nve, k1, phi);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi[0] * (uh1->val[j1]- uh2->val[j1]);
			}
			value[0] = value[0] * value[0];
			errors[0] += 2 * s*gauss0[i1][2] * value[0];
		}

		// H1 norm of u-u_h
		for (i1 = 0; i1 < num_qp1; i1++)
		{
			lambdas[0] = gauss1[i1][0];
			lambdas[1] = gauss1[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			for (k1 = 0; k1 < elementDOF[0].col; k1++)
			{
				morley_basis1(lambdas, s, elen, eta, xi, nv, nve, k1, phi);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi[0] * (uh1->val[j1] - uh2->val[j1]);
				value[1] += phi[1] * (uh1->val[j1] - uh2->val[j1]);
			}
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			errors[1] += 2 * s*gauss1[i1][2] * (value[0] + value[1]);
		}

		// H2 norm of u-u_h
		for (i1 = 0; i1 < num_qp2; i1++)
		{
			lambdas[0] = gauss2[i1][0];
			lambdas[1] = gauss2[i1][1];
			lambdas[2] = 1 - lambdas[0] - lambdas[1];
			value[0] = 0;
			value[1] = 0;
			value[2] = 0;
			for (k1 = 0; k1 < elementDOF[0].col; k1++)
			{
				morley_basis2(lambdas, s, elen, eta, xi, nv, nve, k1, phi);
				j1 = elementDOF[0].val[k][k1];
				value[0] += phi[0] * (uh1->val[j1] - uh2->val[j1]);
				value[1] += phi[1] * (uh1->val[j1] - uh2->val[j1]);
				value[2] += phi[2] * (uh1->val[j1] - uh2->val[j1]);
			}
			value[0] = value[0] * value[0];
			value[1] = value[1] * value[1];
			value[2] = value[2] * value[2];
			errors[2] += 2 * s*gauss1[i1][2] * (value[0] + value[1] + 2 * value[2]);
		}
	}
	errors[1] += errors[0];
	errors[2] += errors[1];

	for (i = 0; i<3; i++)
		errors[i] = sqrt(errors[i]);
}
