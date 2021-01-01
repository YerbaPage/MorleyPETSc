/*
 *  estimate.c
 *
 *  Created by Xuehai Huang on 08/26/2010.
 *  Modified by Xuehai Huang on 10/01/2016.
 *  Copyright 2016 WZU. All rights reserved.
 *
 */

/*! \file estimate.c
 *  \brief Error estimator
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "matvec.h"


/**
 * \fn double estimate(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *sigmah, dvector *uh, dvector *uhstar, double lambda, double mu, dvector *estimator)
 * \brief compute the error estimator for adaptive finite element method
 * \param *elements pointer to the structure of the triangulation
 * \param *elementEdge pointer to relation between tirangles and edges: each row stores 3 edges index
 * \param *edges pointer to edges: the first two columns store the two vertice, the third and fourth columns store the affiliated elements
								   the fourth column stores -1 if the edge is on boundary
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *sigmah pointer to numerical solution
 * \param *uh pointer to numerical solution
 * \param *uhstar pointer to postprocessed displacement
 * \param lambda Lame constant
 * \param mu Lame constant or Poisson ratio of plate
 * \param *estimator pointer to the error estimator
 * \return total posterior error
 */
double estimate(ELEMENT *elements, idenmat *elementEdge, EDGE *edges, dennode *nodes, ELEMENT_DOF *elementDOF, dvector *sigmah, dvector *uh, dvector *uhstar, double lambda, double mu, dvector *estimator)
{
	dvector etaElem, etaEdge;
	int Nt=elements->row;
	int Ne=edges->row;

	create_dvector(Nt, estimator);
	create_dvector(Nt, &etaElem);
	create_dvector(Ne, &etaEdge);

	dvector Qhf;
	create_dvector(uh->row, &Qhf);
	projPiecewiseLagrangeRHS(&Qhf, elements, nodes, &elementDOF[1], lambda, mu);

	int i, j, k;
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
		s=elements->vol[k];
		xi=elements->xi[k];
		eta=elements->eta[k];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		etaElem.val[k]=0;

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
			etaElem.val[k] += 2 * s*gauss0[i1][2] * value[2] * s*s;
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
			etaElem.val[k] += 2 * s*gauss[i1][2] * (value[0] + value[1])*s;
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
			etaElem.val[k] += 2 * s*gauss1[i1][2] * (value[0] + value[1] + 2 * value[2]);
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

				etaEdge.val[edge] += elen*gauss11[i1][1] * value[0] * value[0] * elen;
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

				etaEdge.val[edge] += elen*gauss11[i1][1] * value[0] * value[0] * elen;
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

				etaEdge.val[edge] += elen*gauss12[i1][1] * value[0] * value[0] * elen * elen * elen;
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

				etaEdge.val[edge] += elen*gauss12[i1][1] * value[0] * value[0] * elen * elen * elen;
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
					value[0] += sigmah->val[i] * phi1[0]; // sigmah is wrong, it should be IhgN 
					value[1] += sigmah->val[i] * phi1[1]; // sigmah is wrong, it should be IhgN
				}

				x = xs[0] * lambdas[0] + xs[1] * lambdas[1];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1];
				value[0] -= (sigma11(x, y, lambda, mu)*nv[0] + sigma12(x, y, lambda, mu)*nv[1]);
				value[1] -= (sigma21(x, y, lambda, mu)*nv[0] + sigma22(x, y, lambda, mu)*nv[1]);

				etaEdge.val[edge] += elen*gauss11[i1][1] * value[0] * value[0] * elen;
			}
		}
	} // edge
	free(index);

	double totalError=0;
	for (k = 0; k<Nt; k++)
		totalError += etaElem.val[k];

	for(edge=0;edge<Ne;edge++)
	{
		if (edges->val[edge][3] > -1)
		{
			totalError += etaEdge.val[edge];
			etaEdge.val[edge] *= 0.5;
		}
		else
		{
			etaEdge.val[edge] *= 0.5;
			totalError += etaEdge.val[edge];
		}
	}
	
	/////////////////////////////////////
/*	printf("eta Element\n");
	for(i=0;i<Nt;i++)
	{
		printf("%e,",etaElem.val[i]);
	}
	printf("\n");
	printf("eta Edge\n");
	for(i=0;i<Ne;i++)
	{
		printf("%e,",etaEdge.val[i]);
	}
	printf("\n");*/

	for (k = 0; k<Nt; k++)
	{
		estimator->val[k] = etaElem.val[k];
		for (i = 0; i<3; i++)
		{
			edge = elementEdge->val[k][i];
			estimator->val[k] += etaEdge.val[edge];
		}
	}

	free_dvector(&etaElem);
	free_dvector(&etaEdge);

	//double errors[4];
	//geterrors(errors, sigmah, uh, elements, elementEdge, edges, nodes, elementDOF, nu);
	printf("Posterior Error=%e\n", sqrt(totalError));
//	printf("Enery norm of error=%e, %e\n", errors[0], errors[0] + errors[3]);

	return totalError;
}
