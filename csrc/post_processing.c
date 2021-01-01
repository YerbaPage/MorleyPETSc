/*
 *  post_processing.c
 *  FEM
 *
 *  Created by Xuehai Huang on 03/13/12.
 *  Copyright 2012 WZU. All rights reserved.
 *
 */

/*! \file post_processing.c
 *  \brief Post-processing: Project the discontinuous Galerkin finite element solution to Raviart-Thomas element space
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "precond.h"
#include "matvec.h"


 /**
 * \fn void projPiecewiseLagrangeDisplacement(dvector *Qhu, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
 * \brief L2 projection of the exact solution to the piecewise Lagrangian element space
 * \param *Qhu pointer to the L2 projection of u to the piecewise Lagrangian element space
 * \param *elements pointer to the structure of the triangulation
 * \param *nodes pointer to the nodes location of the triangulation
 * \param *elementDOF pointer to relation between elements and DOFs
 * \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
 */
void projPiecewiseLagrangeDisplacement(dvector *Qhu, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int dop = elementDOF->dop;
	int dof = elementDOF->dof;
	int nt = elementDOF->row;
	ddenmat A;
	dvector b[2], w;
	create_dden_matrix(elementDOF->col, elementDOF->col, &A);
	create_dvector(A.col, &b[0]);
	create_dvector(A.col, &b[1]);
	create_dvector(A.col, &w);

	double phi[2];
	int i, j, k, i1;
	double x, y, lambdas[3], xs[3], ys[3], s;

	int num_qp = getNumQuadPoints(dop * 2, 2); // the number of numerical intergation points
	double gauss[num_qp][3];
	init_Gauss(num_qp, 2, gauss); // gauss intergation initial

	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k < nt; k++)
	{
		s = elements->vol[k];
		for (i = 0; i < A.row; i++)
		{
			for (j = 0; j < A.col; j++)
			{
				A.val[i][j] = 0;

				for (i1 = 0; i1 < num_qp; i1++)
				{
					lambdas[0] = gauss[i1][0];
					lambdas[1] = gauss[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, i, dop, &phi[0]);
					lagrange_basis(lambdas, j, dop, &phi[1]);
					A.val[i][j] += 2 * s*gauss[i1][2] * phi[0] * phi[1];
				}
			}
		}

		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}

		for (i = 0; i<A.row; i++)
		{
			b[0].val[i] = 0;
			b[1].val[i] = 0;
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, i, dop, &phi[0]);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b[0].val[i] += 2 * s*gauss0[i1][2] * u1(x, y, lambda, mu)*phi[0];
				b[1].val[i] += 2 * s*gauss0[i1][2] * u2(x, y, lambda, mu)*phi[0];
			} // i1
		} // k1

		init_dvector(&w, 0.0);
		den_pcg(&A, &b[0], &w, 10000, 1e-14, NULL, 0);
		for (i = 0; i < A.row; i++)
			Qhu->val[elementDOF->val[k][i]] = w.val[i];
		init_dvector(&w, 0.0);
		den_pcg(&A, &b[1], &w, 10000, 1e-14, NULL, 0);
		for (i = 0; i < A.row; i++)
			Qhu->val[elementDOF->val[k][i] + dof] = w.val[i];
	}

	free_dden_matrix(&A);
	free_dvector(&b[0]);
	free_dvector(&b[1]);
	free_dvector(&w);
}

/**
* \fn void projPiecewiseLagrangeRHS(dvector *Qhf, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief L2 projection of the right hand side term f to the piecewise Lagrangian element space
* \param *Qhf pointer to the L2 projection of f to the piecewise Lagrangian element space
* \param *elements pointer to the structure of the triangulation
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void projPiecewiseLagrangeRHS(dvector *Qhf, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int dop = elementDOF->dop;
	int dof = elementDOF->dof;
	int nt = elementDOF->row;
	ddenmat A;
	dvector b[2], w;
	create_dden_matrix(elementDOF->col, elementDOF->col, &A);
	create_dvector(A.col, &b[0]);
	create_dvector(A.col, &b[1]);
	create_dvector(A.col, &w);

	double phi[2];
	int i, j, k, i1;
	double x, y, lambdas[3], xs[3], ys[3], s;

	int num_qp = getNumQuadPoints(dop * 2, 2); // the number of numerical intergation points
	double gauss[num_qp][3];
	init_Gauss(num_qp, 2, gauss); // gauss intergation initial

	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k < nt; k++)
	{
		s = elements->vol[k];
		for (i = 0; i < A.row; i++)
		{
			for (j = 0; j < A.col; j++)
			{
				A.val[i][j] = 0;

				for (i1 = 0; i1 < num_qp; i1++)
				{
					lambdas[0] = gauss[i1][0];
					lambdas[1] = gauss[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, i, dop, &phi[0]);
					lagrange_basis(lambdas, j, dop, &phi[1]);
					A.val[i][j] += 2 * s*gauss[i1][2] * phi[0] * phi[1];
				}
			}
		}

		for (i = 0; i<3; i++)
		{
			j = elements->val[k][i];
			xs[i] = nodes->val[j][0];
			ys[i] = nodes->val[j][1];
		}

		for (i = 0; i<A.row; i++)
		{
			b[0].val[i] = 0;
			b[1].val[i] = 0;
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				lagrange_basis(lambdas, i, dop, &phi[0]);
				x = xs[0] * lambdas[0] + xs[1] * lambdas[1] + xs[2] * lambdas[2];
				y = ys[0] * lambdas[0] + ys[1] * lambdas[1] + ys[2] * lambdas[2];
				b[0].val[i] += 2 * s*gauss0[i1][2] * f1(x, y, lambda, mu)*phi[0];
				b[1].val[i] += 2 * s*gauss0[i1][2] * f2(x, y, lambda, mu)*phi[0];
			} // i1
		} // k1

		init_dvector(&w, 0.0);
		den_pcg(&A, &b[0], &w, 10000, 1e-14, NULL, 0);
		for (i = 0; i < A.row; i++)
			Qhf->val[elementDOF->val[k][i]] = w.val[i];
		init_dvector(&w, 0.0);
		den_pcg(&A, &b[1], &w, 10000, 1e-14, NULL, 0);
		for (i = 0; i < A.row; i++)
			Qhf->val[elementDOF->val[k][i] + dof] = w.val[i];
	}

	free_dden_matrix(&A);
	free_dvector(&b[0]);
	free_dvector(&b[1]);
	free_dvector(&w);
}

/**
* \fn void postprocess2newDisplacement(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
* \brief L2 projection of the exact solution to the piecewise Lagrangian element space
* \param *uhstar pointer to postprocessed displacement
* \param *sigmah pointer to numerical solution
* \param *uh pointer to numerical solution
* \param *elements pointer to the structure of the triangulation
* \param *nodes pointer to the nodes location of the triangulation
* \param *elementDOF pointer to relation between elements and DOFs
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return void
*/
void postprocess2newDisplacement(dvector *uhstar, dvector *sigmah, dvector *uh, ELEMENT *elements, dennode *nodes, ELEMENT_DOF *elementDOF, double lambda, double mu)
{
	int dop1 = elementDOF[1].dop;
	int dop2 = elementDOF[2].dop;
	int ldof1 = elementDOF[1].col;
	int ldof2 = elementDOF[2].col;
	int nt = elements->row;

	ddenmat A;
	dvector b, x;
	create_dden_matrix((ldof1 + ldof2) * 2, (ldof1 + ldof2) * 2, &A);
	create_dvector(A.row, &b);
	create_dvector(A.row, &x);

	create_dvector(elementDOF[2].dof * 2, uhstar);

	double phi[2], phi1[3], phi2[3];
	int i, j, k, i1, j1, k1;
	double lambdas[3], *eta, *xi, s;
	double value[4];
	double **tensorBasis[3];

	int num_qp = getNumQuadPoints(dop1 + dop2, 2); // the number of numerical intergation points
	double gauss[num_qp][3];
	init_Gauss(num_qp, 2, gauss); // gauss intergation initial

	int num_qp0 = 49; // the number of numerical intergation points
	double gauss0[num_qp0][3];
	init_Gauss(num_qp0, 2, gauss0); // gauss intergation initial

	for (k = 0; k < nt; k++)
	{
		// set parameters
		s = elements->vol[k];
		xi = elements->xi[k];
		eta = elements->eta[k];
		for (i = 0; i < 3; i++)
			tensorBasis[i] = nodes->tensorBasis[elements->val[k][i]];
		// end set parameters

		// initial
		for (i = 0; i<A.row; i++)
		{
			for (j = 0; j<A.col; j++)
				A.val[i][j] = 0;
			b.val[i] = 0;
		}

		// A
		for (i = 0; i < ldof2; i++)
		{
			// assemble top-left block of A
			for (j = 0; j < ldof2; j++)
			{
				for (i1 = 0; i1 < num_qp; i1++)
				{
					lambdas[0] = gauss[i1][0];
					lambdas[1] = gauss[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis1(lambdas, s, eta, xi, i, dop2, phi1);
					lagrange_basis1(lambdas, s, eta, xi, j, dop2, phi2);
					A.val[i][j] += 2 * s*gauss[i1][2] * (phi1[0] * phi2[0] + phi1[1] * phi2[1] / 2);
					A.val[i][j + ldof2] += 2 * s*gauss[i1][2] * (phi1[1] * phi2[0] / 2);
					A.val[i + ldof2][j] += 2 * s*gauss[i1][2] * (phi1[0] * phi2[1] / 2);
					A.val[i + ldof2][j + ldof2] += 2 * s*gauss[i1][2] * (phi1[0] * phi2[0] / 2 + phi1[1] * phi2[1]);
				}
			}

			for (j = 0; j < ldof1; j++)
			{
				// assemble top-right block of A
				for (i1 = 0; i1 < num_qp; i1++)
				{
					lambdas[0] = gauss[i1][0];
					lambdas[1] = gauss[i1][1];
					lambdas[2] = 1 - lambdas[0] - lambdas[1];
					lagrange_basis(lambdas, i, dop2, &phi[0]);
					lagrange_basis(lambdas, j, dop1, &phi[1]);
					A.val[i][j + ldof2 * 2] += 2 * s*gauss[i1][2] * phi[0] * phi[1];
					A.val[i + ldof2][j + ldof1 + ldof2 * 2] += 2 * s*gauss[i1][2] * phi[0] * phi[1];
				}
				// assemble bottom-left block of A
				A.val[j + ldof2 * 2][i] = A.val[i][j + ldof2 * 2];
				A.val[j + ldof1 + ldof2 * 2][i + ldof2] = A.val[i + ldof2][j + ldof1 + ldof2 * 2];
			}
		}

		// the first part of b
		for (i = 0; i<ldof2; i++)
		{
			for (i1 = 0; i1<num_qp0; i1++)
			{
				lambdas[0] = gauss0[i1][0];
				lambdas[1] = gauss0[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				value[0] = 0;
				value[1] = 0;
				value[2] = 0;
				lagrange_basis1(lambdas, s, eta, xi, i, dop2, phi2);
				for (k1 = 0; k1<elementDOF[0].col; k1++)
				{
					huzhang_basis(lambdas, elements->nvector[k], elements->tvector[k], tensorBasis, k1, elementDOF[0].dop, phi1);
					j1 = elementDOF[0].val[k][k1];
					value[0] += phi1[0] * sigmah->val[j1];
					value[1] += phi1[1] * sigmah->val[j1];
					value[2] += phi1[2] * sigmah->val[j1];
				}
				b.val[i] += 2 * s*gauss0[i1][2] * (value[0] * phi2[0] + value[2] * phi2[1] - (value[0] + value[1]) * phi2[0] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
				b.val[i + ldof2] += 2 * s*gauss0[i1][2] * (value[1] * phi2[1] + value[2] * phi2[0] - (value[0] + value[1]) * phi2[1] * lambda / (2 * lambda + 2 * mu)) / (2 * mu);
			} // i1
		} // i

		  // the second part of b
		for (i = 0; i<ldof1; i++)
		{
			for (i1 = 0; i1<num_qp; i1++)
			{
				lambdas[0] = gauss[i1][0];
				lambdas[1] = gauss[i1][1];
				lambdas[2] = 1 - lambdas[0] - lambdas[1];
				value[0] = 0;
				value[1] = 0;
				lagrange_basis(lambdas, i, dop1, &phi[1]);
				for (k1 = 0; k1<ldof1; k1++)
				{
					lagrange_basis(lambdas, k1, dop1, &phi[0]);
					j1 = elementDOF[1].val[k][k1];
					value[0] += phi[0] * uh->val[j1];
					j1 += elementDOF[1].dof;
					value[1] += phi[0] * uh->val[j1];
				}
				b.val[i + ldof2 * 2] += 2 * s*gauss[i1][2] * value[0] * phi[1];
				b.val[i + ldof1 + ldof2 * 2] += 2 * s*gauss[i1][2] * value[1] * phi[1];
			} // i1
		} // i

		fgmres_den(&A, &b, &x, 80, 15000, 1e-13, NULL, 0);

		for (i = 0; i < ldof2; i++)
		{
			uhstar->val[elementDOF[2].val[k][i]] = x.val[i];
			uhstar->val[elementDOF[2].val[k][i] + elementDOF[2].dof] = x.val[i + ldof2];
		}
	}

	free_dden_matrix(&A);
	free_dvector(&b);
	free_dvector(&x);
}
