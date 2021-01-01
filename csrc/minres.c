/*
 *  minres.c
 *
 *  Created by Xuehai Huang on 10/22/2009.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file minres.c
 *  \brief Preconditioned Minimal Residual Method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "precond.h"
#include "matvec.h"

 /**
 * \fn int minres(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, precond *pre, int print_level)
 *	 \brief A preconditioned minimal residual (MINRES) method for solving Au=b
 *	 \param *A	 pointer to the coefficient matrix
 *	 \param *b	 pointer to the dvector of right hand side
 *	 \param *x	 pointer to the dvector of DOFs
 *	 \param MaxIt integer, maximal number of iterations
 *	 \param tol double float, the tolerance for stopage
 *	 \param *pre pointer to the structure of precondition (precond)
 *	 \param print_level how much information to print out
 *	 \return the number of iterations
 */
int minres(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, precond *pre, int print_level)
{
	double eps = DBL_EPSILON;
	int i, iter = 0, n;
	int istop = 0;
	double Anorm = 0.0, Acond = 0.0, Arnorm = 0.0;
	double rnorm = 0.0, ynorm = 0.0;
	int done = 0;

	n = b[0].row;

	// Step 1
	/*
	* Set up y and v for the first Lanczos vector v1.
	* y = beta1 P' v1, where P = C^(-1).
	* v is really P'v1
	*/

	dvector y, r1;
	create_dvector(n, &y);
	create_dvector(n, &r1);

	// r = b-A*x
	copy_dvector(b, &r1);
	sparse_mv(-1.0, A, x->val, r1.val);

	// y = B*r1
	if (pre == NULL)
		copy_dvector(&r1, &y); /* No preconditioner, B=I */
	else
		pre->fct_dvec(&r1, &y, pre->data); /* Preconditioning */

	double beta1 = dot_dvector(&r1, &y);


	// Test for an indefined preconditioner
	// If b = 0 exactly stop with x = x0.
	if (beta1 < 0.0)
	{
		istop = 9;
		done = 1;
		// show = true;
	}
	else
	{
		if (beta1 < eps)
		{
			done = 1;
			// show = true;
		}
		else
			beta1 = sqrt(beta1); // Normalize y to get v1 later
	}

	// STEP 2
	/* Initialize other quantities */
	double oldb = 0.0, beta = beta1, dbar = 0.0, epsln = 0.0, oldeps = 0.0;
	double qrnorm = beta1, phi = 0.0, phibar = beta1, rhs1 = beta1;
	double rhs2 = 0.0, tnorm2 = 0.0, ynorm2 = 0.0;
	double cs = -1.0, sn = 0.0;
	double gmax = 0.0, gmin = DBL_MAX;
	double alpha = 0.0, gamma = 0.0;
	double delta = 0.0, gbar = 0.0;
	double z = 0.0;
	double root;

	dvector w, w1, w2, r2, v;
	create_dvector(n, &w);
	create_dvector(n, &w1);
	create_dvector(n, &w2);
	create_dvector(n, &r2);
	create_dvector(n, &v);

	copy_dvector(&r1, &r2);

	/* Main Iteration */
	if (!done)
	{
		for (iter = 0; iter < MaxIt; ++iter)
		{
			// STEP 3
			/*
			-----------------------------------------------------------------
			Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
			The general iteration is similar to the case k = 1 with v0 = 0:
			p1      = A * v1  -  beta1 * v0,
			alpha1  = v1'p1,
			q2      = p2  -  alpha1 * v1,
			beta2^2 = q2'q2,
			v2      = (1/beta2) q2.

			Again, y = betak P vk,  where  P = C**(-1).
			.... more description needed.
			-----------------------------------------------------------------
			*/
			axy_dvector(1. / beta, &y, &v); // Normalize previous vector (in y)

											// y = A*v-r1*beta/oldb
			init_dvector(&y, 0);
			sparse_mv(1.0, A, v.val, y.val);
			if (iter)
				axpy_dvector(-beta / oldb, &r1, &y);

			alpha = dot_dvector(&v, &y);   // alphak
										   // y = y-r2*alpha/beta
			axpy_dvector(-alpha / beta, &r2, &y);
			copy_dvector(&r2, &r1);
			copy_dvector(&y, &r2);

			// y = B*r2
			if (pre == NULL)
				copy_dvector(&r2, &y); /* No preconditioner, B=I */
			else
				pre->fct_dvec(&r2, &y, pre->data); /* Preconditioning */

			oldb = beta; //oldb = betak
			beta = dot_dvector(&r2, &y);
			if (beta < 0)
			{
				istop = 9;
				break;
			}
			beta = sqrt(beta);

			tnorm2 += alpha*alpha + oldb*oldb + beta*beta;

			if (iter == 0)    //Initialize a few things
			{
				if (beta / beta1 < 10.0*eps)
					istop = 10;
			}

			// Apply previous rotation Q_{k-1} to get
			// [delta_k epsln_{k+1}] = [cs sn]  [dbar_k 0]
			// [gbar_k   dbar_{k+1}]   [sn -cs] [alpha_k beta_{k+1}].
			oldeps = epsln;
			delta = cs*dbar + sn*alpha;
			gbar = sn*dbar - cs*alpha;
			epsln = sn*beta;
			dbar = -cs*beta;
			root = sqrt(gbar*gbar + dbar*dbar);
			Arnorm = phibar * root; // ||Ar_{k-1}||

									// Compute next plane rotation Q_k
			gamma = sqrt(gbar*gbar + beta*beta); // gamma_k
			gamma = max(gamma, eps);
			cs = gbar / gamma;                     // c_k
			sn = beta / gamma;                     // s_k
			phi = cs*phibar;                     // phi_k
			phibar = sn*phibar;                  // phibar_{k+1}

												 // Update x
			copy_dvector(&w2, &w1);
			copy_dvector(&w, &w2);
			copy_dvector(&v, &w);
			axpy_dvector(-oldeps, &w1, &w);
			axpy_dvector(-delta, &w2, &w);
			axy_dvector(1. / gamma, &w, &w);
			//			printf("phi=%e, w01=%e, w02=%e, w11=%e, w12=%e\n",phi,w[0].val[1], w[0].val[2], w[1].val[1], w[1].val[2]);
			axpy_dvector(phi, &w, x);

			// go round again
			gmax = max(gmax, gamma);
			gmin = min(gmin, gamma);
			z = rhs1 / gamma;
			rhs1 = rhs2 - delta*z;
			rhs2 = -epsln*z;

			// Estimate various norms
			Anorm = sqrt(tnorm2);
			ynorm2 = dot_dvector(x, x);
			ynorm = sqrt(ynorm2);
			double epsa = Anorm*eps;
			double epsx = epsa*ynorm;
			double epsr = Anorm*ynorm*tol;
			double diag = gbar;
			if (0 == diag)
				diag = epsa;

			qrnorm = phibar;
			rnorm = qrnorm;
			double test1 = 0.0, test2 = 0.0;
			//		test1 = rnorm / (Anorm*ynorm); // ||r||/(||A|| ||x||)
			test1 = rnorm / beta1; // ||r||/(||A|| ||x||)
			test2 = root / Anorm;         // ||A r_{k-1}|| / (||A|| ||r_{k-1}||)

										  //			printf("x=%e, %e, %e, %e, rnorm=%e, beta1=%e\n", x[0].val[1], x[0].val[2], x[1].val[1], x[1].val[2], rnorm, beta1);

										  // Estimate cond(A)
										  /*
										  In this version we look at the diagonals of  R  in the
										  factorization of the lower Hessenberg matrix,  Q * H = R,
										  where H is the tridiagonal matrix from Lanczos with one
										  extra row, beta(k+1) e_k^T.
										  */
			Acond = gmax / gmin;

			//See if any of the stopping criteria is satisfied
			if (0 == istop)
			{
				double t1 = 1.0 + test1, t2 = 1.0 + test2; //This test work if tol < eps
				if (t2 <= 1.) istop = 2;
				if (t1 <= 1.) istop = 1;
				if (iter >= MaxIt - 1) istop = 6;
				if (Acond >= .1 / eps) istop = 4;
				if (epsx >= beta1)   istop = 3;
				if (test2 <= tol)  istop = 2;
				if (test1 <= tol)    istop = 1;
			}

			// output iteration information if needed	
			if (print_level>1) {
				if (iter == 0) {
					printf("It Num |  ||r||/||b|| | ||Ar||/(||A|| ||r||) |     Anorm    |     Acond    |  gbar/Anorm\n");
				}
				printf("%6d | %12.5e | %20.5e | %12.5e | %12.5e | %12.5e\n", iter, test1, test2, Anorm, Acond, gbar / Anorm);
			}


			/*			if (show)
			std::cout << std::setw(6) << itn
			<< std::setw(14) << test1
			<< std::setw(14) << test2
			<< std::setw(14) << Anorm
			<< std::setw(14) << Acond
			<< std::setw(14) << gbar / Anorm << std::endl;*/

			if (0 != istop)
				break;
		}
	}

	// Display final status

	switch (istop)
	{
	case 0:printf(" beta1 = 0.  The exact solution is  x = 0 \n"); break;
	case 1:printf(" A solution to Ax = b was found, given tol \n"); break;
	case 2:printf(" A least-squares solution was found, given tol \n"); break;
	case 3:printf(" Reasonable accuracy achieved, given eps \n"); break;
	case 4:printf(" x has converged to an eigenvector \n"); break;
	case 5:printf(" acond has exceeded 0.1/eps \n"); break;
	case 6:printf(" The iteration limit was reached \n"); break;
	case 7:printf(" A  does not define a symmetric matrix \n"); break;
	case 8:printf(" M  does not define a symmetric matrix \n"); break;
	case 9:printf(" M  does not define a pos-def preconditioner \n"); break;
	default:printf(" beta2 = 0.  If M = I, b and x are eigenvectors \n");
	}
	printf(" Number of iterations: %d\n", iter);
	printf(" Anorm = %e\t Acond = %e\n", Anorm, Acond);
	printf(" rnorm = %e\t ynorm = %e\n", rnorm, ynorm);
	printf(" Arnorm = %e\n", Arnorm);
	/*if (show)
	{
	std::cout << std::setfill('-') << std::setw(80) << "-" << "\n";
	std::cout << msg[istop] << "\n";
	std::cout << " Number of iterations: " << itn << "\n";
	std::cout << " Anorm = " << Anorm << "\t Acond = " << Acond << "\n";
	std::cout << " rnorm = " << rnorm << "\t ynorm = " << ynorm << "\n";
	std::cout << " Arnorm = " << Arnorm << "\n";
	std::cout << std::setfill('-') << std::setw(80) << "-" << std::endl;
	std::cout << std::setfill(' ');
	}*/

	free_dvector(&v);
	free_dvector(&y);
	free_dvector(&r1);
	free_dvector(&r2);

	free_dvector(&w);
	free_dvector(&w1);
	free_dvector(&w2);
	//	return istop;

	return iter;
}

/**
* \fn int minres2b(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, precond *pre, int print_level)
*	 \brief A preconditioned minimal residual (MINRES) method for solving Au=b in 2 blocks
*	 \param *A	 pointer to the coefficient matrix
*	 \param *b	 pointer to the dvector of right hand side
*	 \param *x	 pointer to the dvector of DOFs
*	 \param MaxIt integer, maximal number of iterations
*	 \param tol double float, the tolerance for stopage
*	 \param *pre pointer to the structure of precondition (precond)
*	 \param print_level how much information to print out
*	 \return the number of iterations
*/
int minres2b(dCSRmat *A, dvector *b, dvector *x, int MaxIt, double tol, precond *pre, int print_level)
{
	double eps = DBL_EPSILON;
	int i, iter = 0, n[2];
	int istop = 0;
	double Anorm = 0.0, Acond = 0.0, Arnorm = 0.0;
	double rnorm = 0.0, ynorm = 0.0;
	int done = 0;

	n[0] = b[0].row;
	n[1] = b[1].row;

	// Step 1
	/*
	 * Set up y and v for the first Lanczos vector v1.
	 * y = beta1 P' v1, where P = C^(-1).
	 * v is really P'v1
	 */

	dvector y[2], r1[2];
	for (i = 0; i < 2; i++)
	{
		create_dvector(n[i], &y[i]);
		create_dvector(n[i], &r1[i]);
	}

	// r = b-A*x
	copy_dvector2b(b, r1);
	sparse_mv2b(-1.0, A, x, r1);

	// y = B*r1
	if (pre == NULL)
		copy_dvector2b(r1, y); /* No preconditioner, B=I */
	else
		pre->fct_dvec(r1, y, pre->data); /* Preconditioning */

	double beta1= dot_dvector2b(r1, y);


	// Test for an indefined preconditioner
	// If b = 0 exactly stop with x = x0.
	if (beta1 < 0.0)
	{
		istop = 9;
		done = 1;
		// show = true;
	}
	else
	{
		if (beta1 < eps)
		{
			done = 1;
			// show = true;
		}
		else
			beta1 = sqrt(beta1); // Normalize y to get v1 later
	}

	// STEP 2
	/* Initialize other quantities */
	double oldb = 0.0, beta = beta1, dbar = 0.0, epsln = 0.0, oldeps = 0.0;
	double qrnorm = beta1, phi = 0.0, phibar = beta1, rhs1 = beta1;
	double rhs2 = 0.0, tnorm2 = 0.0, ynorm2 = 0.0;
	double cs = -1.0, sn = 0.0;
	double gmax = 0.0, gmin = DBL_MAX;
	double alpha = 0.0, gamma = 0.0;
	double delta = 0.0, gbar = 0.0;
	double z = 0.0;
	double root;
	
	dvector w[2], w1[2], w2[2], r2[2], v[2];
	for (i = 0; i < 2; i++)
	{
		create_dvector(n[i], &w[i]);
		create_dvector(n[i], &w1[i]);
		create_dvector(n[i], &w2[i]);
		create_dvector(n[i], &r2[i]);
		create_dvector(n[i], &v[i]);
	}

	copy_dvector2b(r1, r2);

	/* Main Iteration */
	if (!done)
	{
		for (iter = 0; iter < MaxIt; ++iter)
		{
			// STEP 3
			/*
			-----------------------------------------------------------------
			Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
			The general iteration is similar to the case k = 1 with v0 = 0:
			  p1      = A * v1  -  beta1 * v0,
			  alpha1  = v1'p1,
			  q2      = p2  -  alpha1 * v1,
			  beta2^2 = q2'q2,
			  v2      = (1/beta2) q2.

			Again, y = betak P vk,  where  P = C**(-1).
			.... more description needed.
			-----------------------------------------------------------------
				*/
			axy_dvector2b(1. / beta, y, v); // Normalize previous vector (in y)

											// y = A*v-r1*beta/oldb
			init_dvector2b(y, 0);
			sparse_mv2b(1.0, A, v, y);
			if(iter)
				axpy_dvector2b(-beta / oldb, r1, y);

			alpha = dot_dvector2b(v, y);   // alphak
										   // y = y-r2*alpha/beta
			axpy_dvector2b(-alpha / beta, r2, y);
			copy_dvector2b(r2, r1);
			copy_dvector2b(y, r2);

			// y = B*r2
			if (pre == NULL)
				copy_dvector2b(r2, y); /* No preconditioner, B=I */
			else
				pre->fct_dvec(r2, y, pre->data); /* Preconditioning */

			oldb = beta; //oldb = betak
			beta = dot_dvector2b(r2, y);
			if (beta < 0)
			{
				istop = 9;
				break;
			}
			beta = sqrt(beta);
			
			tnorm2 += alpha*alpha + oldb*oldb + beta*beta;

			if (iter == 0)    //Initialize a few things
			{
				if (beta / beta1 < 10.0*eps)
					istop = 10;
			}

			// Apply previous rotation Q_{k-1} to get
			// [delta_k epsln_{k+1}] = [cs sn]  [dbar_k 0]
			// [gbar_k   dbar_{k+1}]   [sn -cs] [alpha_k beta_{k+1}].
			oldeps = epsln;
			delta = cs*dbar + sn*alpha;
			gbar = sn*dbar - cs*alpha;
			epsln = sn*beta;
			dbar = -cs*beta;
			root = sqrt(gbar*gbar + dbar*dbar);
			Arnorm = phibar * root; // ||Ar_{k-1}||

									// Compute next plane rotation Q_k
			gamma = sqrt(gbar*gbar + beta*beta); // gamma_k
			gamma = max(gamma, eps);
			cs = gbar / gamma;                     // c_k
			sn = beta / gamma;                     // s_k
			phi = cs*phibar;                     // phi_k
			phibar = sn*phibar;                  // phibar_{k+1}

												 // Update x
			copy_dvector2b(w2, w1);
			copy_dvector2b(w, w2);
			copy_dvector2b(v, w);
			axpy_dvector2b(-oldeps, w1, w);
			axpy_dvector2b(-delta, w2, w);
			axy_dvector2b(1. / gamma, w, w);
//			printf("phi=%e, w01=%e, w02=%e, w11=%e, w12=%e\n",phi,w[0].val[1], w[0].val[2], w[1].val[1], w[1].val[2]);
			axpy_dvector2b(phi, w, x);
			
			// go round again
			gmax = max(gmax, gamma);
			gmin = min(gmin, gamma);
			z = rhs1 / gamma;
			rhs1 = rhs2 - delta*z;
			rhs2 = -epsln*z;

			// Estimate various norms
			Anorm = sqrt(tnorm2);
			ynorm2 = dot_dvector2b(x, x);
			ynorm = sqrt(ynorm2);
			double epsa = Anorm*eps;
			double epsx = epsa*ynorm;
			double epsr = Anorm*ynorm*tol;
			double diag = gbar;
			if (0 == diag)
				diag = epsa;

			qrnorm = phibar;
			rnorm = qrnorm;
			double test1 = 0.0, test2 = 0.0;
			//		test1 = rnorm / (Anorm*ynorm); // ||r||/(||A|| ||x||)
			test1 = rnorm / beta1; // ||r||/(||A|| ||x||)
			test2 = root / Anorm;         // ||A r_{k-1}|| / (||A|| ||r_{k-1}||)

//			printf("x=%e, %e, %e, %e, rnorm=%e, beta1=%e\n", x[0].val[1], x[0].val[2], x[1].val[1], x[1].val[2], rnorm, beta1);

			// Estimate cond(A)
			/*
			In this version we look at the diagonals of  R  in the
			factorization of the lower Hessenberg matrix,  Q * H = R,
			where H is the tridiagonal matrix from Lanczos with one
			extra row, beta(k+1) e_k^T.
			*/
			Acond = gmax / gmin;

			//See if any of the stopping criteria is satisfied
			if (0 == istop)
			{
				double t1 = 1.0 + test1, t2 = 1.0 + test2; //This test work if tol < eps
				if (t2 <= 1.) istop = 2;
				if (t1 <= 1.) istop = 1;
				if (iter >= MaxIt - 1) istop = 6;
				if (Acond >= .1 / eps) istop = 4;
				if (epsx >= beta1)   istop = 3;
				if (test2 <= tol)  istop = 2;
				if (test1 <= tol)    istop = 1;
			}

			// output iteration information if needed	
			if (print_level>1) {
				if (iter == 0) {
					printf("It Num |  ||r||/||b|| | ||Ar||/(||A|| ||r||) |     Anorm    |     Acond    |  gbar/Anorm\n");
				}
				printf("%6d | %12.5e | %20.5e | %12.5e | %12.5e | %12.5e\n", iter, test1, test2, Anorm, Acond, gbar/Anorm);
			}


/*			if (show)
				std::cout << std::setw(6) << itn
				<< std::setw(14) << test1
				<< std::setw(14) << test2
				<< std::setw(14) << Anorm
				<< std::setw(14) << Acond
				<< std::setw(14) << gbar / Anorm << std::endl;*/

			if (0 != istop)
				break;
		}
	}

	// Display final status

	switch (istop)
	{
	case 0:printf(" beta1 = 0.  The exact solution is  x = 0 \n"); break;
	case 1:printf(" A solution to Ax = b was found, given tol \n"); break;
	case 2:printf(" A least-squares solution was found, given tol \n"); break;
	case 3:printf(" Reasonable accuracy achieved, given eps \n"); break;
	case 4:printf(" x has converged to an eigenvector \n"); break;
	case 5:printf(" acond has exceeded 0.1/eps \n"); break;
	case 6:printf(" The iteration limit was reached \n"); break;
	case 7:printf(" A  does not define a symmetric matrix \n"); break;
	case 8:printf(" M  does not define a symmetric matrix \n"); break;
	case 9:printf(" M  does not define a pos-def preconditioner \n"); break;
	default:printf(" beta2 = 0.  If M = I, b and x are eigenvectors \n");
	}
	printf(" Number of iterations: %d\n", iter);
	printf(" Anorm = %e\t Acond = %e\n", Anorm, Acond);
	printf(" rnorm = %e\t ynorm = %e\n", rnorm, ynorm);
	printf(" Arnorm = %e\n", Arnorm);
	/*if (show)
	{
		std::cout << std::setfill('-') << std::setw(80) << "-" << "\n";
		std::cout << msg[istop] << "\n";
		std::cout << " Number of iterations: " << itn << "\n";
		std::cout << " Anorm = " << Anorm << "\t Acond = " << Acond << "\n";
		std::cout << " rnorm = " << rnorm << "\t ynorm = " << ynorm << "\n";
		std::cout << " Arnorm = " << Arnorm << "\n";
		std::cout << std::setfill('-') << std::setw(80) << "-" << std::endl;
		std::cout << std::setfill(' ');
	}*/

	for (i = 0; i < 2; i++)
	{
		free_dvector(&v[i]);
		free_dvector(&y[i]);
		free_dvector(&r1[i]);
		free_dvector(&r2[i]);

		free_dvector(&w[i]);
		free_dvector(&w1[i]);
		free_dvector(&w2[i]);
	}
//	return istop;
	
	return iter;
}
