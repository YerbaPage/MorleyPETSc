/*
 *  basis.c
 *
 *  Created by Xuehai Huang on 3/29/09.
 *  Copyright 2009 PSU. All rights reserved.
 *
 */

/*! \file basis.c
 *  \brief Basis functions
 */
 
#include <math.h>
#include "header.h"

#define omega PI*3.0/4.0
#define z 0.544483736782464 // omega=PI*3.0/4.0

 /**
 * \fn void cart2pol(x,y,r,theta)
 * \brief Transform Cartesian coordinates to polar
 * \param x the x-axis value of the point
 * \param y the y-axis value of the point
 * \param r the distance from the origin to a point in the x-y plane
 * \param theta a counterclockwise angular displacement in radians from the positive x-axis
 * \return void
 */
void cart2pol(double x, double y, double *r, double *theta)
{
	(*r) = sqrt(x*x + y*y);
	(*theta) = atan2(y, x);
}

/**
* \fn double f(double x, double y)
* \brief load f, i.e. right hand side when the true solution u is x*x*(x-1)*(x-1)*y*y*(y-1)*(y-1)
*		  \Delta^2 u = f
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \return function value
*/
double f(double x, double y)
{
	//	return 2*PI*PI*sin(PI*x)*sin(PI*y);

	double a, b, c;
	a = x*x - x;
	b = y*y - y;
	c = a*a + b*b + 12 * a*b + 2 * a + 2 * b;
	return c * 24 + 8;
}

/**
* \fn double u(double x, double y)
* \brief true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \return function value
*/
double u(double x, double y)
{
	return x*x*(x - 1)*(x - 1)*y*y*(y - 1)*(y - 1);
}

/**
* \fn double u_x(double x, double y)
* \brief the x-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \return x-directional partial derivative of true solution u
*/
double u_x(double x, double y)
{
	return (4 * x*x*x - 6 * x*x + 2 * x)*y*y*(y - 1)*(y - 1);
}

/**
* \fn double u_y(double x, double y)
* \brief the y-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \return y-directional partial derivative of true solution u
*/
double u_y(double x, double y)
{
	return x*x*(x - 1)*(x - 1)*(4 * y*y*y - 6 * y*y + 2 * y);
}

/**
* \fn double u_xx(double x, double y)
* \brief the xx-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \return xx-directional partial derivative of true solution u
*/
double u_xx(double x, double y)
{
	return (12 * x*x - 12 * x + 2)*y*y*(y - 1)*(y - 1);
}

/**
* \fn double u_xy(double x, double y)
* \brief the xy-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \return xy-directional partial derivative of true solution u
*/
double u_xy(double x, double y)
{
	return (4 * x*x*x - 6 * x*x + 2 * x)*(4 * y*y*y - 6 * y*y + 2 * y);
}

/**
* \fn double u_yy(double x, double y)
* \brief the yy-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \return yy-directional partial derivative of true solution u
*/
double u_yy(double x, double y)
{
	return x*x*(x - 1)*(x - 1)*(12 * y*y - 12 * y + 2);
}

double f1(double x, double y, double lambda, double mu)
{
	double a = (3 * x*y - 2)*(x*x + y*y) + 5 * pow(x*y - 1, 2) - 2 * x*x*y*y;
	return 1;
	return -8 * (x + y)*a;

}

double f2(double x, double y, double lambda, double mu)
{
	double a = (3 * x*y + 2)*(x*x + y*y) - 5 * pow(x*y + 1, 2) + 2 * x*x*y*y;
	return 1;
	return -8 * (x - y)*a;
}

/**
* \fn double u1(double x, double y, double lambda, double mu)
* \brief true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param nu Lame constant or Poisson ratio of plate
* \return function value
*/
double u1(double x, double y, double lambda, double mu)
{
//	return x*x + y*y;
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t0 = (pow(r, z)*(t2*cos(t*t2) - (cos(omega*t2)*cos(t*t3)*(z - (lambda*2.0 + mu*4.0) / (lambda + mu) + 1.0)) / cos(omega*t3))*(-1.0 / 2.0)) / mu;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u2(double x, double y, double lambda, double mu)
* \brief true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return function value
*/
double u2(double x, double y, double lambda, double mu)
{
//	return x*x*10 + y*y*100;
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t0 = (pow(r, z)*(t2*sin(t*t2) - (cos(omega*t2)*sin(t*t3)*(z + (lambda*2.0 + mu*4.0) / (lambda + mu) - 1.0)) / cos(omega*t3))*(1.0 / 2.0)) / mu;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u1_x(double x, double y, double lambda, double mu)
* \brief the x-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return x-directional partial derivative of true solution u
*/
double u1_x(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = t*t2;
	t6 = omega*t3;
	t7 = cos(t6);
	t8 = 1.0 / t7;
	t9 = omega*t2;
	t10 = cos(t9);
	t11 = t*t3;
	t12 = lambda + mu;
	t13 = 1.0 / t12;
	t14 = lambda*2.0;
	t15 = mu*4.0;
	t16 = t14 + t15;
	t17 = z - t13*t16 + 1.0;
	t0 = pow(r, t3)*t4*z*cos(t)*(t2*cos(t5) - t8*t10*t17*cos(t11))*(-1.0 / 2.0) - (pow(r, z)*t4*sin(t)*((t2*t2)*sin(t5) - t3*t8*t10*t17*sin(t11))*(1.0 / 2.0)) / r;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u1_y(double x, double y, double lambda, double mu)
* \brief the y-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return y-directional partial derivative of true solution u
*/
double u1_y(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = t*t2;
	t6 = omega*t3;
	t7 = cos(t6);
	t8 = 1.0 / t7;
	t9 = omega*t2;
	t10 = cos(t9);
	t11 = t*t3;
	t12 = lambda + mu;
	t13 = 1.0 / t12;
	t14 = lambda*2.0;
	t15 = mu*4.0;
	t16 = t14 + t15;
	t17 = z - t13*t16 + 1.0;
	t0 = pow(r, t3)*t4*z*sin(t)*(t2*cos(t5) - t8*t10*t17*cos(t11))*(-1.0 / 2.0) + (pow(r, z)*t4*cos(t)*((t2*t2)*sin(t5) - t3*t8*t10*t17*sin(t11))*(1.0 / 2.0)) / r;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u2_x(double x, double y, double lambda, double mu)
* \brief the x-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return x-directional partial derivative of true solution u
*/
double u2_x(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = t*t2;
	t6 = omega*t3;
	t7 = cos(t6);
	t8 = 1.0 / t7;
	t9 = omega*t2;
	t10 = cos(t9);
	t11 = t*t3;
	t12 = lambda + mu;
	t13 = 1.0 / t12;
	t14 = lambda*2.0;
	t15 = mu*4.0;
	t16 = t14 + t15;
	t17 = t13*t16;
	t18 = t17 + z - 1.0;
	t0 = pow(r, t3)*t4*z*cos(t)*(t2*sin(t5) - t8*t10*t18*sin(t11))*(1.0 / 2.0) - (pow(r, z)*t4*sin(t)*((t2*t2)*cos(t5) - t3*t8*t10*t18*cos(t11))*(1.0 / 2.0)) / r;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u2_y(double x, double y, double lambda, double mu)
* \brief the y-directional partial derivative of true solution u
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return y-directional partial derivative of true solution u
*/
double u2_y(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = t*t2;
	t6 = omega*t3;
	t7 = cos(t6);
	t8 = 1.0 / t7;
	t9 = omega*t2;
	t10 = cos(t9);
	t11 = t*t3;
	t12 = lambda + mu;
	t13 = 1.0 / t12;
	t14 = lambda*2.0;
	t15 = mu*4.0;
	t16 = t14 + t15;
	t17 = t13*t16;
	t18 = t17 + z - 1.0;
	t0 = pow(r, t3)*t4*z*sin(t)*(t2*sin(t5) - t8*t10*t18*sin(t11))*(1.0 / 2.0) + (pow(r, z)*t4*cos(t)*((t2*t2)*cos(t5) - t3*t8*t10*t18*cos(t11))*(1.0 / 2.0)) / r;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u1_xx(double x, double y, double lambda, double mu)
* \brief the xx-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return xx-directional partial derivative of true solution u1
*/
double u1_xx(double x, double y, double lambda, double mu)
{
	//	return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = cos(t);
	t6 = t*t2;
	t7 = omega*t3;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t2;
	t11 = cos(t10);
	t12 = t*t3;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t24 = t14*t17;
	t18 = -t24 + z + 1.0;
	t19 = sin(t);
	t20 = sin(t6);
	t21 = t2*t2;
	t22 = t20*t21;
	t23 = sin(t12);
	t33 = t3*t9*t11*t18*t23;
	t25 = t22 - t33;
	t26 = 1.0 / r;
	t27 = pow(r, z);
	t28 = cos(t6);
	t29 = cos(t12);
	t30 = pow(r, t3);
	t31 = t2*t28;
	t32 = t31 - t9*t11*t18*t29;
	t0 = -t5*(1.0 / (r*r)*t4*t19*t25*t27*(-1.0 / 2.0) + t4*t19*t25*t26*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t3*t4*t5*t32*z*(1.0 / 2.0)) + t19*t26*(t4*t19*t26*t27*(t2*t21*t28 - (t3*t3)*t9*t11*t18*t29)*(1.0 / 2.0) + t4*t5*t25*t26*t27*(1.0 / 2.0) - t4*t5*t25*t30*z*(1.0 / 2.0) - t4*t19*t30*t32*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u1_xy(double x, double y, double lambda, double mu)
* \brief the xy-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return xy-directional partial derivative of true solution u1
*/
double u1_xy(double x, double y, double lambda, double mu)
{
	//		return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = sin(t);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t23 = t14*t17;
	t18 = -t23 + z + 1.0;
	t19 = sin(t6);
	t20 = t3*t3;
	t21 = t19*t20;
	t22 = sin(t12);
	t33 = t4*t9*t11*t18*t22;
	t24 = t21 - t33;
	t25 = 1.0 / r;
	t26 = cos(t);
	t27 = pow(r, z);
	t28 = cos(t6);
	t29 = cos(t12);
	t30 = pow(r, t4);
	t31 = t3*t28;
	t32 = t31 - t9*t11*t18*t29;
	t0 = -t2*(1.0 / (r*r)*t2*t5*t24*t27*(-1.0 / 2.0) + t2*t5*t24*t25*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t4*t5*t26*t32*z*(1.0 / 2.0)) - t25*t26*(t2*t5*t25*t27*(t3*t20*t28 - (t4*t4)*t9*t11*t18*t29)*(1.0 / 2.0) + t5*t24*t25*t26*t27*(1.0 / 2.0) - t2*t5*t30*t32*z*(1.0 / 2.0) - t5*t24*t26*t30*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u1_yx(double x, double y, double lambda, double mu)
* \brief the yx-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return yx-directional partial derivative of true solution u1
*/
double u1_yx(double x, double y, double lambda, double mu)
{
	//		return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = cos(t);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t23 = t14*t17;
	t18 = -t23 + z + 1.0;
	t19 = sin(t6);
	t20 = t3*t3;
	t21 = t19*t20;
	t22 = sin(t12);
	t33 = t4*t9*t11*t18*t22;
	t24 = t21 - t33;
	t25 = 1.0 / r;
	t26 = sin(t);
	t27 = pow(r, z);
	t28 = cos(t6);
	t29 = cos(t12);
	t30 = pow(r, t4);
	t31 = t3*t28;
	t32 = t31 - t9*t11*t18*t29;
	t0 = -t2*(1.0 / (r*r)*t2*t5*t24*t27*(1.0 / 2.0) - t2*t5*t24*t25*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t4*t5*t26*t32*z*(1.0 / 2.0)) - t25*t26*(t2*t5*t25*t27*(t3*t20*t28 - (t4*t4)*t9*t11*t18*t29)*(1.0 / 2.0) - t5*t24*t25*t26*t27*(1.0 / 2.0) - t2*t5*t30*t32*z*(1.0 / 2.0) + t5*t24*t26*t30*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u1_yy(double x, double y, double lambda, double mu)
* \brief the yy-directional partial derivative of true solution u1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return yy-directional partial derivative of true solution u1
*/
double u1_yy(double x, double y, double lambda, double mu)
{
	//		return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = sin(t);
	t6 = t*t2;
	t7 = omega*t3;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t2;
	t11 = cos(t10);
	t12 = t*t3;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t24 = t14*t17;
	t18 = -t24 + z + 1.0;
	t19 = cos(t);
	t20 = sin(t6);
	t21 = t2*t2;
	t22 = t20*t21;
	t23 = sin(t12);
	t33 = t3*t9*t11*t18*t23;
	t25 = t22 - t33;
	t26 = 1.0 / r;
	t27 = pow(r, z);
	t28 = cos(t6);
	t29 = cos(t12);
	t30 = pow(r, t3);
	t31 = t2*t28;
	t32 = t31 - t9*t11*t18*t29;
	t0 = -t5*(1.0 / (r*r)*t4*t19*t25*t27*(1.0 / 2.0) - t4*t19*t25*t26*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t3*t4*t5*t32*z*(1.0 / 2.0)) + t19*t26*(t4*t19*t26*t27*(t2*t21*t28 - (t3*t3)*t9*t11*t18*t29)*(1.0 / 2.0) - t4*t5*t25*t26*t27*(1.0 / 2.0) + t4*t5*t25*t30*z*(1.0 / 2.0) - t4*t19*t30*t32*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u2_xx(double x, double y, double lambda, double mu)
* \brief the xx-directional partial derivative of true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return xx-directional partial derivative of true solution u2
*/
double u2_xx(double x, double y, double lambda, double mu)
{
	//		return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = cos(t);
	t6 = t*t2;
	t7 = omega*t3;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t2;
	t11 = cos(t10);
	t12 = t*t3;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = sin(t);
	t21 = cos(t6);
	t22 = t2*t2;
	t23 = t21*t22;
	t24 = cos(t12);
	t33 = t3*t9*t11*t19*t24;
	t25 = t23 - t33;
	t26 = 1.0 / r;
	t27 = pow(r, z);
	t28 = sin(t6);
	t29 = sin(t12);
	t30 = pow(r, t3);
	t31 = t2*t28;
	t32 = t31 - t9*t11*t19*t29;
	t0 = t5*(1.0 / (r*r)*t4*t20*t25*t27*(1.0 / 2.0) - t4*t20*t25*t26*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t3*t4*t5*t32*z*(1.0 / 2.0)) - t20*t26*(t4*t20*t26*t27*(t2*t22*t28 - (t3*t3)*t9*t11*t19*t29)*(1.0 / 2.0) - t4*t5*t25*t26*t27*(1.0 / 2.0) + t4*t5*t25*t30*z*(1.0 / 2.0) - t4*t20*t30*t32*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u2_xy(double x, double y, double lambda, double mu)
* \brief the xy-directional partial derivative of true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return xy-directional partial derivative of true solution u2
*/
double u2_xy(double x, double y, double lambda, double mu)
{
	//		return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = sin(t);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = cos(t6);
	t21 = t3*t3;
	t22 = t20*t21;
	t23 = cos(t12);
	t33 = t4*t9*t11*t19*t23;
	t24 = t22 - t33;
	t25 = 1.0 / r;
	t26 = cos(t);
	t27 = pow(r, z);
	t28 = sin(t6);
	t29 = sin(t12);
	t30 = pow(r, t4);
	t31 = t3*t28;
	t32 = t31 - t9*t11*t19*t29;
	t0 = t2*(1.0 / (r*r)*t2*t5*t24*t27*(1.0 / 2.0) - t2*t5*t24*t25*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t4*t5*t26*t32*z*(1.0 / 2.0)) + t25*t26*(t2*t5*t25*t27*(t3*t21*t28 - (t4*t4)*t9*t11*t19*t29)*(1.0 / 2.0) - t5*t24*t25*t26*t27*(1.0 / 2.0) - t2*t5*t30*t32*z*(1.0 / 2.0) + t5*t24*t26*t30*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u2_yx(double x, double y, double lambda, double mu)
* \brief the yx-directional partial derivative of true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return yx-directional partial derivative of true solution u2
*/
double u2_yx(double x, double y, double lambda, double mu)
{
	//		return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = cos(t);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = cos(t6);
	t21 = t3*t3;
	t22 = t20*t21;
	t23 = cos(t12);
	t33 = t4*t9*t11*t19*t23;
	t24 = t22 - t33;
	t25 = 1.0 / r;
	t26 = sin(t);
	t27 = pow(r, z);
	t28 = sin(t6);
	t29 = sin(t12);
	t30 = pow(r, t4);
	t31 = t3*t28;
	t32 = t31 - t9*t11*t19*t29;
	t0 = t2*(1.0 / (r*r)*t2*t5*t24*t27*(-1.0 / 2.0) + t2*t5*t24*t25*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t4*t5*t26*t32*z*(1.0 / 2.0)) + t25*t26*(t2*t5*t25*t27*(t3*t21*t28 - (t4*t4)*t9*t11*t19*t29)*(1.0 / 2.0) + t5*t24*t25*t26*t27*(1.0 / 2.0) - t2*t5*t30*t32*z*(1.0 / 2.0) - t5*t24*t26*t30*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double u2_yy(double x, double y, double lambda, double mu)
* \brief the yy-directional partial derivative of true solution u2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return yy-directional partial derivative of true solution u2
*/
double u2_yy(double x, double y, double lambda, double mu)
{
	//		return 0;
	/**for rotated L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = sin(t);
	t6 = t*t2;
	t7 = omega*t3;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t2;
	t11 = cos(t10);
	t12 = t*t3;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = cos(t);
	t21 = cos(t6);
	t22 = t2*t2;
	t23 = t21*t22;
	t24 = cos(t12);
	t33 = t3*t9*t11*t19*t24;
	t25 = t23 - t33;
	t26 = 1.0 / r;
	t27 = pow(r, z);
	t28 = sin(t6);
	t29 = sin(t12);
	t30 = pow(r, t3);
	t31 = t2*t28;
	t32 = t31 - t9*t11*t19*t29;
	t0 = t5*(1.0 / (r*r)*t4*t20*t25*t27*(-1.0 / 2.0) + t4*t20*t25*t26*t30*z*(1.0 / 2.0) + pow(r, z - 2.0)*t3*t4*t5*t32*z*(1.0 / 2.0)) - t20*t26*(t4*t20*t26*t27*(t2*t22*t28 - (t3*t3)*t9*t11*t19*t29)*(1.0 / 2.0) + t4*t5*t25*t26*t27*(1.0 / 2.0) - t4*t5*t25*t30*z*(1.0 / 2.0) - t4*t20*t30*t32*z*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double sigma11(double x, double y, double lambda, double mu)
* \brief sigma11
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma11
*/
double sigma11(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = t*t2;
	t6 = omega*t3;
	t7 = cos(t6);
	t8 = 1.0 / t7;
	t9 = omega*t2;
	t10 = cos(t9);
	t11 = t*t3;
	t12 = lambda + mu;
	t13 = 1.0 / t12;
	t14 = lambda*2.0;
	t15 = mu*4.0;
	t16 = t14 + t15;
	t23 = t13*t16;
	t17 = -t23 + z + 1.0;
	t18 = pow(r, t3);
	t19 = cos(t);
	t20 = cos(t5);
	t21 = t2*t20;
	t22 = cos(t11);
	t24 = t21 - t8*t10*t17*t22;
	t25 = t4*t18*t19*t24*z*(1.0 / 2.0);
	t26 = sin(t);
	t27 = sin(t5);
	t28 = sin(t11);
	t29 = pow(r, z);
	t30 = 1.0 / r;
	t31 = t2*t2;
	t32 = t23 + z - 1.0;
	t33 = t27*t31;
	t34 = t33 - t3*t8*t10*t17*t28;
	t35 = t4*t26*t29*t30*t34*(1.0 / 2.0);
	t0 = mu*(t25 + t35)*-2.0 - lambda*(t25 + t35 - t4*t18*t26*z*(t2*t27 - t8*t10*t28*t32)*(1.0 / 2.0) - t4*t19*t29*t30*(t20*t31 - t3*t8*t10*t22*t32)*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double sigma22(double x, double y, double lambda, double mu)
* \brief sigma22
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma22
*/
double sigma22(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
	double t61, t62, t63, t64, t65, t66, t67, t68, t69;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = t*t2;
	t6 = omega*t3;
	t7 = cos(t6);
	t8 = 1.0 / t7;
	t9 = omega*t2;
	t10 = cos(t9);
	t11 = t*t3;
	t12 = lambda + mu;
	t13 = 1.0 / t12;
	t14 = lambda*2.0;
	t15 = mu*4.0;
	t16 = t14 + t15;
	t17 = t13*t16;
	t18 = t17 + z - 1.0;
	t19 = pow(r, t3);
	t20 = cos(t);
	t21 = cos(t5);
	t22 = cos(t11);
	t23 = sin(t);
	t24 = sin(t5);
	t25 = t2*t24;
	t26 = sin(t11);
	t27 = t25 - t8*t10*t18*t26;
	t28 = t4*t19*t23*t27*z*(1.0 / 2.0);
	t29 = pow(r, z);
	t30 = 1.0 / r;
	t31 = t2*t2;
	t32 = t21*t31;
	t33 = t32 - t3*t8*t10*t18*t22;
	t34 = t4*t20*t29*t30*t33*(1.0 / 2.0);
	t35 = -t17 + z + 1.0;
	t0 = mu*(t28 + t34)*2.0 + lambda*(t28 + t34 - t4*t19*t20*z*(t2*t21 - t8*t10*t22*t35)*(1.0 / 2.0) - t4*t23*t29*t30*(t24*t31 - t3*t8*t10*t26*t35)*(1.0 / 2.0));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double sigma12(double x, double y, double lambda, double mu)
* \brief sigma12
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma12
*/
double sigma12(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = pow(r, t3);
	t6 = t*t2;
	t7 = omega*t3;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t2;
	t11 = cos(t10);
	t12 = t*t3;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = sin(t);
	t19 = cos(t6);
	t20 = cos(t12);
	t21 = t14*t17;
	t22 = t21 + z - 1.0;
	t23 = pow(r, z);
	t24 = 1.0 / r;
	t25 = cos(t);
	t26 = sin(t6);
	t27 = t2*t2;
	t28 = sin(t12);
	t0 = mu*(t4*t5*t25*z*(t2*t26 - t9*t11*t22*t28)*(1.0 / 4.0) - t4*t18*t23*t24*(t19*t27 - t3*t9*t11*t20*t22)*(1.0 / 4.0) + t4*t23*t24*t25*(t26*t27 - t3*t9*t11*t28*(-t21 + z + 1.0))*(1.0 / 4.0) - t4*t5*t18*z*(t2*t19 - t9*t11*t20*(z - t14*t17 + 1.0))*(1.0 / 4.0))*2.0;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double sigma21(double x, double y, double lambda, double mu)
* \brief sigma21
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return sigma21
*/
double sigma21(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44, t45, t46, t47, t48, t49, t50;
	double t51, t52, t53, t54, t55;
	t2 = z + 1.0;
	t3 = z - 1.0;
	t4 = 1.0 / mu;
	t5 = pow(r, t3);
	t6 = t*t2;
	t7 = omega*t3;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t2;
	t11 = cos(t10);
	t12 = t*t3;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = sin(t);
	t19 = cos(t6);
	t20 = cos(t12);
	t21 = t14*t17;
	t22 = t21 + z - 1.0;
	t23 = pow(r, z);
	t24 = 1.0 / r;
	t25 = cos(t);
	t26 = sin(t6);
	t27 = t2*t2;
	t28 = sin(t12);
	t0 = mu*(t4*t5*t25*z*(t2*t26 - t9*t11*t22*t28)*(1.0 / 4.0) - t4*t18*t23*t24*(t19*t27 - t3*t9*t11*t20*t22)*(1.0 / 4.0) + t4*t23*t24*t25*(t26*t27 - t3*t9*t11*t28*(-t21 + z + 1.0))*(1.0 / 4.0) - t4*t5*t18*z*(t2*t19 - t9*t11*t20*(z - t14*t17 + 1.0))*(1.0 / 4.0))*2.0;

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double ut_t(double x, double y, double lambda, double mu, int edgeFlag)
* \brief the tangential derivative of u \dot t on boundary
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential derivative of u \dot t on boundary
*/
double ut_t(double x, double y, double lambda, double mu, int edgeFlag)
{
	/**for L-shaped domain **/
	double t0;
	switch (edgeFlag)
	{
	case 1: t0 = ut_t_1(x, y, lambda, mu); break;
	case 2: t0 = ut_t_2(x, y, lambda, mu); break;
	case 3: t0 = ut_t_3(x, y, lambda, mu); break;
	case 4: t0 = ut_t_4(x, y, lambda, mu); break;
	default: t0 = 0;
	}
	return t0;
}

/**
* \fn double ut_t_1(double x, double y, double lambda, double mu)
* \brief the tangential derivative of u \dot t on boundary type 1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential derivative of u \dot t on boundary type 1
*/
double ut_t_1(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = sqrt(2.0);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = pow(r, t4);
	t21 = sin(t);
	t22 = sin(t6);
	t23 = t3*t22;
	t24 = sin(t12);
	t25 = t23 - t9*t11*t19*t24;
	t26 = pow(r, z);
	t27 = 1.0 / r;
	t28 = cos(t);
	t29 = cos(t6);
	t30 = t3*t3;
	t31 = t29*t30;
	t32 = cos(t12);
	t33 = t31 - t4*t9*t11*t19*t32;
	t34 = -t18 + z + 1.0;
	t35 = t3*t29;
	t36 = t35 - t9*t11*t32*t34;
	t37 = t22*t30;
	t38 = t37 - t4*t9*t11*t24*t34;
	t0 = t2*(t2*(t5*t21*t26*t27*t33*(1.0 / 2.0) - t5*t20*t25*t28*z*(1.0 / 2.0))*(1.0 / 2.0) + t2*(t5*t26*t27*t28*t33*(1.0 / 2.0) + t5*t20*t21*t25*z*(1.0 / 2.0))*(1.0 / 2.0))*(1.0 / 2.0) - t2*(t2*(t5*t21*t26*t27*t38*(1.0 / 2.0) + t5*t20*t28*t36*z*(1.0 / 2.0))*(1.0 / 2.0) + t2*(t5*t26*t27*t28*t38*(1.0 / 2.0) - t5*t20*t21*t36*z*(1.0 / 2.0))*(1.0 / 2.0))*(1.0 / 2.0);

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double ut_t_2(double x, double y, double lambda, double mu)
* \brief the tangential derivative of u \dot t on boundary type 2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential derivative of u \dot t on boundary type 2
*/
double ut_t_2(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = sqrt(2.0);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = pow(r, t4);
	t21 = sin(t);
	t22 = sin(t6);
	t23 = t3*t22;
	t24 = sin(t12);
	t25 = t23 - t9*t11*t19*t24;
	t26 = pow(r, z);
	t27 = 1.0 / r;
	t28 = cos(t);
	t29 = cos(t6);
	t30 = t3*t3;
	t31 = t29*t30;
	t32 = cos(t12);
	t33 = t31 - t4*t9*t11*t19*t32;
	t34 = -t18 + z + 1.0;
	t35 = t3*t29;
	t36 = t35 - t9*t11*t32*t34;
	t37 = t22*t30;
	t38 = t37 - t4*t9*t11*t24*t34;
	t0 = t2*(t2*(t5*t21*t26*t27*t33*(1.0 / 2.0) - t5*t20*t25*t28*z*(1.0 / 2.0))*(1.0 / 2.0) - t2*(t5*t26*t27*t28*t33*(1.0 / 2.0) + t5*t20*t21*t25*z*(1.0 / 2.0))*(1.0 / 2.0))*(-1.0 / 2.0) - t2*(t2*(t5*t21*t26*t27*t38*(1.0 / 2.0) + t5*t20*t28*t36*z*(1.0 / 2.0))*(1.0 / 2.0) - t2*(t5*t26*t27*t28*t38*(1.0 / 2.0) - t5*t20*t21*t36*z*(1.0 / 2.0))*(1.0 / 2.0))*(1.0 / 2.0);

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double ut_t_3(double x, double y, double lambda, double mu)
* \brief the tangential derivative of u \dot t on boundary type 3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential derivative of u \dot t on boundary type 3
*/
double ut_t_3(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = sqrt(2.0);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = pow(r, t4);
	t21 = sin(t);
	t22 = sin(t6);
	t23 = t3*t22;
	t24 = sin(t12);
	t25 = t23 - t9*t11*t19*t24;
	t26 = pow(r, z);
	t27 = 1.0 / r;
	t28 = cos(t);
	t29 = cos(t6);
	t30 = t3*t3;
	t31 = t29*t30;
	t32 = cos(t12);
	t33 = t31 - t4*t9*t11*t19*t32;
	t34 = -t18 + z + 1.0;
	t35 = t3*t29;
	t36 = t35 - t9*t11*t32*t34;
	t37 = t22*t30;
	t38 = t37 - t4*t9*t11*t24*t34;
	t0 = t2*(t2*(t5*t21*t26*t27*t33*(1.0 / 2.0) - t5*t20*t25*t28*z*(1.0 / 2.0))*(1.0 / 2.0) + t2*(t5*t26*t27*t28*t33*(1.0 / 2.0) + t5*t20*t21*t25*z*(1.0 / 2.0))*(1.0 / 2.0))*(1.0 / 2.0) - t2*(t2*(t5*t21*t26*t27*t38*(1.0 / 2.0) + t5*t20*t28*t36*z*(1.0 / 2.0))*(1.0 / 2.0) + t2*(t5*t26*t27*t28*t38*(1.0 / 2.0) - t5*t20*t21*t36*z*(1.0 / 2.0))*(1.0 / 2.0))*(1.0 / 2.0);

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double ut_t_4(double x, double y, double lambda, double mu)
* \brief the tangential derivative of u \dot t on boundary type 4
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential derivative of u \dot t on boundary type 4
*/
double ut_t_4(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = sqrt(2.0);
	t3 = z + 1.0;
	t4 = z - 1.0;
	t5 = 1.0 / mu;
	t6 = t*t3;
	t7 = omega*t4;
	t8 = cos(t7);
	t9 = 1.0 / t8;
	t10 = omega*t3;
	t11 = cos(t10);
	t12 = t*t4;
	t13 = lambda + mu;
	t14 = 1.0 / t13;
	t15 = lambda*2.0;
	t16 = mu*4.0;
	t17 = t15 + t16;
	t18 = t14*t17;
	t19 = t18 + z - 1.0;
	t20 = pow(r, t4);
	t21 = sin(t);
	t22 = sin(t6);
	t23 = t3*t22;
	t24 = sin(t12);
	t25 = t23 - t9*t11*t19*t24;
	t26 = pow(r, z);
	t27 = 1.0 / r;
	t28 = cos(t);
	t29 = cos(t6);
	t30 = t3*t3;
	t31 = t29*t30;
	t32 = cos(t12);
	t33 = t31 - t4*t9*t11*t19*t32;
	t34 = -t18 + z + 1.0;
	t35 = t3*t29;
	t36 = t35 - t9*t11*t32*t34;
	t37 = t22*t30;
	t38 = t37 - t4*t9*t11*t24*t34;
	t0 = t2*(t2*(t5*t21*t26*t27*t33*(1.0 / 2.0) - t5*t20*t25*t28*z*(1.0 / 2.0))*(1.0 / 2.0) - t2*(t5*t26*t27*t28*t33*(1.0 / 2.0) + t5*t20*t21*t25*z*(1.0 / 2.0))*(1.0 / 2.0))*(-1.0 / 2.0) - t2*(t2*(t5*t21*t26*t27*t38*(1.0 / 2.0) + t5*t20*t28*t36*z*(1.0 / 2.0))*(1.0 / 2.0) - t2*(t5*t26*t27*t28*t38*(1.0 / 2.0) - t5*t20*t21*t36*z*(1.0 / 2.0))*(1.0 / 2.0))*(1.0 / 2.0);

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double un_tt(double x, double y, double lambda, double mu, int edgeFlag)
* \brief the tangential-tangential derivative of u \dot n on boundary
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential-tangential derivative of u \dot n on boundary
*/
double un_tt(double x, double y, double lambda, double mu, int edgeFlag)
{
	/**for L-shaped domain **/
	double t0;
	switch (edgeFlag)
	{
	case 1: t0 = un_tt_1(x, y, lambda, mu); break;
	case 2: t0 = un_tt_2(x, y, lambda, mu); break;
	case 3: t0 = un_tt_3(x, y, lambda, mu); break;
	case 4: t0 = un_tt_4(x, y, lambda, mu); break;
	default: t0 = 0;
	}
	return t0;
}

/**
* \fn double un_tt_1(double x, double y, double lambda, double mu)
* \brief the tangential-tangential derivative of u \dot n on boundary type 1
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential-tangential derivative of u \dot n on boundary type 1
*/
double un_tt_1(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = z - 1.0;
	t3 = omega*t2;
	t4 = cos(t3);
	t5 = cos(t);
	t6 = z + 1.0;
	t7 = t5*t5;
	t8 = t*t6;
	t9 = cos(t8);
	t10 = omega*t6;
	t11 = cos(t10);
	t12 = t*t2;
	t13 = cos(t12);
	t14 = sin(t);
	t15 = t14*t14;
	t16 = sin(t8);
	t17 = sin(t12);
	t18 = z*z;
	t0 = (sqrt(2.0)*pow(r, z - 2.0)*(lambda*t4*t7*t9*3.0 + lambda*t4*t7*t16 - lambda*t4*t9*t15 + lambda*t7*t11*t13*3.0 - lambda*t4*t15*t16*3.0 - lambda*t7*t11*t17 - lambda*t11*t13*t15 + lambda*t11*t15*t17*3.0 + mu*t4*t7*t9*3.0 + mu*t4*t7*t16 - mu*t4*t9*t15 + mu*t7*t11*t13*9.0 - mu*t4*t15*t16*3.0 - mu*t7*t11*t17*3.0 - mu*t11*t13*t15*3.0 + mu*t11*t15*t17*9.0 + lambda*t4*t5*t9*t14*2.0 - lambda*t4*t5*t14*t16*2.0 + lambda*t5*t11*t13*t14*2.0 - lambda*t4*t7*t16*t18*4.0 + lambda*t4*t9*t15*t18*4.0 + lambda*t5*t11*t14*t17*2.0 - lambda*t7*t11*t17*t18*8.0 + lambda*t11*t13*t15*t18*4.0 + lambda*t11*t15*t17*t18*4.0 + mu*t4*t5*t9*t14*2.0 - mu*t4*t5*t14*t16*2.0 + mu*t5*t11*t13*t14*6.0 - mu*t4*t7*t16*t18*4.0 + mu*t4*t9*t15*t18*4.0 + mu*t5*t11*t14*t17*6.0 + mu*t7*t11*t13*t18*4.0 - mu*t7*t11*t17*t18*1.2E1 + mu*t11*t15*t17*t18*8.0 + lambda*t4*t7*t9*z*5.0 - lambda*t4*t7*t16*z + lambda*t4*t9*t15*z - lambda*t7*t11*t13*z*5.0 - lambda*t4*t15*t16*z*5.0 + lambda*t7*t11*t17*z*5.0 - lambda*t11*t13*t15*z - lambda*t11*t15*t17*z*7.0 + mu*t4*t7*t9*z*5.0 - mu*t4*t7*t16*z + mu*t4*t9*t15*z - mu*t7*t11*t13*z*1.7E1 - mu*t4*t15*t16*z*5.0 + mu*t7*t11*t17*z*9.0 + mu*t11*t13*t15*z*3.0 - mu*t11*t15*t17*z*1.9E1 + lambda*t4*t5*t9*t14*t18*4.0 - lambda*t4*t5*t14*t16*t18*4.0 + lambda*t5*t11*t13*t14*t18*1.2E1 - lambda*t5*t11*t14*t17*t18*4.0 + mu*t4*t5*t9*t14*t18*4.0 - mu*t4*t5*t14*t16*t18*4.0 + mu*t5*t11*t13*t14*t18*2.0E1 + mu*t5*t11*t14*t17*t18*4.0 + lambda*t4*t5*t9*t14*z*2.0 - lambda*t4*t7*t9*t18*z*2.0 - lambda*t4*t5*t14*t16*z*2.0 - lambda*t5*t11*t13*t14*z*1.0E1 - lambda*t4*t7*t16*t18*z*2.0 + lambda*t4*t9*t15*t18*z*2.0 - lambda*t5*t11*t14*t17*z*6.0 + lambda*t7*t11*t13*t18*z*2.0 + lambda*t4*t15*t16*t18*z*2.0 + lambda*t7*t11*t17*t18*z*2.0 - lambda*t11*t13*t15*t18*z*2.0 - lambda*t11*t15*t17*t18*z*2.0 + mu*t4*t5*t9*t14*z*2.0 - mu*t4*t7*t9*t18*z*2.0 - mu*t4*t5*t14*t16*z*2.0 - mu*t5*t11*t13*t14*z*2.6E1 - mu*t4*t7*t16*t18*z*2.0 + mu*t4*t9*t15*t18*z*2.0 - mu*t5*t11*t14*t17*z*2.2E1 + mu*t7*t11*t13*t18*z*2.0 + mu*t4*t15*t16*t18*z*2.0 + mu*t7*t11*t17*t18*z*2.0 - mu*t11*t13*t15*t18*z*2.0 - mu*t11*t15*t17*t18*z*2.0 + lambda*t4*t5*t9*t14*t18*z*4.0 - lambda*t4*t5*t14*t16*t18*z*4.0 - lambda*t5*t11*t13*t14*t18*z*4.0 + lambda*t5*t11*t14*t17*t18*z*4.0 + mu*t4*t5*t9*t14*t18*z*4.0 - mu*t4*t5*t14*t16*t18*z*4.0 - mu*t5*t11*t13*t14*t18*z*4.0 + mu*t5*t11*t14*t17*t18*z*4.0)*(-1.0 / 8.0)) / (mu*t4*(lambda + mu));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double un_tt_2(double x, double y, double lambda, double mu)
* \brief the tangential-tangential derivative of u \dot n on boundary type 2
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential-tangential derivative of u \dot n on boundary type 2
*/
double un_tt_2(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = z - 1.0;
	t3 = omega*t2;
	t4 = cos(t3);
	t5 = cos(t);
	t6 = z + 1.0;
	t7 = t5*t5;
	t8 = t*t6;
	t9 = cos(t8);
	t10 = omega*t6;
	t11 = cos(t10);
	t12 = t*t2;
	t13 = cos(t12);
	t14 = sin(t);
	t15 = t14*t14;
	t16 = sin(t8);
	t17 = sin(t12);
	t18 = z*z;
	t0 = (sqrt(2.0)*pow(r, z - 2.0)*(lambda*t4*t7*t9*-3.0 + lambda*t4*t7*t16 + lambda*t4*t9*t15 - lambda*t7*t11*t13*3.0 - lambda*t4*t15*t16*3.0 - lambda*t7*t11*t17 + lambda*t11*t13*t15 + lambda*t11*t15*t17*3.0 - mu*t4*t7*t9*3.0 + mu*t4*t7*t16 + mu*t4*t9*t15 - mu*t7*t11*t13*9.0 - mu*t4*t15*t16*3.0 - mu*t7*t11*t17*3.0 + mu*t11*t13*t15*3.0 + mu*t11*t15*t17*9.0 + lambda*t4*t5*t9*t14*2.0 + lambda*t4*t5*t14*t16*2.0 + lambda*t5*t11*t13*t14*2.0 - lambda*t4*t7*t16*t18*4.0 - lambda*t4*t9*t15*t18*4.0 - lambda*t5*t11*t14*t17*2.0 - lambda*t7*t11*t17*t18*8.0 - lambda*t11*t13*t15*t18*4.0 + lambda*t11*t15*t17*t18*4.0 + mu*t4*t5*t9*t14*2.0 + mu*t4*t5*t14*t16*2.0 + mu*t5*t11*t13*t14*6.0 - mu*t4*t7*t16*t18*4.0 - mu*t4*t9*t15*t18*4.0 - mu*t5*t11*t14*t17*6.0 - mu*t7*t11*t13*t18*4.0 - mu*t7*t11*t17*t18*1.2E1 + mu*t11*t15*t17*t18*8.0 - lambda*t4*t7*t9*z*5.0 - lambda*t4*t7*t16*z - lambda*t4*t9*t15*z + lambda*t7*t11*t13*z*5.0 - lambda*t4*t15*t16*z*5.0 + lambda*t7*t11*t17*z*5.0 + lambda*t11*t13*t15*z - lambda*t11*t15*t17*z*7.0 - mu*t4*t7*t9*z*5.0 - mu*t4*t7*t16*z - mu*t4*t9*t15*z + mu*t7*t11*t13*z*1.7E1 - mu*t4*t15*t16*z*5.0 + mu*t7*t11*t17*z*9.0 - mu*t11*t13*t15*z*3.0 - mu*t11*t15*t17*z*1.9E1 + lambda*t4*t5*t9*t14*t18*4.0 + lambda*t4*t5*t14*t16*t18*4.0 + lambda*t5*t11*t13*t14*t18*1.2E1 + lambda*t5*t11*t14*t17*t18*4.0 + mu*t4*t5*t9*t14*t18*4.0 + mu*t4*t5*t14*t16*t18*4.0 + mu*t5*t11*t13*t14*t18*2.0E1 - mu*t5*t11*t14*t17*t18*4.0 + lambda*t4*t5*t9*t14*z*2.0 + lambda*t4*t7*t9*t18*z*2.0 + lambda*t4*t5*t14*t16*z*2.0 - lambda*t5*t11*t13*t14*z*1.0E1 - lambda*t4*t7*t16*t18*z*2.0 - lambda*t4*t9*t15*t18*z*2.0 + lambda*t5*t11*t14*t17*z*6.0 - lambda*t7*t11*t13*t18*z*2.0 + lambda*t4*t15*t16*t18*z*2.0 + lambda*t7*t11*t17*t18*z*2.0 + lambda*t11*t13*t15*t18*z*2.0 - lambda*t11*t15*t17*t18*z*2.0 + mu*t4*t5*t9*t14*z*2.0 + mu*t4*t7*t9*t18*z*2.0 + mu*t4*t5*t14*t16*z*2.0 - mu*t5*t11*t13*t14*z*2.6E1 - mu*t4*t7*t16*t18*z*2.0 - mu*t4*t9*t15*t18*z*2.0 + mu*t5*t11*t14*t17*z*2.2E1 - mu*t7*t11*t13*t18*z*2.0 + mu*t4*t15*t16*t18*z*2.0 + mu*t7*t11*t17*t18*z*2.0 + mu*t11*t13*t15*t18*z*2.0 - mu*t11*t15*t17*t18*z*2.0 + lambda*t4*t5*t9*t14*t18*z*4.0 + lambda*t4*t5*t14*t16*t18*z*4.0 - lambda*t5*t11*t13*t14*t18*z*4.0 - lambda*t5*t11*t14*t17*t18*z*4.0 + mu*t4*t5*t9*t14*t18*z*4.0 + mu*t4*t5*t14*t16*t18*z*4.0 - mu*t5*t11*t13*t14*t18*z*4.0 - mu*t5*t11*t14*t17*t18*z*4.0)*(-1.0 / 8.0)) / (mu*t4*(lambda + mu));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double un_tt_3(double x, double y, double lambda, double mu)
* \brief the tangential-tangential derivative of u \dot n on boundary type 3
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential-tangential derivative of u \dot n on boundary type 3
*/
double un_tt_3(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = z - 1.0;
	t3 = omega*t2;
	t4 = cos(t3);
	t5 = cos(t);
	t6 = z + 1.0;
	t7 = t5*t5;
	t8 = t*t6;
	t9 = cos(t8);
	t10 = omega*t6;
	t11 = cos(t10);
	t12 = t*t2;
	t13 = cos(t12);
	t14 = sin(t);
	t15 = t14*t14;
	t16 = sin(t8);
	t17 = sin(t12);
	t18 = z*z;
	t0 = (sqrt(2.0)*pow(r, z - 2.0)*(lambda*t4*t7*t9*3.0 + lambda*t4*t7*t16 - lambda*t4*t9*t15 + lambda*t7*t11*t13*3.0 - lambda*t4*t15*t16*3.0 - lambda*t7*t11*t17 - lambda*t11*t13*t15 + lambda*t11*t15*t17*3.0 + mu*t4*t7*t9*3.0 + mu*t4*t7*t16 - mu*t4*t9*t15 + mu*t7*t11*t13*9.0 - mu*t4*t15*t16*3.0 - mu*t7*t11*t17*3.0 - mu*t11*t13*t15*3.0 + mu*t11*t15*t17*9.0 + lambda*t4*t5*t9*t14*2.0 - lambda*t4*t5*t14*t16*2.0 + lambda*t5*t11*t13*t14*2.0 - lambda*t4*t7*t16*t18*4.0 + lambda*t4*t9*t15*t18*4.0 + lambda*t5*t11*t14*t17*2.0 - lambda*t7*t11*t17*t18*8.0 + lambda*t11*t13*t15*t18*4.0 + lambda*t11*t15*t17*t18*4.0 + mu*t4*t5*t9*t14*2.0 - mu*t4*t5*t14*t16*2.0 + mu*t5*t11*t13*t14*6.0 - mu*t4*t7*t16*t18*4.0 + mu*t4*t9*t15*t18*4.0 + mu*t5*t11*t14*t17*6.0 + mu*t7*t11*t13*t18*4.0 - mu*t7*t11*t17*t18*1.2E1 + mu*t11*t15*t17*t18*8.0 + lambda*t4*t7*t9*z*5.0 - lambda*t4*t7*t16*z + lambda*t4*t9*t15*z - lambda*t7*t11*t13*z*5.0 - lambda*t4*t15*t16*z*5.0 + lambda*t7*t11*t17*z*5.0 - lambda*t11*t13*t15*z - lambda*t11*t15*t17*z*7.0 + mu*t4*t7*t9*z*5.0 - mu*t4*t7*t16*z + mu*t4*t9*t15*z - mu*t7*t11*t13*z*1.7E1 - mu*t4*t15*t16*z*5.0 + mu*t7*t11*t17*z*9.0 + mu*t11*t13*t15*z*3.0 - mu*t11*t15*t17*z*1.9E1 + lambda*t4*t5*t9*t14*t18*4.0 - lambda*t4*t5*t14*t16*t18*4.0 + lambda*t5*t11*t13*t14*t18*1.2E1 - lambda*t5*t11*t14*t17*t18*4.0 + mu*t4*t5*t9*t14*t18*4.0 - mu*t4*t5*t14*t16*t18*4.0 + mu*t5*t11*t13*t14*t18*2.0E1 + mu*t5*t11*t14*t17*t18*4.0 + lambda*t4*t5*t9*t14*z*2.0 - lambda*t4*t7*t9*t18*z*2.0 - lambda*t4*t5*t14*t16*z*2.0 - lambda*t5*t11*t13*t14*z*1.0E1 - lambda*t4*t7*t16*t18*z*2.0 + lambda*t4*t9*t15*t18*z*2.0 - lambda*t5*t11*t14*t17*z*6.0 + lambda*t7*t11*t13*t18*z*2.0 + lambda*t4*t15*t16*t18*z*2.0 + lambda*t7*t11*t17*t18*z*2.0 - lambda*t11*t13*t15*t18*z*2.0 - lambda*t11*t15*t17*t18*z*2.0 + mu*t4*t5*t9*t14*z*2.0 - mu*t4*t7*t9*t18*z*2.0 - mu*t4*t5*t14*t16*z*2.0 - mu*t5*t11*t13*t14*z*2.6E1 - mu*t4*t7*t16*t18*z*2.0 + mu*t4*t9*t15*t18*z*2.0 - mu*t5*t11*t14*t17*z*2.2E1 + mu*t7*t11*t13*t18*z*2.0 + mu*t4*t15*t16*t18*z*2.0 + mu*t7*t11*t17*t18*z*2.0 - mu*t11*t13*t15*t18*z*2.0 - mu*t11*t15*t17*t18*z*2.0 + lambda*t4*t5*t9*t14*t18*z*4.0 - lambda*t4*t5*t14*t16*t18*z*4.0 - lambda*t5*t11*t13*t14*t18*z*4.0 + lambda*t5*t11*t14*t17*t18*z*4.0 + mu*t4*t5*t9*t14*t18*z*4.0 - mu*t4*t5*t14*t16*t18*z*4.0 - mu*t5*t11*t13*t14*t18*z*4.0 + mu*t5*t11*t14*t17*t18*z*4.0)*(1.0 / 8.0)) / (mu*t4*(lambda + mu));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}

/**
* \fn double un_tt_4(double x, double y, double lambda, double mu)
* \brief the tangential-tangential derivative of u \dot n on boundary type 4
* \param x the x-axis value of the point
* \param y the y-axis value of the point
* \param lambda Lame constant
* \param mu Lame constant or Poisson ratio of plate
* \return the tangential-tangential derivative of u \dot n on boundary type 4
*/
double un_tt_4(double x, double y, double lambda, double mu)
{
	/**for L-shaped domain **/
	double r, t;
	cart2pol(x, y, &r, &t);

	double t0, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
	double t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
	double t31, t32, t33, t34, t35, t36, t37, t38, t39, t40;
	double t41, t42, t43, t44;
	t2 = z - 1.0;
	t3 = omega*t2;
	t4 = cos(t3);
	t5 = cos(t);
	t6 = z + 1.0;
	t7 = t5*t5;
	t8 = t*t6;
	t9 = cos(t8);
	t10 = omega*t6;
	t11 = cos(t10);
	t12 = t*t2;
	t13 = cos(t12);
	t14 = sin(t);
	t15 = t14*t14;
	t16 = sin(t8);
	t17 = sin(t12);
	t18 = z*z;
	t0 = (sqrt(2.0)*pow(r, z - 2.0)*(lambda*t4*t7*t9*-3.0 + lambda*t4*t7*t16 + lambda*t4*t9*t15 - lambda*t7*t11*t13*3.0 - lambda*t4*t15*t16*3.0 - lambda*t7*t11*t17 + lambda*t11*t13*t15 + lambda*t11*t15*t17*3.0 - mu*t4*t7*t9*3.0 + mu*t4*t7*t16 + mu*t4*t9*t15 - mu*t7*t11*t13*9.0 - mu*t4*t15*t16*3.0 - mu*t7*t11*t17*3.0 + mu*t11*t13*t15*3.0 + mu*t11*t15*t17*9.0 + lambda*t4*t5*t9*t14*2.0 + lambda*t4*t5*t14*t16*2.0 + lambda*t5*t11*t13*t14*2.0 - lambda*t4*t7*t16*t18*4.0 - lambda*t4*t9*t15*t18*4.0 - lambda*t5*t11*t14*t17*2.0 - lambda*t7*t11*t17*t18*8.0 - lambda*t11*t13*t15*t18*4.0 + lambda*t11*t15*t17*t18*4.0 + mu*t4*t5*t9*t14*2.0 + mu*t4*t5*t14*t16*2.0 + mu*t5*t11*t13*t14*6.0 - mu*t4*t7*t16*t18*4.0 - mu*t4*t9*t15*t18*4.0 - mu*t5*t11*t14*t17*6.0 - mu*t7*t11*t13*t18*4.0 - mu*t7*t11*t17*t18*1.2E1 + mu*t11*t15*t17*t18*8.0 - lambda*t4*t7*t9*z*5.0 - lambda*t4*t7*t16*z - lambda*t4*t9*t15*z + lambda*t7*t11*t13*z*5.0 - lambda*t4*t15*t16*z*5.0 + lambda*t7*t11*t17*z*5.0 + lambda*t11*t13*t15*z - lambda*t11*t15*t17*z*7.0 - mu*t4*t7*t9*z*5.0 - mu*t4*t7*t16*z - mu*t4*t9*t15*z + mu*t7*t11*t13*z*1.7E1 - mu*t4*t15*t16*z*5.0 + mu*t7*t11*t17*z*9.0 - mu*t11*t13*t15*z*3.0 - mu*t11*t15*t17*z*1.9E1 + lambda*t4*t5*t9*t14*t18*4.0 + lambda*t4*t5*t14*t16*t18*4.0 + lambda*t5*t11*t13*t14*t18*1.2E1 + lambda*t5*t11*t14*t17*t18*4.0 + mu*t4*t5*t9*t14*t18*4.0 + mu*t4*t5*t14*t16*t18*4.0 + mu*t5*t11*t13*t14*t18*2.0E1 - mu*t5*t11*t14*t17*t18*4.0 + lambda*t4*t5*t9*t14*z*2.0 + lambda*t4*t7*t9*t18*z*2.0 + lambda*t4*t5*t14*t16*z*2.0 - lambda*t5*t11*t13*t14*z*1.0E1 - lambda*t4*t7*t16*t18*z*2.0 - lambda*t4*t9*t15*t18*z*2.0 + lambda*t5*t11*t14*t17*z*6.0 - lambda*t7*t11*t13*t18*z*2.0 + lambda*t4*t15*t16*t18*z*2.0 + lambda*t7*t11*t17*t18*z*2.0 + lambda*t11*t13*t15*t18*z*2.0 - lambda*t11*t15*t17*t18*z*2.0 + mu*t4*t5*t9*t14*z*2.0 + mu*t4*t7*t9*t18*z*2.0 + mu*t4*t5*t14*t16*z*2.0 - mu*t5*t11*t13*t14*z*2.6E1 - mu*t4*t7*t16*t18*z*2.0 - mu*t4*t9*t15*t18*z*2.0 + mu*t5*t11*t14*t17*z*2.2E1 - mu*t7*t11*t13*t18*z*2.0 + mu*t4*t15*t16*t18*z*2.0 + mu*t7*t11*t17*t18*z*2.0 + mu*t11*t13*t15*t18*z*2.0 - mu*t11*t15*t17*t18*z*2.0 + lambda*t4*t5*t9*t14*t18*z*4.0 + lambda*t4*t5*t14*t16*t18*z*4.0 - lambda*t5*t11*t13*t14*t18*z*4.0 - lambda*t5*t11*t14*t17*t18*z*4.0 + mu*t4*t5*t9*t14*t18*z*4.0 + mu*t4*t5*t14*t16*t18*z*4.0 - mu*t5*t11*t13*t14*t18*z*4.0 - mu*t5*t11*t14*t17*t18*z*4.0)*(1.0 / 8.0)) / (mu*t4*(lambda + mu));

	return t0;
	//	return t0 / ((lambda + mu)*(lambda + mu));
}


/** 
 * \fn void morley_basis(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double *phi)
 * \brief basis function of Morley element
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param **nv the unit normal vectors of the three edges
 * \param **nve the unit normal vectors of the three edges (fixed for each edge)
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void morley_basis(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double *phi)
{
	if (index >= 6 || index<0)
	{
		*phi = 0;
		return;
	}

	int i, i1, i2;
	double sij[3], orient[3];
	for (i = 0; i < 3; i++)
	{
		sij[i] = -(xi[(i + 1) % 3] * xi[(i + 2) % 3] + eta[(i + 1) % 3] * eta[(i + 2) % 3]);
		orient[i] = nv[i][0] * nve[i][0]+ nv[i][1] * nve[i][1];
	}

	if (index<3)
	{
		i = index;
		i1 = (i + 1) % 3;
		i2 = (i + 2) % 3;
		*phi = lambda[i] * lambda[i] + sij[i2] / (elen[i1] * elen[i1])*lambda[i1] * (1 - lambda[i1]) + sij[i1] / (elen[i2] * elen[i2])*lambda[i2] * (1 - lambda[i2]);
	}
	else
	{ 
		i = index - 3;
		*phi = orient[i] * 2 * s / elen[i] * lambda[i] * (lambda[i] - 1);
	}
}

/** 
 * \fn void morley_basis1(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double phi[2])
 * \brief the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param **nv the unit normal vectors of the three edges
 * \param **nve the unit normal vectors of the three edges (fixed for each edge)
 * \param index the indicator of the basis function
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \return void
 */
void morley_basis1(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double phi[2])
{
	if (index >= 6 || index<0)
	{
		phi[0]=0;
		phi[1]=0;
		return;
	}

	int i, i1, i2;
	double sij[3], orient[3];
	for (i = 0; i < 3; i++)
	{
		sij[i] = -(xi[(i + 1) % 3] * xi[(i + 2) % 3] + eta[(i + 1) % 3] * eta[(i + 2) % 3]);
		orient[i] = nv[i][0] * nve[i][0] + nv[i][1] * nve[i][1];
	}

	if (index<3)
	{
		i = index;
		i1 = (i + 1) % 3;
		i2 = (i + 2) % 3;
		phi[0] = (2 * lambda[i] * eta[i] + sij[i2] / (elen[i1] * elen[i1])* (1 - 2 * lambda[i1])*eta[i1] + sij[i1] / (elen[i2] * elen[i2])* (1 - 2 * lambda[i2])*eta[i2]) / (2 * s);
		phi[1] = -(2 * lambda[i] * xi[i] + sij[i2] / (elen[i1] * elen[i1])* (1 - 2 * lambda[i1])*xi[i1] + sij[i1] / (elen[i2] * elen[i2])* (1 - 2 * lambda[i2])*xi[i2]) / (2 * s);
	}
	else
	{
		i = index - 3;
		phi[0] = orient[i] / elen[i] * (2 * lambda[i] - 1) * eta[i];
		phi[1] = -orient[i] / elen[i] * (2 * lambda[i] - 1) * xi[i];
	}
}

/** 
 * \fn void morley_basis2(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double phi[3])
 * \brief the second order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param **nv the unit normal vectors of the three edges
 * \param **nve the unit normal vectors of the three edges (fixed for each edge)
 * \param index the indicator of the basis function
 * \param phi[3] the second order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \return void
 */
void morley_basis2(double *lambda, double s, double elen[3], double eta[3], double xi[3], double **nv, double **nve, int index, double phi[3])
{
	if (index >= 6 || index<0)
	{
		phi[0] = 0;
		phi[1] = 0;
		phi[2] = 0;
		return;
	}

	int i, i1, i2;
	double sij[3], orient[3];
	for (i = 0; i < 3; i++)
	{
		sij[i] = -(xi[(i + 1) % 3] * xi[(i + 2) % 3] + eta[(i + 1) % 3] * eta[(i + 2) % 3]);
		orient[i] = nv[i][0] * nve[i][0] + nv[i][1] * nve[i][1];
	}
	
	if (index<3)
	{
		i = index;
		i1 = (i + 1) % 3;
		i2 = (i + 2) % 3;
		phi[0] = (eta[i] * eta[i] - sij[i2] / (elen[i1] * elen[i1])*eta[i1] * eta[i1] - sij[i1] / (elen[i2] * elen[i2])*eta[i2] * eta[i2]) / (2 * s*s);
		phi[1] = (xi[i] * xi[i] - sij[i2] / (elen[i1] * elen[i1])*xi[i1] * xi[i1] - sij[i1] / (elen[i2] * elen[i2])*xi[i2] * xi[i2]) / (2 * s*s);
		phi[2] = -(xi[i] * eta[i] - sij[i2] / (elen[i1] * elen[i1])*xi[i1] * eta[i1] - sij[i1] / (elen[i2] * elen[i2])*xi[i2] * eta[i2]) / (2 * s*s);
	}
	else
	{
		i = index - 3;
		phi[0] = eta[i] * eta[i] / (s*elen[i])*orient[i];
		phi[1] = xi[i] * xi[i] / (s*elen[i])*orient[i];
		phi[2] = -xi[i] * eta[i] / (s*elen[i])*orient[i];
	}
}

/** 
 * \fn void lagrange1D_basis(double lambda, int index, int dop, double *phi)
 * \brief basis function of Lagrange element in 1d space
 * \param lambda  area coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void lagrange1D_basis(double lambda, int index, int dop, double *phi)
{
	int in, ie, ii;
	int dofs = dop+1; // degrees of freedom
	if(index>= dofs || index<0)
	{
		*phi=0;
		return;
	}

	if(dop==0)
	{
		*phi = 1;
	}

	else if(dop==1)
	{
		switch(index)
		{
		case 0: *phi = lambda; break;
		case 1: *phi = 1-lambda; break;
		default: *phi = 0;
		}
	} // dop=1
	
	else if(dop==2)
	{
		switch(index)
		{
		case 0: *phi = lambda*(2*lambda-1); break;
		case 1: *phi = (1-lambda)*(1-2*lambda); break;
		case 2: *phi = 4*lambda*(1-lambda); break;
		default: *phi = 0;
		}
	} // dop=2
	
	else if(dop==3)
	{
		switch(index)
		{
		case 0: *phi = lambda*(3*lambda-1)*(3*lambda-2)/2.0; break;
		case 1: *phi = -(3*lambda-1)*(3*lambda-2)*(lambda-1)/2.0; break;
		case 2: *phi = -9*lambda*(3*lambda-1)*(lambda-1)/2.0; break;
		case 3: *phi = 9*lambda*(3*lambda-2)*(lambda-1)/2.0; break;
		default: *phi = 0;
		}
	} // dop=3

	else if(dop==4)
	{
		switch(index)
		{
		case 0: *phi = lambda*(4*lambda-1)*(2*lambda-1)*(4*lambda-3)/3.0; break;
		case 1: *phi = (4*lambda-1)*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0; break;
		case 2: *phi = -16*lambda*(4*lambda-1)*(2*lambda-1)*(lambda-1)/3.0; break;
		case 3: *phi = 4*lambda*(4*lambda-1)*(4*lambda-3)*(lambda-1); break;
		case 4: *phi = -16*lambda*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0; break;
		default: *phi = 0;
		}
	} // dop=4

	else if(dop==5)
	{
		switch(index)
		{
		case 0: *phi = lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)/24.0; break;
		case 1: *phi = -(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0; break;
		case 2: *phi = -25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(lambda-1)/24.0; break;
		case 3: *phi = 25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-4)*(lambda-1)/12.0; break;
		case 4: *phi = -25*lambda*(5*lambda-1)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/12.0; break;
		case 5: *phi = 25*lambda*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0; break;
		default: *phi = 0;
		}
	} // dop=5
}

/** 
 * \fn void lagrange1D_basis1(double lambda, int index, int dop, double h, double *phi)
 * \brief the first order derivative of basis function of Lagrange element in 1d space
 * \param lambda  area coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param h length of edge
 * \param *phi basis function
 * \return void
 */
void lagrange1D_basis1(double lambda, int index, int dop, double h, double *phi)
{
	int in, ie, ii;
	int dofs = dop+1; // degrees of freedom
	if(index>= dofs || index<0)
	{
		*phi=0;
		return;
	}

	double slp=-1.0/h; // slope

	if(dop==0)
	{
		*phi = 0;
	}

	else if(dop==1)
	{
		switch(index)
		{
		case 0: *phi = slp; break;
		case 1: *phi = -slp; break;
		default: *phi = 0;
		}
	} // dop=1
	
	else if(dop==2)
	{
		switch(index)
		{
		case 0: *phi = slp*(2*lambda-1) + lambda*2*slp; break;
		case 1: *phi = -slp*(1-2*lambda) - (1-lambda)*2*slp; break;
		case 2: *phi = 4*slp*(1-lambda) - 4*lambda*slp; break;
		default: *phi = 0;
		}
	} // dop=2
	
	else if(dop==3)
	{
		switch(index)
		{
		case 0: *phi = slp*(3*lambda-1)*(3*lambda-2)/2.0 + lambda*3*slp*(3*lambda-2)/2.0 + lambda*(3*lambda-1)*3*slp/2.0; break;
		case 1: *phi = -3*slp*(3*lambda-2)*(lambda-1)/2.0 - (3*lambda-1)*3*slp*(lambda-1)/2.0 - (3*lambda-1)*(3*lambda-2)*slp/2.0; break;
		case 2: *phi = -9*slp*(3*lambda-1)*(lambda-1)/2.0 - 9*lambda*3*slp*(lambda-1)/2.0 - 9*lambda*(3*lambda-1)*slp/2.0; break;
		case 3: *phi = 9*slp*(3*lambda-2)*(lambda-1)/2.0 + 9*lambda*3*slp*(lambda-1)/2.0 + 9*lambda*(3*lambda-2)*slp/2.0; break;
		default: *phi = 0;
		}
	} // dop=3

	else if(dop==4)
	{
		switch(index)
		{
		case 0: *phi = slp*(4*lambda-1)*(2*lambda-1)*(4*lambda-3)/3.0 + lambda*4*slp*(2*lambda-1)*(4*lambda-3)/3.0 + lambda*(4*lambda-1)*2*slp*(4*lambda-3)/3.0 + lambda*(4*lambda-1)*(2*lambda-1)*4*slp/3.0; break;
		case 1: *phi = 4*slp*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0 + (4*lambda-1)*2*slp*(4*lambda-3)*(lambda-1)/3.0 + (4*lambda-1)*(2*lambda-1)*4*slp*(lambda-1)/3.0 + (4*lambda-1)*(2*lambda-1)*(4*lambda-3)*slp/3.0; break;
		case 2: *phi = -16*slp*(4*lambda-1)*(2*lambda-1)*(lambda-1)/3.0 - 16*lambda*4*slp*(2*lambda-1)*(lambda-1)/3.0 - 16*lambda*(4*lambda-1)*2*slp*(lambda-1)/3.0 - 16*lambda*(4*lambda-1)*(2*lambda-1)*slp/3.0; break;
		case 3: *phi = 4*slp*(4*lambda-1)*(4*lambda-3)*(lambda-1) + 4*lambda*4*slp*(4*lambda-3)*(lambda-1) + 4*lambda*(4*lambda-1)*4*slp*(lambda-1) + 4*lambda*(4*lambda-1)*(4*lambda-3)*slp; break;
		case 4: *phi = -16*slp*(2*lambda-1)*(4*lambda-3)*(lambda-1)/3.0 - 16*lambda*2*slp*(4*lambda-3)*(lambda-1)/3.0 - 16*lambda*(2*lambda-1)*4*slp*(lambda-1)/3.0 - 16*lambda*(2*lambda-1)*(4*lambda-3)*slp/3.0; break;
		default: *phi = 0;
		}
	} // dop=4

	else if(dop==5)
	{
		switch(index)
		{
		case 0: *phi = slp*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)/24.0 + lambda*5*slp*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)/24.0 + lambda*(5*lambda-1)*5*slp*(5*lambda-3)*(5*lambda-4)/24.0 + lambda*(5*lambda-1)*(5*lambda-2)*5*slp*(5*lambda-4)/24.0 + lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*5*slp/24.0; break;
		case 1: *phi = -5*slp*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 - (5*lambda-1)*5*slp*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 - (5*lambda-1)*(5*lambda-2)*5*slp*(5*lambda-4)*(lambda-1)/24.0 - (5*lambda-1)*(5*lambda-2)*(5*lambda-3)*5*slp*(lambda-1)/24.0 - (5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*slp/24.0; break;
		case 2: *phi = -25*slp*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*(lambda-1)/24.0 - 25*lambda*5*slp*(5*lambda-2)*(5*lambda-3)*(lambda-1)/24.0 - 25*lambda*(5*lambda-1)*5*slp*(5*lambda-3)*(lambda-1)/24.0 - 25*lambda*(5*lambda-1)*(5*lambda-2)*5*slp*(lambda-1)/24.0 - 25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-3)*slp/24.0; break;
		case 3: *phi = 25*slp*(5*lambda-1)*(5*lambda-2)*(5*lambda-4)*(lambda-1)/12.0 + 25*lambda*5*slp*(5*lambda-2)*(5*lambda-4)*(lambda-1)/12.0 + 25*lambda*(5*lambda-1)*5*slp*(5*lambda-4)*(lambda-1)/12.0 + 25*lambda*(5*lambda-1)*(5*lambda-2)*5*slp*(lambda-1)/12.0 + 25*lambda*(5*lambda-1)*(5*lambda-2)*(5*lambda-4)*slp/12.0; break;
		case 4: *phi = -25*slp*(5*lambda-1)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/12.0 - 25*lambda*5*slp*(5*lambda-3)*(5*lambda-4)*(lambda-1)/12.0 - 25*lambda*(5*lambda-1)*5*slp*(5*lambda-4)*(lambda-1)/12.0 - 25*lambda*(5*lambda-1)*(5*lambda-3)*5*slp*(lambda-1)/12.0 - 25*lambda*(5*lambda-1)*(5*lambda-3)*(5*lambda-4)*slp/12.0; break;
		case 5: *phi = 25*slp*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 + 25*lambda*5*slp*(5*lambda-3)*(5*lambda-4)*(lambda-1)/24.0 + 25*lambda*(5*lambda-2)*5*slp*(5*lambda-4)*(lambda-1)/24.0 + 25*lambda*(5*lambda-2)*(5*lambda-3)*5*slp*(lambda-1)/24.0 + 25*lambda*(5*lambda-2)*(5*lambda-3)*(5*lambda-4)*slp/24.0; break;
		default: *phi = 0;
		}
	} // dop=5
}

/** 
 * \fn void lagrange_basis(double *lambda, int index, int dop, double *phi)
 * \brief basis function of Lagrange element
 * \param *lambda pointer to the area coordiante
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void lagrange_basis(double *lambda, int index, int dop, double *phi)
{
	int in, ie, ii;
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0)
	{
		*phi=0;
		return;
	}

	if(dop==0)
	{
		*phi = 1;
	}

	else if(dop==1)
	{
		*phi = lambda[index];
	} // dop=1
	
	else if(dop==2)
	{
		if(index<3)
		{
			*phi = lambda[index]*(2*lambda[index]-1);
		}
		else
		{
			*phi = 4.0*lambda[(index+1)%3]*lambda[(index+2)%3];
		}
	} // dop=2
	
	else if(dop==3)
	{
		if(index<3)
		{
			*phi = lambda[index]*(3*lambda[index]-1)*(3*lambda[index]-2)/2.0;
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(3*lambda[(ie+1+ii)%3]-1)*9.0/2.0;
		}
		else
		{
			*phi = 27.0*lambda[0]*lambda[1]*lambda[2];
		}
	} // dop=3

	else if(dop==4)
	{
		if(index<3)
		{
			*phi = lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3)/6.0;
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			switch(ii)
			{
			case 0:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(4*lambda[(ie+1)%3]-1)*(4*lambda[(ie+1)%3]-2)*8.0/3.0;
				break; 
			case 1:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(4*lambda[(ie+1)%3]-1)*(4*lambda[(ie+2)%3]-1)*4.0;
				break;
			case 2:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(4*lambda[(ie+2)%3]-1)*(4*lambda[(ie+2)%3]-2)*8.0/3.0;
				break;
			default:
				*phi = 0;
			}
		}
		else
		{
			in = index-3*dop;
			*phi = 32.0*lambda[0]*lambda[1]*lambda[2]*(4*lambda[in]-1);
		}
	} // dop=4

	else if(dop==5)
	{
		if(index<3)
		{
			*phi = lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4)/24.0;
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			switch(ii)
			{
			case 0:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+1)%3]-1)*(5*lambda[(ie+1)%3]-2)*(5*lambda[(ie+1)%3]-3)*25.0/24.0;
				break; 
			case 1:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+1)%3]-1)*(5*lambda[(ie+2)%3]-1)*(5*lambda[(ie+1)%3]-2)*25.0/12.0;
				break;
			case 2:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+1)%3]-1)*(5*lambda[(ie+2)%3]-1)*(5*lambda[(ie+2)%3]-2)*25.0/12.0;
				break;
			case 3:
				*phi = lambda[(ie+1)%3]*lambda[(ie+2)%3]*(5*lambda[(ie+2)%3]-1)*(5*lambda[(ie+2)%3]-2)*(5*lambda[(ie+2)%3]-3)*25.0/24.0;
				break;
			default:
				*phi = 0;
			}
		}
		else if(index < 3*dop+3)
		{
			in = index-3*dop;
			*phi = lambda[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2)*125.0/6.0;
		}
		else
		{
			in = index-3*dop-3;
			*phi = lambda[0]*lambda[1]*lambda[2]*(5*lambda[(in+1)%3]-1)*(5*lambda[(in+2)%3]-1)*125.0/4.0;
		}
	} // dop=5
}

/** 
 * \fn void lagrange_basis1(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[2])
 * \brief the first order derivative of Lagrange element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
 * \return void
 */
void lagrange_basis1(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[2])
{
	int in, ie, ii, i1, i2, i3;
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0)
	{
		phi[0]=0;
		phi[1]=0;
		return;
	}

	if(dop==0)
	{
		phi[0]=0;
		phi[1]=0;
	} // dop=0

	else if(dop==1)
	{
		phi[0]=eta[index]/(2.0*s);
		phi[1]=-xi[index]/(2.0*s);
	} // dop=1

	else if(dop==2)
	{
		if(index<3)
		{
			phi[0]=(4*lambda[index]-1)*eta[index]/(2.0*s);
			phi[1]=-(4*lambda[index]-1)*xi[index]/(2.0*s);
		}
		else
		{
			i1=(index+1)%3;
			i2=(index+2)%3;
			phi[0]=4*(lambda[i1]*eta[i2]+lambda[i2]*eta[i1])/(2.0*s);
			phi[1]=-4*(lambda[i1]*xi[i2]+lambda[i2]*xi[i1])/(2.0*s);
		}
	} // dop=2

	else if(dop==3)
	{
		if(index<3)
		{
			phi[0]=((3*lambda[index]-1)*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-1))*eta[index]/(4.0*s);
			phi[1]=-((3*lambda[index]-1)*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-2) + 3*lambda[index]*(3*lambda[index]-1))*xi[index]/(4.0*s);
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			i3 = (ie+1+ii)%3;
			phi[0]=(lambda[i2]*(3*lambda[i3]-1)*eta[i1] + lambda[i1]*(3*lambda[i3]-1)*eta[i2] + 3*lambda[i1]*lambda[i2]*eta[i3])*9.0/(4.0*s);
			phi[1]=-(lambda[i2]*(3*lambda[i3]-1)*xi[i1] + lambda[i1]*(3*lambda[i3]-1)*xi[i2] + 3*lambda[i1]*lambda[i2]*xi[i3])*9.0/(4.0*s);
		}
		else
		{
			phi[0]=27*(lambda[1]*lambda[2]*eta[0] + lambda[2]*lambda[0]*eta[1] + lambda[0]*lambda[1]*eta[2])/(2.0*s);
			phi[1]=-27*(lambda[1]*lambda[2]*xi[0] + lambda[2]*lambda[0]*xi[1] + lambda[0]*lambda[1]*xi[2])/(2.0*s);
		}
	} // dop=3

	else if(dop==4)
	{
		if(index<3)
		{
			phi[0]=((4*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-2))*eta[index]/(12.0*s);
			phi[1]=-((4*lambda[index]-1)*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-2)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-3) + 4*lambda[index]*(4*lambda[index]-1)*(4*lambda[index]-2))*xi[index]/(12.0*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]=(eta[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*eta[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*4*eta[i1]*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*eta[i1])*4.0/(3.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*xi[i2]*(4*lambda[i1]-1)*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*4*xi[i1]*(4*lambda[i1]-2) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*xi[i1])*4.0/(3.0*s);
				break; 
			case 1:
				phi[0]=(eta[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*eta[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*4*eta[i1]*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*eta[i2])*2.0/s;
				phi[1]=-(xi[i1]*lambda[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*xi[i2]*(4*lambda[i1]-1)*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*4*xi[i1]*(4*lambda[i2]-1) + lambda[i1]*lambda[i2]*(4*lambda[i1]-1)*4*xi[i2])*2.0/s;
				break;
			case 2:
				phi[0]=(eta[i1]*lambda[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*eta[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*4*eta[i2]*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*(4*lambda[i2]-1)*4*eta[i2])*4.0/(3.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*xi[i2]*(4*lambda[i2]-1)*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*4*xi[i2]*(4*lambda[i2]-2) + lambda[i1]*lambda[i2]*(4*lambda[i2]-1)*4*xi[i2])*4.0/(3.0*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
			}
		}
		else
		{
			in = index-3*dop;
			phi[0]=(eta[0]*lambda[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*eta[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*eta[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*lambda[2]*4*eta[in])*16.0/s;
			phi[1]=-(xi[0]*lambda[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*xi[1]*lambda[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*xi[2]*(4*lambda[in]-1) + lambda[0]*lambda[1]*lambda[2]*4*xi[in])*16.0/s;
		}
	} // dop=4

	else if(dop==5)
	{
		if(index<3)
		{
			phi[0]=(eta[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*5*eta[index]*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*5*eta[index]*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*5*eta[index]*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*5*eta[index])/(48.0*s);
			phi[1]=-(xi[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*5*xi[index]*(5*lambda[index]-2)*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*5*xi[index]*(5*lambda[index]-3)*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*5*xi[index]*(5*lambda[index]-4) + lambda[index]*(5*lambda[index]-1)*(5*lambda[index]-2)*(5*lambda[index]-3)*5*xi[index])/(48.0*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*eta[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*5*eta[i1]*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*eta[i1]*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*5*eta[i1])*25.0/(48.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*xi[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*5*xi[i1]*(5*lambda[i1]-2)*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*xi[i1]*(5*lambda[i1]-3) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i1]-2)*5*xi[i1])*25.0/(48.0*s);
				break; 
			case 1:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*eta[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*5*eta[i1]*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*eta[i2]*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*eta[i1])*25.0/(24.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*xi[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*5*xi[i1]*(5*lambda[i2]-1)*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*xi[i2]*(5*lambda[i1]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*xi[i1])*25.0/(24.0*s);
				break;
			case 2:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*eta[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*5*eta[i1]*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*eta[i2]*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*eta[i2])*25.0/(24.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*xi[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*5*xi[i1]*(5*lambda[i2]-1)*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*5*xi[i2]*(5*lambda[i2]-2) + lambda[i1]*lambda[i2]*(5*lambda[i1]-1)*(5*lambda[i2]-1)*5*xi[i2])*25.0/(24.0*s);
				break;
			case 3:
				phi[0]=(eta[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*eta[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*5*eta[i2]*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*5*eta[i2]*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*5*eta[i2])*25.0/(48.0*s);
				phi[1]=-(xi[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*xi[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*5*xi[i2]*(5*lambda[i2]-2)*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*5*xi[i2]*(5*lambda[i2]-3) + lambda[i1]*lambda[i2]*(5*lambda[i2]-1)*(5*lambda[i2]-2)*5*xi[i2])*25.0/(48.0*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
			}
		}
		else if(index < 3*dop+3)
		{
			in = index-3*dop;
			phi[0]=(eta[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*eta[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*eta[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*5*eta[in]*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*5*eta[in])*125.0/(12.0*s);
			phi[1]=-(xi[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*xi[1]*lambda[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*xi[2]*(5*lambda[in]-1)*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*5*xi[in]*(5*lambda[in]-2) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[in]-1)*5*xi[in])*125.0/(12.0*s);
		}
		else
		{
			in = index-3*dop-3;
			i1 = (in+1)%3;
			i2 = (in+2)%3;
			phi[0]=(eta[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*eta[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*eta[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*5*eta[i1]*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*5*eta[i2])*125.0/(8.0*s);
			phi[1]=-(xi[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*xi[1]*lambda[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*xi[2]*(5*lambda[i1]-1)*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*5*xi[i1]*(5*lambda[i2]-1) + lambda[0]*lambda[1]*lambda[2]*(5*lambda[i1]-1)*5*xi[i2])*125.0/(8.0*s);
		}
	} // dop=5
}

/** 
 * \fn void lagrange_basis2(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[3])
 * \brief the second order derivative of Lagrange element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] the first order derivative of Morley element basis function: (\partial_{xx}phi, \partial_{yy}phi, \partial_{xy}phi)
 * \return void
 */
void lagrange_basis2(double *lambda, double s, double eta[3], double xi[3], int index, int dop, double phi[3])
{
	int in, ie, ii, i1, i2, i3;
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom
	if(index>= dofs || index<0)
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
		return;
	}

	if(dop==0)
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
	}

	if(dop==1)
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
	}

	else if(dop==2)
	{
		if(index<3)
		{
			phi[0]=eta[index]*eta[index]/(s*s);
			phi[1]=xi[index]*xi[index]/(s*s);
			phi[2]=-eta[index]*xi[index]/(s*s);
		}
		else
		{
			i1=(index+1)%3;
			i2=(index+2)%3;
			phi[0]=2*eta[i1]*eta[i2]/(s*s);
			phi[1]=2*xi[i1]*xi[i2]/(s*s);
			phi[2]=-(eta[i1]*xi[i2]+eta[i2]*xi[i1])/(s*s);
		}
	} // dop=2

	else if(dop==3)
	{
		if(index<3)
		{
			phi[0]=9*(3*lambda[index]-1)*eta[index]*eta[index]/(4*s*s);
			phi[1]=9*(3*lambda[index]-1)*xi[index]*xi[index]/(4*s*s);
			phi[2]=-9*(3*lambda[index]-1)*eta[index]*xi[index]/(4*s*s);
		}
		else if(index < 3*dop) 
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			i3 = (ie+1+ii)%3;
			phi[0]=9*(3*lambda[i1]*eta[i2]*eta[i3]+3*lambda[i2]*eta[i3]*eta[i1]+(3*lambda[i3]-1)*eta[i1]*eta[i2])/(4*s*s);
			phi[1]=9*(3*lambda[i1]*xi[i2]*xi[i3]+3*lambda[i2]*xi[i3]*xi[i1]+(3*lambda[i3]-1)*xi[i1]*xi[i2])/(4*s*s);
			phi[2]=-9*(3*lambda[i1]*(eta[i2]*xi[i3]+eta[i3]*xi[i2]) + 3*lambda[i2]*(eta[i3]*xi[i1]+eta[i1]*xi[i3]) + (3*lambda[i3]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]))/(8*s*s);
		}
		else
		{
			phi[0]=27*(lambda[0]*eta[1]*eta[2]+lambda[1]*eta[2]*eta[0]+lambda[2]*eta[0]*eta[1])/(2*s*s);
			phi[1]=27*(lambda[0]*xi[1]*xi[2]+lambda[1]*xi[2]*xi[0]+lambda[2]*xi[0]*xi[1])/(2*s*s);
			phi[2]=-27*(lambda[0]*(eta[1]*xi[2]+eta[2]*xi[1])+lambda[1]*(eta[2]*xi[0]+eta[0]*xi[2])+lambda[2]*(eta[0]*xi[1]+eta[1]*xi[0]))/(4*s*s);
		}
	} // dop=3

	else if(dop==4)
	{
		if(index<3)
		{
			phi[0]=(128*lambda[index]*lambda[index]-96*lambda[index]+44.0/3.0)*eta[index]*eta[index]/(4*s*s);
			phi[1]=(128*lambda[index]*lambda[index]-96*lambda[index]+44.0/3.0)*xi[index]*xi[index]/(4*s*s);
			phi[2]=-(128*lambda[index]*lambda[index]-96*lambda[index]+44.0/3.0)*eta[index]*xi[index]/(4*s*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]=(32.0/3.0*(24*lambda[i1]*lambda[i1]-12*lambda[i1]+1)*eta[i1]*eta[i2] + 64*lambda[i2]*(4*lambda[i1]-1)*eta[i1]*eta[i1])/(4*s*s);
				phi[1]=(32.0/3.0*(24*lambda[i1]*lambda[i1]-12*lambda[i1]+1)*xi[i1]*xi[i2] + 64*lambda[i2]*(4*lambda[i1]-1)*xi[i1]*xi[i1])/(4*s*s);
				phi[2]=-(16.0/3.0*(24*lambda[i1]*lambda[i1]-12*lambda[i1]+1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + 64*lambda[i2]*(4*lambda[i1]-1)*eta[i1]*xi[i1])/(4*s*s);
				break; 
			case 1:
				phi[0]=(32*(4*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*eta[i1] + 8*(8*lambda[i1]-1)*(8*lambda[i2]-1)*eta[i1]*eta[i2] + 32*(4*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*eta[i2])/(4*s*s);
				phi[1]=(32*(4*lambda[i2]*lambda[i2]-lambda[i2])*xi[i1]*xi[i1] + 8*(8*lambda[i1]-1)*(8*lambda[i2]-1)*xi[i1]*xi[i2] + 32*(4*lambda[i1]*lambda[i1]-lambda[i1])*xi[i2]*xi[i2])/(4*s*s);
				phi[2]=-(32*(4*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*xi[i1] + 4*(8*lambda[i1]-1)*(8*lambda[i2]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + 32*(4*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*xi[i2])/(4*s*s);
				break;
			case 2:
				phi[0]=(32.0/3.0*(24*lambda[i2]*lambda[i2]-12*lambda[i2]+1)*eta[i1]*eta[i2] + 64*lambda[i1]*(4*lambda[i2]-1)*eta[i2]*eta[i2])/(4*s*s);
				phi[1]=(32.0/3.0*(24*lambda[i2]*lambda[i2]-12*lambda[i2]+1)*xi[i1]*xi[i2] + 64*lambda[i1]*(4*lambda[i2]-1)*xi[i2]*xi[i2])/(4*s*s);
				phi[2]=-(16.0/3.0*(24*lambda[i2]*lambda[i2]-12*lambda[i2]+1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + 64*lambda[i1]*(4*lambda[i2]-1)*eta[i2]*xi[i2])/(4*s*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
				phi[2] = 0;
			}
		}
		else
		{
			in = index-3*dop;
			phi[0]=(eta[1]*eta[2]*lambda[0]*(4*lambda[in]-1) + eta[2]*eta[0]*lambda[1]*(4*lambda[in]-1) + eta[0]*eta[1]*lambda[2]*(4*lambda[in]-1) + 4*eta[in]*eta[0]*lambda[1]*lambda[2] + 4*eta[in]*eta[1]*lambda[2]*lambda[0] + 4*eta[in]*eta[2]*lambda[0]*lambda[1])*16.0/(s*s);
			phi[1]=(xi[1]*xi[2]*lambda[0]*(4*lambda[in]-1) + xi[2]*xi[0]*lambda[1]*(4*lambda[in]-1) + xi[0]*xi[1]*lambda[2]*(4*lambda[in]-1) + 4*xi[in]*xi[0]*lambda[1]*lambda[2] + 4*xi[in]*xi[1]*lambda[2]*lambda[0] + 4*xi[in]*xi[2]*lambda[0]*lambda[1])*16.0/(s*s);
			phi[2]=-((eta[1]*xi[2]+eta[2]*xi[1])*lambda[0]*(4*lambda[in]-1) + (eta[2]*xi[0]+eta[0]*xi[2])*lambda[1]*(4*lambda[in]-1) + (eta[0]*xi[1]+eta[1]*xi[0])*lambda[2]*(4*lambda[in]-1) + 4*(eta[in]*xi[0]+eta[0]*xi[in])*lambda[1]*lambda[2] + 4*(eta[in]*xi[1]+eta[1]*xi[in])*lambda[2]*lambda[0] + 4*(eta[in]*xi[2]+eta[2]*xi[in])*lambda[0]*lambda[1])*8.0/(s*s);
		}
	} // dop=4

	else if(dop==5)
	{
		if(index<3)
		{
			phi[0]=(3125.0/6.0*lambda[index]*lambda[index]*lambda[index]-625*lambda[index]*lambda[index]+875.0/4.0*lambda[index]-125.0/6.0)*eta[index]*eta[index]/(4*s*s);
			phi[1]=(3125.0/6.0*lambda[index]*lambda[index]*lambda[index]-625*lambda[index]*lambda[index]+875.0/4.0*lambda[index]-125.0/6.0)*xi[index]*xi[index]/(4*s*s);
			phi[2]=-(3125.0/6.0*lambda[index]*lambda[index]*lambda[index]-625*lambda[index]*lambda[index]+875.0/4.0*lambda[index]-125.0/6.0)*eta[index]*xi[index]/(4*s*s);
		}
		else if(index < 3*dop)
		{
			in = index-3;
			ie = in/(dop-1);
			ii = in%(dop-1);
			i1 = (ie+1)%3;
			i2 = (ie+2)%3;
			switch(ii)
			{
			case 0:
				phi[0]= 25.0/12.0*((500*lambda[i1]*lambda[i1]*lambda[i1]-450*lambda[i1]*lambda[i1]+110*lambda[i1]-6)*eta[i1]*eta[i2] + (750*lambda[i1]*lambda[i1]-450*lambda[i1]+55)*lambda[i2]*eta[i1]*eta[i1])/(4*s*s);
				phi[1]= 25.0/12.0*((500*lambda[i1]*lambda[i1]*lambda[i1]-450*lambda[i1]*lambda[i1]+110*lambda[i1]-6)*xi[i1]*xi[i2] + (750*lambda[i1]*lambda[i1]-450*lambda[i1]+55)*lambda[i2]*xi[i1]*xi[i1])/(4*s*s);
				phi[2]= -25.0/12.0*((250*lambda[i1]*lambda[i1]*lambda[i1]-225*lambda[i1]*lambda[i1]+55*lambda[i1]-3)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (750*lambda[i1]*lambda[i1]-450*lambda[i1]+55)*lambda[i2]*eta[i1]*xi[i1])/(4*s*s);
				break; 
			case 1:
				phi[0]= 25.0/12.0*((150*lambda[i1]-30)*(5*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*eta[i1] + (75*lambda[i1]*lambda[i1]-30*lambda[i1]+2)*(10*lambda[i2]-1)*2*eta[i1]*eta[i2] + (25*lambda[i1]*lambda[i1]*lambda[i1]-15*lambda[i1]*lambda[i1]+2*lambda[i1])*10*eta[i2]*eta[i2])/(4*s*s);
				phi[1]= 25.0/12.0*((150*lambda[i1]-30)*(5*lambda[i2]*lambda[i2]-lambda[i2])*xi[i1]*xi[i1] + (75*lambda[i1]*lambda[i1]-30*lambda[i1]+2)*(10*lambda[i2]-1)*2*xi[i1]*xi[i2] + (25*lambda[i1]*lambda[i1]*lambda[i1]-15*lambda[i1]*lambda[i1]+2*lambda[i1])*10*xi[i2]*xi[i2])/(4*s*s);
				phi[2]= -25.0/12.0*((150*lambda[i1]-30)*(5*lambda[i2]*lambda[i2]-lambda[i2])*eta[i1]*xi[i1] + (75*lambda[i1]*lambda[i1]-30*lambda[i1]+2)*(10*lambda[i2]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (25*lambda[i1]*lambda[i1]*lambda[i1]-15*lambda[i1]*lambda[i1]+2*lambda[i1])*10*eta[i2]*xi[i2])/(4*s*s);
				break;
			case 2:
				phi[0]= 25.0/12.0*((150*lambda[i2]-30)*(5*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*eta[i2] + (75*lambda[i2]*lambda[i2]-30*lambda[i2]+2)*(10*lambda[i1]-1)*2*eta[i1]*eta[i2] + (25*lambda[i2]*lambda[i2]*lambda[i2]-15*lambda[i2]*lambda[i2]+2*lambda[i2])*10*eta[i1]*eta[i1])/(4*s*s);
				phi[1]= 25.0/12.0*((150*lambda[i2]-30)*(5*lambda[i1]*lambda[i1]-lambda[i1])*xi[i2]*xi[i2] + (75*lambda[i2]*lambda[i2]-30*lambda[i2]+2)*(10*lambda[i1]-1)*2*xi[i1]*xi[i2] + (25*lambda[i2]*lambda[i2]*lambda[i2]-15*lambda[i2]*lambda[i2]+2*lambda[i2])*10*xi[i1]*xi[i1])/(4*s*s);
				phi[2]= -25.0/12.0*((150*lambda[i2]-30)*(5*lambda[i1]*lambda[i1]-lambda[i1])*eta[i2]*xi[i2] + (75*lambda[i2]*lambda[i2]-30*lambda[i2]+2)*(10*lambda[i1]-1)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (25*lambda[i2]*lambda[i2]*lambda[i2]-15*lambda[i2]*lambda[i2]+2*lambda[i2])*10*eta[i1]*xi[i1])/(4*s*s);
				break;
			case 3:
				phi[0]= 25.0/12.0*((500*lambda[i2]*lambda[i2]*lambda[i2]-450*lambda[i2]*lambda[i2]+110*lambda[i2]-6)*eta[i1]*eta[i2] + (750*lambda[i2]*lambda[i2]-450*lambda[i2]+55)*lambda[i1]*eta[i2]*eta[i2])/(4*s*s);
				phi[1]= 25.0/12.0*((500*lambda[i2]*lambda[i2]*lambda[i2]-450*lambda[i2]*lambda[i2]+110*lambda[i2]-6)*xi[i1]*xi[i2] + (750*lambda[i2]*lambda[i2]-450*lambda[i2]+55)*lambda[i1]*xi[i2]*xi[i2])/(4*s*s);
				phi[2]= -25.0/12.0*((250*lambda[i2]*lambda[i2]*lambda[i2]-225*lambda[i2]*lambda[i2]+55*lambda[i2]-3)*(eta[i1]*xi[i2]+eta[i2]*xi[i1]) + (750*lambda[i2]*lambda[i2]-450*lambda[i2]+55)*lambda[i1]*eta[i2]*xi[i2])/(4*s*s);
				break;
			default:
				phi[0] = 0;
				phi[1] = 0;
				phi[2] = 0;
			}
		}
		else if(index < 3*dop+3)
		{
			in = index-3*dop;
			phi[0]= 125.0/6.0*(2*(eta[0]*eta[1]*lambda[2]+eta[1]*eta[2]*lambda[0]+eta[2]*eta[0]*lambda[1])*(25*lambda[in]*lambda[in]-15*lambda[in]+2) + 2*(eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*eta[in] + lambda[0]*lambda[1]*lambda[2]*50*eta[in]*eta[in])/(4*s*s);
			phi[1]= 125.0/6.0*(2*(xi[0]*xi[1]*lambda[2]+xi[1]*xi[2]*lambda[0]+xi[2]*xi[0]*lambda[1])*(25*lambda[in]*lambda[in]-15*lambda[in]+2) + 2*(xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*xi[in] + lambda[0]*lambda[1]*lambda[2]*50*xi[in]*xi[in])/(4*s*s);
			phi[2]= -125.0/6.0*(((eta[0]*xi[1]+eta[1]*xi[0])*lambda[2]+(eta[1]*xi[2]+eta[2]*xi[1])*lambda[0]+(eta[2]*xi[0]+eta[0]*xi[2])*lambda[1])*(25*lambda[in]*lambda[in]-15*lambda[in]+2) + (xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*eta[in] + (eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*(50*lambda[in]-15)*xi[in] + lambda[0]*lambda[1]*lambda[2]*50*eta[in]*xi[in])/(4*s*s);
		}
		else
		{
			in = index-3*dop-3;
			i1 = (in+1)%3;
			i2 = (in+2)%3;
			phi[0]= 125.0/4.0*(2*(eta[0]*eta[1]*lambda[2]+eta[1]*eta[2]*lambda[0]+eta[2]*eta[0]*lambda[1])*(5*lambda[i1]-1)*(5*lambda[i2]-1) + 2*(eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*eta[i1]+(25*lambda[i1]-5)*eta[i2]) + lambda[0]*lambda[1]*lambda[2]*50*eta[i1]*eta[i2])/(4*s*s);
			phi[1]= 125.0/4.0*(2*(xi[0]*xi[1]*lambda[2]+xi[1]*xi[2]*lambda[0]+xi[2]*xi[0]*lambda[1])*(5*lambda[i1]-1)*(5*lambda[i2]-1) + 2*(xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*xi[i1]+(25*lambda[i1]-5)*xi[i2]) + lambda[0]*lambda[1]*lambda[2]*50*xi[i1]*xi[i2])/(4*s*s);
			phi[2]= -125.0/4.0*(((eta[0]*xi[1]+eta[1]*xi[0])*lambda[2]+(eta[1]*xi[2]+eta[2]*xi[1])*lambda[0]+(eta[2]*xi[0]+eta[0]*xi[2])*lambda[1])*(5*lambda[i1]-1)*(5*lambda[i2]-1) + (xi[0]*lambda[1]*lambda[2]+xi[1]*lambda[2]*lambda[0]+xi[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*eta[i1]+(25*lambda[i1]-5)*eta[i2]) + (eta[0]*lambda[1]*lambda[2]+eta[1]*lambda[2]*lambda[0]+eta[2]*lambda[0]*lambda[1])*((25*lambda[i2]-5)*xi[i1]+(25*lambda[i1]-5)*xi[i2]) + lambda[0]*lambda[1]*lambda[2]*25*(eta[i1]*xi[i2]+eta[i2]*xi[i1]))/(4*s*s);
		}
	} // dop=5
}

/**
* \fn void lagrangeSymTensor_basis(double *lambda, int index, int dop, double *phi)
* \brief basis function of Lagrange element
* \param *lambda pointer to the area coordiante
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi basis function
* \return void
*/
void lagrangeSymTensor_basis(double *lambda, int index, int dop, double phi[3])
{
	int in, ie, ii;
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom
	
	phi[0] = 0;
	phi[1] = 0;
	phi[2] = 0;
	if (index >= dofs*3 || index<0)
		return;

	lagrange_basis(lambda, index%dofs, dop, phi + index / dofs);
}

/**
* \fn void ncp1_basis(double *lambda, int index, double *phi)
* \brief basis function of nonconforming P1 element
* \param *lambda pointer to the area coordiante
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void ncp1_basis(double *lambda, int index, double *phi)
{
	if (index >= 3 || index<0)
	{
		*phi = 0;
		return;
	}

	*phi = 1 - 2 * lambda[index];
}

/**
* \fn void ncp1_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2])
* \brief the first order derivative of nonconforming P1 element basis function: (\partial_{x}phi, \partial_{y}phi)
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
* \return void
*/
void ncp1_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2])
{
	if (index >= 3 || index<0)
	{
		phi[0] = 0;
		phi[1] = 0;
		return;
	}

	phi[0] = -eta[index] / s;
	phi[1] = xi[index] / s;
}

/**
* \fn void mini_basis(double *lambda, int index, double *phi)
* \brief basis function of MINI element for Stokes equation
* \param *lambda pointer to the area coordiante
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void mini_basis(double *lambda, int index, double *phi)
{
	if (index >= 4 || index<0)
	{
		*phi = 0;
		return;
	}

	// the last dof is the function value at the centroid of the triangle
	if (index < 3)
		*phi = lambda[index] - 9.0*lambda[0] * lambda[1] * lambda[2];
	else
		*phi = 27.0*lambda[0] * lambda[1] * lambda[2];
	/********* the last dof is the mean value of the function on triangle
	if (index < 3)
		*phi = lambda[index] - 20.0*lambda[0] * lambda[1] * lambda[2];
	else
		*phi = 60.0*lambda[0] * lambda[1] * lambda[2];
		*/
}

/**
* \fn void mini_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2])
* \brief the first order derivative of MINI element basis function: (\partial_{x}phi, \partial_{y}phi)
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
* \return void
*/
void mini_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2])
{
	if (index >= 4 || index<0)
	{
		phi[0] = 0;
		phi[1] = 0;
		return;
	}

	// the last dof is the function value at the centroid of the triangle
	if (index < 3)
	{
		phi[0] = eta[index] / (2.0*s) - 9 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -xi[index] / (2.0*s) + 9 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}
	else
	{
		phi[0] = 27 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -27 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}

	/********* the last dof is the mean value of the function on triangle
	if (index < 3)
	{
		phi[0] = eta[index] / (2.0*s) - 20 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -xi[index] / (2.0*s) + 20 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}
	else
	{
		phi[0] = 60 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -60 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}*/
}

/**
* \fn void miniSymTensor_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3])
* \brief basis function of symmtric stress of MINI element for Stokes equation
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void miniSymTensor_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3])
{
	int dofs = 3 + 1 * 3 + 2; // degrees of freedom

	phi[0] = 0;
	phi[1] = 0;
	phi[2] = 0;
	if (index >= dofs || index<0)
		return;

	if (index < 3)
	{
		phi[0] = lambda[index];
		phi[1] = lambda[index];
	}
	else if (index < 6)
	{
		phi[index - 3] = 1;
	}
	else if (index == 6)
	{
		phi[0] = 27 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / sqrt(2.0*s);
		phi[2] = -13.5 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / sqrt(2.0*s);
	}
	else
	{
		phi[1] = -27 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / sqrt(2.0*s);
		phi[2] = 13.5 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / sqrt(2.0*s);
	}
}

/**
* \fn void miniSymTensorP1P0_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3])
* \brief basis function of symmtric stress of MINI element for Stokes equation
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void miniSymTensorP1P0_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3])
{
	int dofs = 3 * 2 + 1 + 2; // degrees of freedom

	phi[0] = 0;
	phi[1] = 0;
	phi[2] = 0;
	if (index >= dofs || index<0)
		return;

	if (index < 6)
		phi[index / 3] = lambda[index % 3];
	else if (index < 7)
		phi[2] = 1;
	else if (index == 7)
	{
		phi[0] = 27 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / sqrt(2.0*s);
		phi[2] = -13.5 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / sqrt(2.0*s);
	}
	else
	{
		phi[1] = -27 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / sqrt(2.0*s);
		phi[2] = 13.5 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / sqrt(2.0*s);
	}
}

/**
* \fn void cr_basis(double *lambda, int index, double *phi)
* \brief basis function of Crouzeix¨CRaviart element for Stokes equation
* \param *lambda pointer to the area coordiante
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void cr_basis(double *lambda, int index, double *phi)
{
	if (index >= 7 || index<0)
	{
		*phi = 0;
		return;
	}

	// the last dof is the function value at the centroid of the triangle
	if (index<3)
	{
		*phi = lambda[index] * (2 * lambda[index] - 1) + 3.0*lambda[0] * lambda[1] * lambda[2];
	}
	else if (index < 6)
	{
		*phi = 4.0*lambda[(index + 1) % 3] * lambda[(index + 2) % 3] - 12.0*lambda[0] * lambda[1] * lambda[2];
	}
	else
		*phi = 27.0*lambda[0] * lambda[1] * lambda[2];
	/********* the last dof is the mean value of the function on triangle
	if (index<3)
	{
		*phi = lambda[index] * (2 * lambda[index] - 1);
	}
	else if (index < 6)
	{
		*phi = 4.0*lambda[(index + 1) % 3] * lambda[(index + 2) % 3] - 20.0*lambda[0] * lambda[1] * lambda[2];
	}
	else
		*phi = 60.0*lambda[0] * lambda[1] * lambda[2];*/
}

/**
* \fn void cr_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2])
* \brief the first order derivative of Crouzeix¨CRaviart element basis function: (\partial_{x}phi, \partial_{y}phi)
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param phi[2] the first order derivative of Morley element basis function: (\partial_{x}phi, \partial_{y}phi)
* \return void
*/
void cr_basis1(double *lambda, double s, double eta[3], double xi[3], int index, double phi[2])
{
	int i1, i2;

	if (index >= 7 || index<0)
	{
		phi[0] = 0;
		phi[1] = 0;
		return;
	}

	// the last dof is the function value at the centroid of the triangle
	if (index < 3)
	{
		phi[0] = (4 * lambda[index] - 1)*eta[index] / (2.0*s) + 3 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -(4 * lambda[index] - 1)*xi[index] / (2.0*s) - 3 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}
	else if (index < 6)
	{
		i1 = (index + 1) % 3;
		i2 = (index + 2) % 3;
		phi[0] = 4 * (lambda[i1] * eta[i2] + lambda[i2] * eta[i1]) / (2.0*s) - 12 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -4 * (lambda[i1] * xi[i2] + lambda[i2] * xi[i1]) / (2.0*s) + 12 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}
	else
	{
		phi[0] = 27 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -27 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}

	/********* the last dof is the mean value of the function on triangle
	if (index < 3)
	{
		phi[0] = (4 * lambda[index] - 1)*eta[index] / (2.0*s);
		phi[1] = -(4 * lambda[index] - 1)*xi[index] / (2.0*s);
	}
	else if (index < 6)
	{
		i1 = (index + 1) % 3;
		i2 = (index + 2) % 3;
		phi[0] = 4 * (lambda[i1] * eta[i2] + lambda[i2] * eta[i1]) / (2.0*s) - 20 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -4 * (lambda[i1] * xi[i2] + lambda[i2] * xi[i1]) / (2.0*s) + 20 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}
	else
	{
		phi[0] = 60 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / (2.0*s);
		phi[1] = -60 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / (2.0*s);
	}*/
}

/**
* \fn void crSymTensor_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3])
* \brief basis function of symmtric stress of Crouzeix¨CRaviart element for Stokes equation
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param index the indicator of the basis function
* \param *phi basis function
* \return void
*/
void crSymTensor_basis(double *lambda, double s, double eta[3], double xi[3], int index, double phi[3])
{
	int dofs = 3 * 3 + 2; // degrees of freedom

	phi[0] = 0;
	phi[1] = 0;
	phi[2] = 0;
	if (index >= dofs || index<0)
		return;

	if (index < 9)
		phi[index / 3] = lambda[index % 3];
	else if (index == 9)
	{
		phi[0] = 27 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / sqrt(2.0*s);
		phi[2] = -13.5 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / sqrt(2.0*s);
	}
	else
	{
		phi[1] = -27 * (lambda[1] * lambda[2] * xi[0] + lambda[2] * lambda[0] * xi[1] + lambda[0] * lambda[1] * xi[2]) / sqrt(2.0*s);
		phi[2] = 13.5 * (lambda[1] * lambda[2] * eta[0] + lambda[2] * lambda[0] * eta[1] + lambda[0] * lambda[1] * eta[2]) / sqrt(2.0*s);
	}
}

/** 
 * \fn void rt_basis(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[2])
 * \brief basis function of Raviart-Thomas element: (phi1, phi2)
 * \param x the x-axis coordiante
 * \param y the y-axis coordiante
 * \param (*T)[2] point the coordiantes of all vertices of current element
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param orient[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[2] basis function of Raviart-Thomas element: (phi1, phi2)
 * \return void
 */
void rt_basis(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[2])
{
	double c1, c2, d1, d2, d3, a1, a2, a3, b1, b2, b3, aa1, aa2, aa3, bb1, bb2, bb3;
	double x1 = T[0][0];
	double x2 = T[1][0];
	double x3 = T[2][0];
	double y1 = T[0][1];
	double y2 = T[1][1];
	double y3 = T[2][1];
	if(dop==1)
	{
		if(index==0)
		{
			c1 = -elen[0]*eta[0]/(s*s);
			c2 = elen[0]*xi[0]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2-elen[0]*xi[2]/(2*s);
			a3 = -d3*x3+elen[0]*xi[1]/(2*s);
			b1 = -d1*y1;
			b2 = -d2*y2-elen[0]*eta[2]/(2*s);
			b3 = -d3*y3+elen[0]*eta[1]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[0];
			c2 *= orient[0];
			aa1 *= orient[0];
			aa2 *= orient[0];
			aa3 *= orient[0];
			bb1 *= orient[0];
			bb2 *= orient[0];
			bb3 *= orient[0];
		}
		else if(index==1)
		{
			c1 = -2*elen[0]*(eta[1]-eta[2])/(s*s);
			c2 = 2*elen[0]*(xi[1]-xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2+3*elen[0]*xi[2]/s;
			a3 = -d3*x3+3*elen[0]*xi[1]/s;
			b1 = -d1*y1;
			b2 = -d2*y2+3*elen[0]*eta[2]/s;
			b3 = -d3*y3+3*elen[0]*eta[1]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==2)
		{
			c1 = -elen[1]*eta[1]/(s*s);
			c2 = elen[1]*xi[1]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+elen[1]*xi[2]/(2*s);
			a2 = -d2*x2;
			a3 = -d3*x3-elen[1]*xi[0]/(2*s);
			b1 = -d1*y1+elen[1]*eta[2]/(2*s);
			b2 = -d2*y2;
			b3 = -d3*y3-elen[1]*eta[0]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[1];
			c2 *= orient[1];
			aa1 *= orient[1];
			aa2 *= orient[1];
			aa3 *= orient[1];
			bb1 *= orient[1];
			bb2 *= orient[1];
			bb3 *= orient[1];
		}
		else if(index==3)
		{
			c1 = -2*elen[1]*(eta[2]-eta[0])/(s*s);
			c2 = 2*elen[1]*(xi[2]-xi[0])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[1]*xi[2]/s;
			a2 = -d2*x2;
			a3 = -d3*x3+3*elen[1]*xi[0]/s;
			b1 = -d1*y1+3*elen[1]*eta[2]/s;
			b2 = -d2*y2;
			b3 = -d3*y3+3*elen[1]*eta[0]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==4)
		{
			c1 = -elen[2]*eta[2]/(s*s);
			c2 = elen[2]*xi[2]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1-elen[2]*xi[1]/(2*s);
			a2 = -d2*x2+elen[2]*xi[0]/(2*s);
			a3 = -d3*x3;
			b1 = -d1*y1-elen[2]*eta[1]/(2*s);
			b2 = -d2*y2+elen[2]*eta[0]/(2*s);
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[2];
			c2 *= orient[2];
			aa1 *= orient[2];
			aa2 *= orient[2];
			aa3 *= orient[2];
			bb1 *= orient[2];
			bb2 *= orient[2];
			bb3 *= orient[2];
		}
		else if(index==5)
		{
			c1 = -2*elen[2]*(eta[0]-eta[1])/(s*s);
			c2 = 2*elen[2]*(xi[0]-xi[1])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[2]*xi[1]/s;
			a2 = -d2*x2+3*elen[2]*xi[0]/s;
			a3 = -d3*x3;
			b1 = -d1*y1+3*elen[2]*eta[1]/s;
			b2 = -d2*y2+3*elen[2]*eta[0]/s;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==6)
		{
			c1 = -(eta[0]*eta[0]+eta[1]*eta[1]+eta[2]*eta[2])/(s*s);
			c2 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==7)
		{
			c1 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			c2 = -(xi[0]*xi[0]+xi[1]*xi[1]+xi[2]*xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
			aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
			bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else
		{
			c1 = 0;
			c2 = 0;
			d1 = 0;
			d2 = 0;
			d3 = 0;
			a1 = 0;
			a2 = 0;
			a3 = 0;
			b1 = 0;
			b2 = 0;
			b3 = 0;
			aa1 = 0;
			aa2 = 0;
			aa3 = 0;
			bb1 = 0;
			bb2 = 0;
			bb3 = 0;
		}
		phi[0] = aa1*x+aa2*y+aa3+x*(c1*x+c2*y);
		phi[1] = bb1*x+bb2*y+bb3+y*(c1*x+c2*y);
	} // dop=1
	else
	{
		phi[0]=0;
		phi[1]=0;
	}
}

/** 
 * \fn void rt_basis1(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[4])
 * \brief the first order derivative of Raviart-Thomas element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
 * \param x the x-axis coordiante
 * \param y the y-axis coordiante
 * \param (*T)[2] point the coordiantes of all vertices of current element
 * \param s the area of the triangule
 * \param elen[3] length of three edge
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param orient[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param phi[4] the first order derivative of Raviart-Thomas element basis function: (\partial_{x}phi1, \partial_{y}phi1, \partial_{x}phi2, \partial_{y}phi2)
 * \return void
 */
void rt_basis1(double x, double y, double (*T)[2], double s, double elen[3], double eta[3], double xi[3], double orient[3], int index, int dop, double phi[4])
{
	double c1, c2, d1, d2, d3, a1, a2, a3, b1, b2, b3, aa1, aa2, aa3, bb1, bb2, bb3;
	double x1 = T[0][0];
	double x2 = T[1][0];
	double x3 = T[2][0];
	double y1 = T[0][1];
	double y2 = T[1][1];
	double y3 = T[2][1];
	if(dop==1)
	{
		if(index==0)
		{
			c1 = -elen[0]*eta[0]/(s*s);
			c2 = elen[0]*xi[0]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2-elen[0]*xi[2]/(2*s);
			a3 = -d3*x3+elen[0]*xi[1]/(2*s);
			b1 = -d1*y1;
			b2 = -d2*y2-elen[0]*eta[2]/(2*s);
			b3 = -d3*y3+elen[0]*eta[1]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[0];
			c2 *= orient[0];
			aa1 *= orient[0];
			aa2 *= orient[0];
	//		aa3 *= orient[0];
			bb1 *= orient[0];
			bb2 *= orient[0];
	//		bb3 *= orient[0];
		}
		else if(index==1)
		{
			c1 = -2*elen[0]*(eta[1]-eta[2])/(s*s);
			c2 = 2*elen[0]*(xi[1]-xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2+3*elen[0]*xi[2]/s;
			a3 = -d3*x3+3*elen[0]*xi[1]/s;
			b1 = -d1*y1;
			b2 = -d2*y2+3*elen[0]*eta[2]/s;
			b3 = -d3*y3+3*elen[0]*eta[1]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==2)
		{
			c1 = -elen[1]*eta[1]/(s*s);
			c2 = elen[1]*xi[1]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+elen[1]*xi[2]/(2*s);
			a2 = -d2*x2;
			a3 = -d3*x3-elen[1]*xi[0]/(2*s);
			b1 = -d1*y1+elen[1]*eta[2]/(2*s);
			b2 = -d2*y2;
			b3 = -d3*y3-elen[1]*eta[0]/(2*s);
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[1];
			c2 *= orient[1];
			aa1 *= orient[1];
			aa2 *= orient[1];
	//		aa3 *= orient[1];
			bb1 *= orient[1];
			bb2 *= orient[1];
	//		bb3 *= orient[1];
		}
		else if(index==3)
		{
			c1 = -2*elen[1]*(eta[2]-eta[0])/(s*s);
			c2 = 2*elen[1]*(xi[2]-xi[0])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[1]*xi[2]/s;
			a2 = -d2*x2;
			a3 = -d3*x3+3*elen[1]*xi[0]/s;
			b1 = -d1*y1+3*elen[1]*eta[2]/s;
			b2 = -d2*y2;
			b3 = -d3*y3+3*elen[1]*eta[0]/s;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==4)
		{
			c1 = -elen[2]*eta[2]/(s*s);
			c2 = elen[2]*xi[2]/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1-elen[2]*xi[1]/(2*s);
			a2 = -d2*x2+elen[2]*xi[0]/(2*s);
			a3 = -d3*x3;
			b1 = -d1*y1-elen[2]*eta[1]/(2*s);
			b2 = -d2*y2+elen[2]*eta[0]/(2*s);
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
			c1 *= orient[2];
			c2 *= orient[2];
			aa1 *= orient[2];
			aa2 *= orient[2];
	//		aa3 *= orient[2];
			bb1 *= orient[2];
			bb2 *= orient[2];
	//		bb3 *= orient[2];
		}
		else if(index==5)
		{
			c1 = -2*elen[2]*(eta[0]-eta[1])/(s*s);
			c2 = 2*elen[2]*(xi[0]-xi[1])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1+3*elen[2]*xi[1]/s;
			a2 = -d2*x2+3*elen[2]*xi[0]/s;
			a3 = -d3*x3;
			b1 = -d1*y1+3*elen[2]*eta[1]/s;
			b2 = -d2*y2+3*elen[2]*eta[0]/s;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==6)
		{
			c1 = -(eta[0]*eta[0]+eta[1]*eta[1]+eta[2]*eta[2])/(s*s);
			c2 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else if(index==7)
		{
			c1 = (xi[0]*eta[0]+xi[1]*eta[1]+xi[2]*eta[2])/(s*s);
			c2 = -(xi[0]*xi[0]+xi[1]*xi[1]+xi[2]*xi[2])/(s*s);
			d1 = c1*x1+c2*y1;
			d2 = c1*x2+c2*y2;
			d3 = c1*x3+c2*y3;
			a1 = -d1*x1;
			a2 = -d2*x2;
			a3 = -d3*x3;
			b1 = -d1*y1;
			b2 = -d2*y2;
			b3 = -d3*y3;
			aa1 = (a1*eta[0]+a2*eta[1]+a3*eta[2])/(2*s);
			aa2 = (-a1*xi[0]-a2*xi[1]-a3*xi[2])/(2*s);
	//		aa3 = (a1*(x2*y3-y2*x3)+a2*(x3*y1-y3*x1)+a3*(x1*y2-y1*x2))/(2*s);
			bb1 = (b1*eta[0]+b2*eta[1]+b3*eta[2])/(2*s);
			bb2 = (-b1*xi[0]-b2*xi[1]-b3*xi[2])/(2*s);
	//		bb3 = (b1*(x2*y3-y2*x3)+b2*(x3*y1-y3*x1)+b3*(x1*y2-y1*x2))/(2*s);
		}
		else
		{
			c1 = 0;
			c2 = 0;
			d1 = 0;
			d2 = 0;
			d3 = 0;
			a1 = 0;
			a2 = 0;
			a3 = 0;
			b1 = 0;
			b2 = 0;
			b3 = 0;
			aa1 = 0;
			aa2 = 0;
	//		aa3 = 0;
			bb1 = 0;
			bb2 = 0;
	//		bb3 = 0;
		}
		phi[0] = aa1+2*c1*x+c2*y;
		phi[1] = aa2+c2*x;
		phi[2] = bb1+c1*y;
		phi[3] = bb2+c1*x+2*c2*y;
	} // dop=1
	else
	{
		phi[0]=0;
		phi[1]=0;
		phi[2]=0;
		phi[3]=0;
	}
}

/** 
 * \fn void arnoldwinther_basis(double *lambda, double *x, double *y, ddenmat3 *basisCoeffs, int element, int index, double *phi)
 * \brief basis function of Arnold-Winther element
 * \param *lambda pointer to the area coordiante
 * \param *x pointer to the horizontal ordinates of three vertices
 * \param *y pointer to the longitudinal ordinates of three vertices
 * \param *basisCoeffs pointer to coefficients of basis functions
 * \param element the current element of triangulation
 * \param index the indicator of the basis function
 * \param *phi basis function
 * \return void
 */
void arnoldwinther_basis(double *lambda, double *x, double *y, ddenmat3 *basisCoeffs, int element, int index, double *phi)
{
	double xx=lambda[0]*x[0]+lambda[1]*x[1]+lambda[2]*x[2];
	double yy=lambda[0]*y[0]+lambda[1]*y[1]+lambda[2]*y[2];
	double *coeffs=basisCoeffs->val[element][index];
	
	phi[0] = coeffs[0]*lambda[0] + coeffs[3]*lambda[1] + coeffs[6]*lambda[2] + coeffs[9]*lambda[1]*lambda[2] + coeffs[12]*lambda[2]*lambda[0] + coeffs[15]*lambda[0]*lambda[1];
	phi[1] = coeffs[1]*lambda[0] + coeffs[4]*lambda[1] + coeffs[7]*lambda[2] + coeffs[10]*lambda[1]*lambda[2] + coeffs[13]*lambda[2]*lambda[0] + coeffs[16]*lambda[0]*lambda[1];
	phi[2] = coeffs[2]*lambda[0] + coeffs[5]*lambda[1] + coeffs[8]*lambda[2] + coeffs[11]*lambda[1]*lambda[2] + coeffs[14]*lambda[2]*lambda[0] + coeffs[17]*lambda[0]*lambda[1];

	phi[0] += coeffs[18]*pow(yy,3) - coeffs[21]*3*xx*yy*yy - coeffs[22]*pow(xx,3)/3.0 - coeffs[23]*xx*xx*yy;
	phi[1] += coeffs[19]*pow(xx,3) - coeffs[20]*3*xx*xx*yy - coeffs[22]*xx*yy*yy - coeffs[23]*pow(yy,3)/3.0;
	phi[2] += coeffs[20]*pow(xx,3) + coeffs[21]*pow(yy,3) + coeffs[22]*xx*xx*yy + coeffs[23]*xx*yy*yy;
}

/** 
 * \fn void arnoldwinther_basisDIV(double *lambda, ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
 * \brief divergence of basis function of Arnold-Winther element
 * \param *lambda pointer to the area coordiante
 * \param *basisCoeffs pointer to coefficients of basis functions
 * \param element the current element of triangulation
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param *phi divergence of basis function
 * \return void
 */
void arnoldwinther_basisDIV(double *lambda, ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
{
	double *coeffs=basisCoeffs->val[element][index];
	double lambdaGrad[3][2];
	int i;
	for(i=0;i<3;i++)
	{
		lambdaGrad[i][0]=eta[i]/(2.0*s);
		lambdaGrad[i][1]=-xi[i]/(2.0*s);
	}

	phi[0] = coeffs[0]*lambdaGrad[0][0] + coeffs[3]*lambdaGrad[1][0] + coeffs[6]*lambdaGrad[2][0];
	phi[0] += coeffs[9]*(lambda[1]*lambdaGrad[2][0]+lambda[2]*lambdaGrad[1][0]) + coeffs[12]*(lambda[2]*lambdaGrad[0][0]+lambda[0]*lambdaGrad[2][0]) + coeffs[15]*(lambda[0]*lambdaGrad[1][0]+lambda[1]*lambdaGrad[0][0]);
	phi[0] += coeffs[2]*lambdaGrad[0][1] + coeffs[5]*lambdaGrad[1][1] + coeffs[8]*lambdaGrad[2][1];
	phi[0] += coeffs[11]*(lambda[1]*lambdaGrad[2][1]+lambda[2]*lambdaGrad[1][1]) + coeffs[14]*(lambda[2]*lambdaGrad[0][1]+lambda[0]*lambdaGrad[2][1]) + coeffs[17]*(lambda[0]*lambdaGrad[1][1]+lambda[1]*lambdaGrad[0][1]);

	phi[1] = coeffs[1]*lambdaGrad[0][1] + coeffs[4]*lambdaGrad[1][1] + coeffs[7]*lambdaGrad[2][1];
	phi[1] += coeffs[10]*(lambda[1]*lambdaGrad[2][1]+lambda[2]*lambdaGrad[1][1]) + coeffs[13]*(lambda[2]*lambdaGrad[0][1]+lambda[0]*lambdaGrad[2][1]) + coeffs[16]*(lambda[0]*lambdaGrad[1][1]+lambda[1]*lambdaGrad[0][1]);
	phi[1] += coeffs[2]*lambdaGrad[0][0] + coeffs[5]*lambdaGrad[1][0] + coeffs[8]*lambdaGrad[2][0];
	phi[1] += coeffs[11]*(lambda[1]*lambdaGrad[2][0]+lambda[2]*lambdaGrad[1][0]) + coeffs[14]*(lambda[2]*lambdaGrad[0][0]+lambda[0]*lambdaGrad[2][0]) + coeffs[17]*(lambda[0]*lambdaGrad[1][0]+lambda[1]*lambdaGrad[0][0]);
}

/** 
 * \fn void arnoldwinther_basisDIV2(ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
 * \brief divergence of divergence of basis function of Arnold-Winther element
 * \param *basisCoeffs pointer to coefficients of basis functions
 * \param element the current element of triangulation
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param index the indicator of the basis function
 * \param *phi divergence of divergence of basis function
 * \return void
 */
void arnoldwinther_basisDIV2(ddenmat3 *basisCoeffs, int element, double s, double eta[3], double xi[3], int index, double *phi)
{
	double *coeffs=basisCoeffs->val[element][index];
	double lambdaGrad[3][2];
	int i;
	for(i=0;i<3;i++)
	{
		lambdaGrad[i][0]=eta[i]/(2.0*s);
		lambdaGrad[i][1]=-xi[i]/(2.0*s);
	}

	*phi = coeffs[9]*2*lambdaGrad[1][0]*lambdaGrad[2][0] + coeffs[12]*2*lambdaGrad[2][0]*lambdaGrad[0][0] + coeffs[15]*2*lambdaGrad[0][0]*lambdaGrad[1][0];
	*phi += coeffs[10]*2*lambdaGrad[1][1]*lambdaGrad[2][1] + coeffs[13]*2*lambdaGrad[2][1]*lambdaGrad[0][1] + coeffs[16]*2*lambdaGrad[0][1]*lambdaGrad[1][1];
	*phi += coeffs[11]*2*(lambdaGrad[1][0]*lambdaGrad[2][1]+lambdaGrad[1][1]*lambdaGrad[2][0]);
	*phi += coeffs[14]*2*(lambdaGrad[2][0]*lambdaGrad[0][1]+lambdaGrad[2][1]*lambdaGrad[0][0]);
	*phi += coeffs[17]*2*(lambdaGrad[0][0]*lambdaGrad[1][1]+lambdaGrad[0][1]*lambdaGrad[1][0]);
}

/** 
 * \fn void huzhang_basis(double *lambda, double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[3])
 * \brief basis function of Hu-Zhang element
 * \param *lambda pointer to the area coordiante
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param **tensorBasis[] point to local tensor basis
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi basis function
 * \return void
 */
void huzhang_basis(double *lambda, double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[3])
{
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom

	phi[0]=0;
	phi[1]=0;
	phi[2]=0;
	if(index>= dofs*3 || index<0)
		return;

	double val;
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for(i=0;i<3;i++)
	{
		nn[i][0]=nv[i][0]*nv[i][0]; nn[i][1]=nv[i][1]*nv[i][1]; nn[i][2]=nv[i][0]*nv[i][1];
		nt[i][0]=nv[i][0]*tv[i][0]*sqrt(2); nt[i][1]=nv[i][1]*tv[i][1]*sqrt(2); nt[i][2]=(nv[i][0]*tv[i][1]+nv[i][1]*tv[i][0])/sqrt(2);
		tt[i][0]=tv[i][0]*tv[i][0]; tt[i][1]=tv[i][1]*tv[i][1]; tt[i][2]=tv[i][0]*tv[i][1];
	}

	if(dop==0)
	{
		phi[index] = 1;
	}
	
	else // dop>=1
	{
		if(index<9)
		{
			lagrange_basis(lambda, index%3, dop, &val);
			phi[0] = val*tensorBasis[index % 3][index / 3][0];
			phi[1] = val*tensorBasis[index % 3][index / 3][1];
			phi[2] = val*tensorBasis[index % 3][index / 3][2];
		}
		else if(index<9+(dop-1)*3)
		{
			i=index-9;
			lagrange_basis(lambda, 3+i, dop, &val);
			phi[0]=val*nn[i/(dop-1)][0];
			phi[1]=val*nn[i/(dop-1)][1];
			phi[2]=val*nn[i/(dop-1)][2];
		}
		else if(index<9+(dop-1)*6)
		{
			i=index-9-(dop-1)*3;
			lagrange_basis(lambda, 3+i, dop, &val);
			phi[0]=val*nt[i/(dop-1)][0];
			phi[1]=val*nt[i/(dop-1)][1];
			phi[2]=val*nt[i/(dop-1)][2];
		}
		else if(index<9+(dop-1)*9)
		{
			i=index-9-(dop-1)*6;
			lagrange_basis(lambda, 3+i, dop, &val);
			phi[0]=val*tt[i/(dop-1)][0];
			phi[1]=val*tt[i/(dop-1)][1];
			phi[2]=val*tt[i/(dop-1)][2];
		}
		else
		{
			i=index-dop*9;
			lagrange_basis(lambda, dop*3+i%(dofs-dop*3), dop, &val);
			phi[i/(dofs-dop*3)]=val;
		}
	} 
	
}

/**
* \fn void huzhang_basis1(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double(*phi)[2])
* \brief the first order derivative of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param **tensorBasis[] point to local tensor basis
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param double(*phi)[2] the first order derivative of basis function of Hu-Zhang element
* \return void
*/
void huzhang_basis1(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double(*phi)[2])
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	phi[0][0] = 0; phi[0][1] = 0;
	phi[1][0] = 0; phi[1][1] = 0;
	phi[2][0] = 0; phi[2][1] = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop == 0)
	{
		return;
	}

	else // dop>=1
	{
		if (index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0][0] = val[0] * tensorBasis[index % 3][index / 3][0];
			phi[0][1] = val[1] * tensorBasis[index % 3][index / 3][0];
			phi[1][0] = val[0] * tensorBasis[index % 3][index / 3][1];
			phi[1][1] = val[1] * tensorBasis[index % 3][index / 3][1];
			phi[2][0] = val[0] * tensorBasis[index % 3][index / 3][2];
			phi[2][1] = val[1] * tensorBasis[index % 3][index / 3][2];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0][0] = val[0] * nn[i / (dop - 1)][0];
			phi[0][1] = val[1] * nn[i / (dop - 1)][0];
			phi[1][0] = val[0] * nn[i / (dop - 1)][1];
			phi[1][1] = val[1] * nn[i / (dop - 1)][1];
			phi[2][0] = val[0] * nn[i / (dop - 1)][2];
			phi[2][1] = val[1] * nn[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0][0] = val[0] * nt[i / (dop - 1)][0];
			phi[0][1] = val[1] * nt[i / (dop - 1)][0];
			phi[1][0] = val[0] * nt[i / (dop - 1)][1];
			phi[1][1] = val[1] * nt[i / (dop - 1)][1];
			phi[2][0] = val[0] * nt[i / (dop - 1)][2];
			phi[2][1] = val[1] * nt[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0][0] = val[0] * tt[i / (dop - 1)][0];
			phi[0][1] = val[1] * tt[i / (dop - 1)][0];
			phi[1][0] = val[0] * tt[i / (dop - 1)][1];
			phi[1][1] = val[1] * tt[i / (dop - 1)][1];
			phi[2][0] = val[0] * tt[i / (dop - 1)][2];
			phi[2][1] = val[1] * tt[i / (dop - 1)][2];
		}
		else
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i % (dofs - dop * 3), dop, val);
			phi[i / (dofs - dop * 3)][0] = val[0];
			phi[i / (dofs - dop * 3)][1] = val[1];
		}
	}

}

/** 
 * \fn void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2])
 * \brief divergence of basis function of Hu-Zhang element
 * \param *lambda pointer to the area coordiante
 * \param s the area of the triangule
 * \param eta[3] some auxiliary parameter
 * \param xi[3] some auxiliary parameter
 * \param **nv the unit normal vectors of the three edges
 * \param **tv the unit tangential vectors of the three edges
 * \param **tensorBasis[] point to local tensor basis
 * \param index the indicator of the basis function
 * \param dop degree of polynomial
 * \param *phi divergence of basis function
 * \return void
 */
void huzhang_basisDIV(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2])
{
	int dofs = (dop+1)*(dop+2)/2; // degrees of freedom

	phi[0]=0;
	phi[1]=0;
	if(index>= dofs*3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for(i=0;i<3;i++)
	{
		nn[i][0]=nv[i][0]*nv[i][0]; nn[i][1]=nv[i][1]*nv[i][1]; nn[i][2]=nv[i][0]*nv[i][1];
		nt[i][0]=nv[i][0]*tv[i][0]*sqrt(2); nt[i][1]=nv[i][1]*tv[i][1]*sqrt(2); nt[i][2]=(nv[i][0]*tv[i][1]+nv[i][1]*tv[i][0])/sqrt(2);
		tt[i][0]=tv[i][0]*tv[i][0]; tt[i][1]=tv[i][1]*tv[i][1]; tt[i][2]=tv[i][0]*tv[i][1];
	}

	if(dop==0)
	{
		return;
	}
	
	else // dop>=1
	{
		if(index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0] = val[0] * tensorBasis[index % 3][index / 3][0] + val[1] * tensorBasis[index % 3][index / 3][2];
			phi[1] = val[0] * tensorBasis[index % 3][index / 3][2] + val[1] * tensorBasis[index % 3][index / 3][1];
		}
		else if(index<9+(dop-1)*3)
		{
			i=index-9;
			lagrange_basis1(lambda, s, eta, xi, 3+i, dop, val);
			phi[0]=val[0]*nn[i/(dop-1)][0]+val[1]*nn[i/(dop-1)][2];
			phi[1]=val[0]*nn[i/(dop-1)][2]+val[1]*nn[i/(dop-1)][1];
		}
		else if(index<9+(dop-1)*6)
		{
			i=index-9-(dop-1)*3;
			lagrange_basis1(lambda, s, eta, xi, 3+i, dop, val);
			phi[0]=val[0]*nt[i/(dop-1)][0]+val[1]*nt[i/(dop-1)][2];
			phi[1]=val[0]*nt[i/(dop-1)][2]+val[1]*nt[i/(dop-1)][1];
		}
		else if(index<9+(dop-1)*9)
		{
			i=index-9-(dop-1)*6;
			lagrange_basis1(lambda, s, eta, xi, 3+i, dop, val);
			phi[0]=val[0]*tt[i/(dop-1)][0]+val[1]*tt[i/(dop-1)][2];
			phi[1]=val[0]*tt[i/(dop-1)][2]+val[1]*tt[i/(dop-1)][1];
		}
		else if(index<dop*9+(dofs-dop*3))
		{
			i=index-dop*9;
			lagrange_basis1(lambda, s, eta, xi, dop*3+i, dop, val);
			phi[0] = val[0];
			phi[1] = 0;
		}
		else if(index<dop*9+(dofs-dop*3)*2)
		{
			i=index-dop*9-(dofs-dop*3);
			lagrange_basis1(lambda, s, eta, xi, dop*3+i, dop, val);
			phi[0] = 0;
			phi[1] = val[1];
		}
		else
		{
			i=index-dop*9-(dofs-dop*3)*2;
			lagrange_basis1(lambda, s, eta, xi, dop*3+i, dop, val);
			phi[0] = val[1];
			phi[1] = val[0];
		}
	} 
}

/**
* \fn void huzhang_basisROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2])
* \brief rotation of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param **tensorBasis[] point to local tensor basis
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2])
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	phi[0] = 0;
	phi[1] = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop == 0)
	{
		return;
	}

	else // dop>=1
	{
		if (index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0] = -val[1] * tensorBasis[index % 3][index / 3][0] + val[0] * tensorBasis[index % 3][index / 3][2];
			phi[1] = -val[1] * tensorBasis[index % 3][index / 3][2] + val[0] * tensorBasis[index % 3][index / 3][1];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * nn[i / (dop - 1)][0] + val[0] * nn[i / (dop - 1)][2];
			phi[1] = -val[1] * nn[i / (dop - 1)][2] + val[0] * nn[i / (dop - 1)][1];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * nt[i / (dop - 1)][0] + val[0] * nt[i / (dop - 1)][2];
			phi[1] = -val[1] * nt[i / (dop - 1)][2] + val[0] * nt[i / (dop - 1)][1];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * tt[i / (dop - 1)][0] + val[0] * tt[i / (dop - 1)][2];
			phi[1] = -val[1] * tt[i / (dop - 1)][2] + val[0] * tt[i / (dop - 1)][1];
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = 0;
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = 0;
			phi[1] = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = val[0];
			phi[1] = -val[1];
		}
	}
}

/**
* \fn void huzhang_basisCurlTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2])
* \brief curl of the trace of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param **tensorBasis[] point to local tensor basis
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisCurlTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double phi[2])
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	phi[0] = 0;
	phi[1] = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[2];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop == 0)
	{
		return;
	}

	else // dop>=1
	{
		if (index<9)
		{
			lagrange_basis1(lambda, s, eta, xi, index % 3, dop, val);
			phi[0] = -val[1] * (tensorBasis[index % 3][index / 3][0] + tensorBasis[index % 3][index / 3][1]);
			phi[1] = val[0] * (tensorBasis[index % 3][index / 3][0] + tensorBasis[index % 3][index / 3][1]);
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
			phi[1] = val[0] * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
			phi[1] = val[0] * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis1(lambda, s, eta, xi, 3 + i, dop, val);
			phi[0] = -val[1] * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
			phi[1] = val[0] * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = -val[1];
			phi[1] = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis1(lambda, s, eta, xi, dop * 3 + i, dop, val);
			phi[0] = 0;
			phi[1] = 0;
		}
	}
}

/**
* \fn void huzhang_basisROTROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double *phi)
* \brief rotrot of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param **tensorBasis[] point to local tensor basis
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisROTROT(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double *phi)
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	*phi = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[3];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop < 2)
	{
		return;
	}

	else // dop>=2
	{
		if (index<9)
		{
			lagrange_basis2(lambda, s, eta, xi, index % 3, dop, val);
			*phi = val[1] * tensorBasis[index % 3][index / 3][0] + val[0] * tensorBasis[index % 3][index / 3][1] - 2 * val[2] * tensorBasis[index % 3][index / 3][2];
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = val[1] * nn[i / (dop - 1)][0] + val[0] * nn[i / (dop - 1)][1] - 2 * val[2] * nn[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = val[1] * nt[i / (dop - 1)][0] + val[0] * nt[i / (dop - 1)][1] - 2 * val[2] * nt[i / (dop - 1)][2];
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = val[1] * tt[i / (dop - 1)][0] + val[0] * tt[i / (dop - 1)][1] - 2 * val[2] * tt[i / (dop - 1)][2];
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[1];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[0];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = -2 * val[2];
		}
	}
}

/**
* \fn void huzhang_basisLaplaceTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double *phi)
* \brief Laplace of the trace of basis function of Hu-Zhang element
* \param *lambda pointer to the area coordiante
* \param s the area of the triangule
* \param eta[3] some auxiliary parameter
* \param xi[3] some auxiliary parameter
* \param **nv the unit normal vectors of the three edges
* \param **tv the unit tangential vectors of the three edges
* \param **tensorBasis[] point to local tensor basis
* \param index the indicator of the basis function
* \param dop degree of polynomial
* \param *phi divergence of basis function
* \return void
*/
void huzhang_basisLaplaceTrace(double *lambda, double s, double eta[3], double xi[3], double **nv, double **tv, double **tensorBasis[], int index, int dop, double *phi)
{
	int dofs = (dop + 1)*(dop + 2) / 2; // degrees of freedom

	*phi = 0;
	if (index >= dofs * 3 || index<0)
		return;

	double val[3];
	double nn[3][3], nt[3][3], tt[3][3];
	int i;
	for (i = 0; i<3; i++)
	{
		nn[i][0] = nv[i][0] * nv[i][0]; nn[i][1] = nv[i][1] * nv[i][1]; nn[i][2] = nv[i][0] * nv[i][1];
		nt[i][0] = nv[i][0] * tv[i][0] * sqrt(2); nt[i][1] = nv[i][1] * tv[i][1] * sqrt(2); nt[i][2] = (nv[i][0] * tv[i][1] + nv[i][1] * tv[i][0]) / sqrt(2);
		tt[i][0] = tv[i][0] * tv[i][0]; tt[i][1] = tv[i][1] * tv[i][1]; tt[i][2] = tv[i][0] * tv[i][1];
	}

	if (dop < 2)
	{
		return;
	}

	else // dop>=2
	{
		if (index<9)
		{
			lagrange_basis2(lambda, s, eta, xi, index % 3, dop, val);
			*phi = (val[0] + val[1]) * (tensorBasis[index % 3][index / 3][0] + tensorBasis[index % 3][index / 3][1]);
		}
		else if (index<9 + (dop - 1) * 3)
		{
			i = index - 9;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (nn[i / (dop - 1)][0] + nn[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 6)
		{
			i = index - 9 - (dop - 1) * 3;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (nt[i / (dop - 1)][0] + nt[i / (dop - 1)][1]);
		}
		else if (index<9 + (dop - 1) * 9)
		{
			i = index - 9 - (dop - 1) * 6;
			lagrange_basis2(lambda, s, eta, xi, 3 + i, dop, val);
			*phi = (val[0] + val[1]) * (tt[i / (dop - 1)][0] + tt[i / (dop - 1)][1]);
		}
		else if (index<dop * 9 + (dofs - dop * 3))
		{
			i = index - dop * 9;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[0] + val[1];
		}
		else if (index<dop * 9 + (dofs - dop * 3) * 2)
		{
			i = index - dop * 9 - (dofs - dop * 3);
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = val[0] + val[1];
		}
		else
		{
			i = index - dop * 9 - (dofs - dop * 3) * 2;
			lagrange_basis2(lambda, s, eta, xi, dop * 3 + i, dop, val);
			*phi = 0;
		}
	}
}

/**
 * \fn double area(double x1,double x2,double x3,double y1,double y2,double y3)
 * \brief get area for triangle p1(x1,y1),p2(x2,y2),p3(x3,y3)
 * \param x1 the x-axis value of the point p1
 * \param x2 the x-axis value of the point p2
 * \param x3 the x-axis value of the point p3
 * \param y1 the y-axis value of the point p1
 * \param y2 the y-axis value of the point p2
 * \param y3 the y-axis value of the point p3
 * \return area of the trianle
 */
double area(double x1,double x2,double x3,double y1,double y2,double y3)
{
	return ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))/2;
}

/** 
 * \fn void localb(double (*nodes)[2],double *b)
 * \brief get local right-hand side b from triangle nodes
 * \param (*nodes)[2] the vertice of the triangule
 * \param *b local right-hand side
 * \return void
 */
void localb(double (*nodes)[2],double *b)
{
	int i;
	double x,y,a;
	double s=area(nodes[0][0],nodes[1][0],nodes[2][0],nodes[0][1],nodes[1][1],nodes[2][1]);
	
	int num_qp=49; // the number of numerical intergation points
	double gauss[num_qp][3];
	init_Gauss(num_qp, 2, gauss); // gauss intergation initial	
	
	for(i=0;i<3;i++)
		b[i]=0;
	
	for(i=0;i<num_qp;i++)
	{
		x=nodes[0][0]*gauss[i][0]+nodes[1][0]*gauss[i][1]+nodes[2][0]*(1-gauss[i][0]-gauss[i][1]);
		y=nodes[0][1]*gauss[i][0]+nodes[1][1]*gauss[i][1]+nodes[2][1]*(1-gauss[i][0]-gauss[i][1]);
		a=f(x,y);
		
		b[0]+=2*s*gauss[i][2]*a*gauss[i][0];
		b[1]+=2*s*gauss[i][2]*a*gauss[i][1];
		b[2]+=2*s*gauss[i][2]*a*(1-gauss[i][0]-gauss[i][1]);
		
	}
}
