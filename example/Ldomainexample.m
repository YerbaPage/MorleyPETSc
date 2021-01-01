clear;
syms r t z omega lambda mu;
x = r*cos(t);
y = r*sin(t);
C1 = -cos((z+1)*omega)/cos((z-1)*omega);
C2 = 2*(lambda+2*mu)/(lambda+mu);

u1 = (-(z+1)*cos((z+1)*t) + (C2-z-1)*C1*cos((z-1)*t))*r^z/(2*mu);
u2 = ((z+1)*sin((z+1)*t) + (C2+z-1)*C1*sin((z-1)*t))*r^z/(2*mu);

[u1x, u1y] = gradientofpolarfun(u1, r, t);
[u2x, u2y] = gradientofpolarfun(u2, r, t);
[u1xx, u1xy] = gradientofpolarfun(u1x, r, t);
[u1yx, u1yy] = gradientofpolarfun(u1y, r, t);
[u2xx, u2xy] = gradientofpolarfun(u2x, r, t);
[u2yx, u2yy] = gradientofpolarfun(u2y, r, t);

epsu = [u1x (u1y+u2x)/2;(u1y+u2x)/2 u2y];
sigma = 2*mu*epsu + lambda*trace(epsu)*eye(2);

[sigma11x, sigma11y] = gradientofpolarfun(sigma(1,1), r, t);
[sigma12x, sigma12y] = gradientofpolarfun(sigma(1,2), r, t);
[sigma21x, sigma21y] = gradientofpolarfun(sigma(2,1), r, t);
[sigma22x, sigma22y] = gradientofpolarfun(sigma(2,2), r, t);

f1 = -(sigma11x + sigma12y); 
f2 = -(sigma21x + sigma22y);

f1 = simplify(f1);
f2 = simplify(f2);

t1=sym([1;-1]/sqrt(2));
t2=sym([1;1]/sqrt(2));
t3=sym([-1;1]/sqrt(2));
t4=sym([-1;-1]/sqrt(2));
n1=[0 1;-1 0]*t1;
n2=[0 1;-1 0]*t2;
n3=[0 1;-1 0]*t3;
n4=[0 1;-1 0]*t4;

u1Grad=[u1x, u1y];
u2Grad=[u2x, u2y];
ut_t_1 = [u1Grad*t1, u2Grad*t1]*t1;
ut_t_2 = [u1Grad*t2, u2Grad*t2]*t2;
ut_t_3 = [u1Grad*t3, u2Grad*t3]*t3;
ut_t_4 = [u1Grad*t4, u2Grad*t4]*t4;
u1Hess=[u1xx, u1xy;u1yx, u1yy];
u2Hess=[u2xx, u2xy;u2yx, u2yy];
un_tt_1 = [t1'*u1Hess*t1, t1'*u2Hess*t1]*n1;
un_tt_2 = [t2'*u1Hess*t2, t2'*u2Hess*t2]*n2;
un_tt_3 = [t3'*u1Hess*t3, t3'*u2Hess*t3]*n3;
un_tt_4 = [t4'*u1Hess*t4, t4'*u2Hess*t4]*n4;

ut_t_1 = simplify(ut_t_1);
ut_t_2 = simplify(ut_t_2);
ut_t_3 = simplify(ut_t_3);
ut_t_4 = simplify(ut_t_4);

un_tt_1 = simplify(un_tt_1);
un_tt_2 = simplify(un_tt_2);
un_tt_3 = simplify(un_tt_3);
un_tt_4 = simplify(un_tt_4);

% write to files
ccode(f1,'file','f1.c')
ccode(f2,'file','f2.c')
ccode(u1,'file','u1.c')
ccode(u2,'file','u2.c')
ccode(u1x,'file','u1x.c')
ccode(u1y,'file','u1y.c')
ccode(u2x,'file','u2x.c')
ccode(u2y,'file','u2y.c')
ccode(sigma(1,1),'file','simga11.c')
ccode(sigma(1,2),'file','simga12.c')
ccode(sigma(2,1),'file','simga21.c')
ccode(sigma(2,2),'file','simga22.c')
ccode(u1xx,'file','u1xx.c')
ccode(u1xy,'file','u1xy.c')
ccode(u1yx,'file','u1yx.c')
ccode(u1yy,'file','u1yy.c')
ccode(u2xx,'file','u2xx.c')
ccode(u2xy,'file','u2xy.c')
ccode(u2yx,'file','u2yx.c')
ccode(u2yy,'file','u2yy.c')

ccode(ut_t_1,'file','ut_t_1.c')
ccode(ut_t_2,'file','ut_t_2.c')
ccode(ut_t_3,'file','ut_t_3.c')
ccode(ut_t_4,'file','ut_t_4.c')

ccode(un_tt_1,'file','un_tt_1.c')
ccode(un_tt_2,'file','un_tt_2.c')
ccode(un_tt_3,'file','un_tt_3.c')
ccode(un_tt_4,'file','un_tt_4.c')

