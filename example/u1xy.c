  t2 = sin(t);
  t3 = z+1.0;
  t4 = z-1.0;
  t5 = 1.0/mu;
  t6 = t*t3;
  t7 = omega*t4;
  t8 = cos(t7);
  t9 = 1.0/t8;
  t10 = omega*t3;
  t11 = cos(t10);
  t12 = t*t4;
  t13 = lambda+mu;
  t14 = 1.0/t13;
  t15 = lambda*2.0;
  t16 = mu*4.0;
  t17 = t15+t16;
  t23 = t14*t17;
  t18 = -t23+z+1.0;
  t19 = sin(t6);
  t20 = t3*t3;
  t21 = t19*t20;
  t22 = sin(t12);
  t33 = t4*t9*t11*t18*t22;
  t24 = t21-t33;
  t25 = 1.0/r;
  t26 = cos(t);
  t27 = pow(r,z);
  t28 = cos(t6);
  t29 = cos(t12);
  t30 = pow(r,t4);
  t31 = t3*t28;
  t32 = t31-t9*t11*t18*t29;
  t0 = -t2*(1.0/(r*r)*t2*t5*t24*t27*(-1.0/2.0)+t2*t5*t24*t25*t30*z*(1.0/2.0)+pow(r,z-2.0)*t4*t5*t26*t32*z*(1.0/2.0))-t25*t26*(t2*t5*t25*t27*(t3*t20*t28-(t4*t4)*t9*t11*t18*t29)*(1.0/2.0)+t5*t24*t25*t26*t27*(1.0/2.0)-t2*t5*t30*t32*z*(1.0/2.0)-t5*t24*t26*t30*z*(1.0/2.0));
