  t2 = z+1.0;
  t3 = z-1.0;
  t4 = 1.0/mu;
  t5 = t*t2;
  t6 = omega*t3;
  t7 = cos(t6);
  t8 = 1.0/t7;
  t9 = omega*t2;
  t10 = cos(t9);
  t11 = t*t3;
  t12 = lambda+mu;
  t13 = 1.0/t12;
  t14 = lambda*2.0;
  t15 = mu*4.0;
  t16 = t14+t15;
  t23 = t13*t16;
  t17 = -t23+z+1.0;
  t18 = pow(r,t3);
  t19 = cos(t);
  t20 = cos(t5);
  t21 = t2*t20;
  t22 = cos(t11);
  t24 = t21-t8*t10*t17*t22;
  t25 = t4*t18*t19*t24*z*(1.0/2.0);
  t26 = sin(t);
  t27 = sin(t5);
  t28 = sin(t11);
  t29 = pow(r,z);
  t30 = 1.0/r;
  t31 = t2*t2;
  t32 = t23+z-1.0;
  t33 = t27*t31;
  t34 = t33-t3*t8*t10*t17*t28;
  t35 = t4*t26*t29*t30*t34*(1.0/2.0);
  t0 = mu*(t25+t35)*-2.0-lambda*(t25+t35-t4*t18*t26*z*(t2*t27-t8*t10*t28*t32)*(1.0/2.0)-t4*t19*t29*t30*(t20*t31-t3*t8*t10*t22*t32)*(1.0/2.0));