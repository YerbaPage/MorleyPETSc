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
  t17 = t13*t16;
  t18 = t17+z-1.0;
  t0 = pow(r,t3)*t4*z*cos(t)*(t2*sin(t5)-t8*t10*t18*sin(t11))*(1.0/2.0)-(pow(r,z)*t4*sin(t)*((t2*t2)*cos(t5)-t3*t8*t10*t18*cos(t11))*(1.0/2.0))/r;