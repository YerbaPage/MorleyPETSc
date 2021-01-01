clear all
clc
load A0.dat;
load A1.dat;
load A2.dat;
load b.dat;
nr=max(A0(:,2));
i=[A0(:,1);A1(:,1);A2(:,1)+nr];
j=[A0(:,2);A1(:,2)+nr;A2(:,2)];
s=[A0(:,3);A1(:,3);A2(:,3)];
clear A0 A1 A2;
A=sparse(i,j,s);
clear i j s;
b=[zeros(nr,1);b];
u=A\b;
dlmwrite('u.dat',u,'\n')
