function [ux uy] = gradientofpolarfun(u, r ,t)
% The gradient in Cartesian coordinate system of polar function u=u(r,t)
% All the functions are expressed in symbolic form
ur=diff(u,r,1); %ur=simplify(ur);
ut=diff(u,t,1); %ut=simplify(ut);
ux=ur*cos(t)-ut*sin(t)/r; %ux=simplify(ux);
uy=ur*sin(t)+ut*cos(t)/r; %uy=simplify(uy);
end

