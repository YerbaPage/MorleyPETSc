%p=load('example2\nodestheta05iter13.dat');
%T=load('example2\elementstheta05iter13.dat');
%uh=load('example2\uhtheta05iter13.dat');
p=load('nodes.dat');
T=load('elements.dat');
uh=load('uh.dat');
showmesh(p,T,uh)

%text(p(:,1),p(:,2),num2str(p(:,3)),'Color','red','FontSize',14)