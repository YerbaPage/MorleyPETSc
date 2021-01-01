function showmesh(node,elem,uh)
% SHOWMESH plots the current mesh determined by node and elem
%
% USAGE
%    showmesh(node,elem)
%
% INPUT 
%     node: coordinate array of all nodes
%     elem: elements 
%

% L. Chen & C. Zhang 11-11-2006

trisurf(elem,node(:,1),node(:,2),zeros(size(node,1),1),'facecolor','none')
view(2); axis equal; axis off; drawnow;
%figure;
%trisurf(elem,node(:,1),node(:,2),uh)