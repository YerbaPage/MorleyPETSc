elem=[3 6 1
4 1 6
6 8 4
7 4 8
4 7 2
5 2 7];
node=[-1 -1 12
-1 1 12
0 -2 1
0 0 22
0 2 1
1 -1 1
1 1 1
2 0 1
];
trisurf(elem,node(:,1),node(:,2),zeros(size(node,1),1),'facecolor','none')
view(2); axis equal; axis off; drawnow;
axis on
