clear;
% theta=0.1
dofs1=[743	857	1776	3770	6755	10801	21373	30293	50589	70545	126510	205513	360599	576462 950787];
sigmmaerr1=[9.3499E+00	6.6077E+00	2.3187E+00	5.9492E-01	1.9327E-01	7.8936E-02	2.6084E-02	1.1805E-02	3.8398E-03	2.1104E-03	6.3940E-04	2.3971E-04	8.2934E-05	3.5378E-05 1.2151E-05];
eta1=[3.8275E+02	2.2496E+02	7.6670E+01	1.7305E+01	6.0373E+00	2.3945E+00	8.3277E-01	3.5031E-01	1.2395E-01	6.7237E-02	2.1504E-02	6.9832E-03	2.6974E-03	1.0319E-03 3.9208E-04];

% theta=0.2
dofs2=[743	971	2862	8800	14612	21993	37995	58711	100405	138323	372505	587961	931507];
sigmmaerr2=[9.3499E+00	5.7291E+00	7.8895E-01	1.7232E-01	5.9566E-02	2.8496E-02	9.5411E-03	4.0095E-03	1.3510E-03	6.6553E-04	1.7188E-04	7.7219E-05	2.6162E-05];
eta2=[3.8275E+02	1.9788E+02	2.5454E+01	4.2727E+00	1.8549E+00	7.4281E-01	2.9971E-01	1.0289E-01	4.1408E-02	1.8734E-02	4.9523E-03	1.6077E-03	7.6118E-04];

dofs1=log(dofs1);
sigmmaerr1=log(sigmmaerr1);
eta1=log(eta1);
dofs2=log(dofs2);
sigmmaerr2=log(sigmmaerr2);
eta2=log(eta2);

plot(dofs1,sigmmaerr1,'o-',dofs1,eta1,'<-','LineWidth',1);
hold on
plot(dofs2,sigmmaerr2,'+-',dofs2,eta2,'*-','LineWidth',1);
hold off
leg=legend('$\ln\|\sigma-\sigma_h\|_A~(\theta=0.1)$ ','$\ln\eta(\sigma_h, \mathcal{T}_h)~(\theta=0.1)$', '$\ln\|\sigma-\sigma_h\|_A~(\theta=0.2)$ ','$\ln\eta(\sigma_h, \mathcal{T}_h)~(\theta=0.2)$');
%plot(dofs,sigmmaerr,'o-',dofs,eta,'<-',h,erroruE,'*-','LineWidth',2);
%leg=legend('$\ln\|\mathbf{\sigma}-\mathbf{\sigma}_h\|_0$','$\ln|u-u_h|_{1}$','$\ln|\!|\!|u-u_h|\!|\!|$');
%leg=legend('$\|\mathbf{\sigma}-\mathbf{\sigma}_h\|_0$','$|u-u_h|_{1}$','$|\!|\!|u-u_h|\!|\!|$','Location','southwest');
set(leg,'Interpreter','latex');
set(xlabel('$\ln(\#\textrm{dofs})$'),'Interpreter','latex');
%set(ylabel('$\ln$Errors'),'Interpreter','latex');

%set (gcf,'Position',[200, 100, 720, 540])
%set (gcf,'Position',[200, 100, 800, 600], 'color','w')
% reset: set (gcf,'Position',[232, 246, 560, 420], 'color','w')


%line([1 1 2 1]+3,[-4 -5 -5 -4]-2.5,'LineWidth',2,'Color',[0 0 0]);

line([1 1 2 1]+9.5,[0 -2 -2 0]-6,'LineWidth',1,'Color',[0 0 0]);
%line([3.5 2.5 3.5 3.5],[-3 -3 -4 -3],'LineWidth',2,'Color',[0 0 0]);

%line([2.5 2.5 3.5 2.5],[-1 -2 -2 -1],'LineWidth',2,'Color',[0 0 0]);

%line([6 6 8 6],[0 -1 -1 0],'LineWidth',2,'Color',[0 0 0]);



%plot(h,errorEk3l0,'p-',h,errorEk3l1,'<-',h,errorEk3l2,'o-',h,errorEk3l3,'s-','LineWidth',2);
%line([1.5 1.5 2.5 1.5],[-6 -8 -8 -6],'LineWidth',2,'Color',[0 0 0]);
%legend('l=0', 'l=1', 'l=2', 'l=3');

%subplot(2,1,1)
%plot(h,error1k2l1,'p-',h,error1k2l0,'<-',h,error1k1l0,'o-','LineWidth',2);
%legend('k=2,l=1', 'k=2,l=0', 'k=1,l=0', 'Location', 'SouthWest');
%grid on;
%subplot(2,1,2)
%plot(h,errorEk2l1,'p-',h,errorEk2l0,'<-',h,errorEk1l0,'o-','LineWidth',2);
%legend('k=2,l=1', 'k=2,l=0', 'k=1,l=0', 'Location', 'SouthWest');

grid on;

