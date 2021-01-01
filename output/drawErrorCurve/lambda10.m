clear;
% theta=0.1
dofs1=[743	800	1498	2285	3243	4055	6406	10075	18082	30019	45931	63389	107067	198098	317003	602845];
sigmmaerr1=[7.5707E-01	7.4546E-01	3.4679E-01	1.3765E-01	5.4404E-02	3.0466E-02	1.3036E-02	4.6526E-03	1.5892E-03	5.8739E-04	2.4482E-04	1.1203E-04	4.0790E-05	1.2099E-05	4.9006E-06	1.4347E-06];
eta1=[2.4645E+01	1.9972E+01	8.5752E+00	3.4154E+00	1.6051E+00	1.0796E+00	4.5738E-01	1.7434E-01	5.9425E-02	2.1987E-02	8.9456E-03	4.4044E-03	1.5194E-03	4.5714E-04	1.6753E-04	5.1426E-05];

% theta=0.2
dofs2=[743	1028	2011	3628	5238	7824	10403	17327	24252	42341	63222	129624	263097	379509 654913];
sigmmaerr2=[7.5707E-01	7.2958E-01	3.4319E-01	1.1087E-01	5.0939E-02	1.6900E-02	7.8885E-03	2.5824E-03	1.2168E-03	3.9682E-04	1.7643E-04	5.1662E-05	2.0922E-05	7.0698E-06 3.0248E-06];
eta2=[2.4645E+01	1.3851E+01	6.3325E+00	3.3486E+00	1.1054E+00	5.5673E-01	2.1398E-01	8.8305E-02	3.5031E-02	1.3529E-02	4.9920E-03	1.6655E-03	4.4608E-04	2.2155E-04 6.6240E-05];

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

line([1 1 2 1]+9,[0 -2 -2 0]-8,'LineWidth',1,'Color',[0 0 0]);
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

