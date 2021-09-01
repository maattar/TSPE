clc
clear
close all

%% Figure 6

load('bm.mat')
load('Cy.mat')
load('Cw.mat')
load('A1.mat')
load('A2.mat')
load('alfa.mat')
load('gama.mat')
load('pssi.mat')

t = (1200:10:1650)';

figure(6)
subplot(2,1,1)
plot(t,100*((ACy./Abm)-1),'vr','LineWidth',1.25)
hold on
plot(t,100*((Aa./Abm)-1),'xm','LineWidth',1.25)
plot(t,100*((Ag./Abm)-1),'+c','LineWidth',1.25)
plot(t,100*((Ap./Abm)-1),'>b','LineWidth',1.25)
plot(t,100*((ACw./Abm)-1),'-.','Color',[0.5 0.5 0.5],'LineWidth',1.5)
plot(t(11:46,1),100*((AA1./Abm(11:46,1))-1),'-','Color',[0.5 0.5 0.5],'LineWidth',1.5)
plot(t(11:46,1),100*((AA2./Abm(11:46,1))-1),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
hold off
grid on
title('Technology')
ylabel('percentage deviation from the benchmark')
subplot(2,1,2)
plot(t,100*((sCy./sbm)-1),'vr','LineWidth',1.25)
hold on
plot(t,100*((sa./sbm)-1),'xm','LineWidth',1.25)
plot(t,100*((sg./sbm)-1),'+c','LineWidth',1.25)
plot(t,100*((sp./sbm)-1),'>b','LineWidth',1.25)
plot(t,100*((sCw./sbm)-1),'-.','Color',[0.5 0.5 0.5],'LineWidth',1.5)
plot(t(11:46,1),100*((sA1./sbm(11:46,1))-1),'-','Color',[0.5 0.5 0.5],'LineWidth',1.5)
plot(t(11:46,1),100*((sA2./sbm(11:46,1))-1),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
hold off
ylim([-40 80])
legend('Clark: Income p.c. (DE)','\alpha, 25% higher','\gamma, 25% higher','\psi=1','Clark: real wage','Allen: real wage (unskilled)','Allen: real wage (skilled)','Location','NorthEast','Orientation','vertical','FontSize',8)
grid on
title('Survival')

%% Figure 7

load('enddate.mat')
ted = (1200:10:1800)';

load('enddateQ.mat')
tedQ = (1200:25:1800)';

figure(7)
subplot(2,1,1)
plot(ted,Aed,'o','MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
hold on
plot(t,Abm,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1.5,'MarkerSize',5)
plot(tedQ,AQ,'o-','LineWidth',1.5,'Color',[153/256 255/256 153/256])
hold off
grid on
title('Technology')
ylabel('1200=1')
legend('Decennial, T=1800','Decennial, T=1650','Quadranscentennial, T=1800')
subplot(2,1,2)
plot(ted,sed,'s','MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
hold on
plot(t,sbm,'s','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'LineWidth',1.5,'MarkerSize',5)
plot(tedQ,sQ,'s-','LineWidth',1.5,'Color',[255/256 153/256 255/256])
hold off
grid on
ylim([0.3 1])
title('Survival')
ylabel('1600=0.67')
legend('Decennial, T=1800','Decennial, T=1650','Quadranscentennial, T=1800')



