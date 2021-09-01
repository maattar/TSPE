%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%  
%%%%%                   Calibration and Identification                %%%%%
%%%%%                                                                 %%%%%
%%%%%                 with Broadberry et al. (2015) data              %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning:

clc; clear; close all;

%% Data

Z = xlsread('Broadberry_Annual.xls'); % Annual data

Qd = a2d(Z);
clear Z
Qd(end,:) = [];
d = (1270:10:1860)';

Qe(:,2) = (4.87/Qd(1,2))*Qd(:,2);
Qe(:,3) = (100/Qd(60,3))*Qd(:,3); % 1st re-scale: for 1860s = 100

Qe(:,1) = [];
clear Qd

%% Horizon

Thor = 40;                     % Tend = 1650
%Thor = 54;                     % Tend = 1800

D = Qe(1:Thor,:); 
d = d(1:Thor,:);
clear Qe                                   
          
T = size(d,1);                     % horizon

tau1 = 33;                         % auxiliary times indexes for figures
% tau2 = 36;                         
% tau3 = 39;
% tau4 = 26;
% tau5 = 8;

%% Observed Variables

P = D(:,1);                        % Population                  (millions)
y = D(:,2)/100;                    % 2nd re-scale: for 1860s = 1

%% Population Growth Rate 

gP = zeros(T-1,1);

for t=1:T-1
    gP(t,1) = P(t+1,1)/P(t,1);
end

%% Horizon (adjusted)              

T  = T-1;
d  = d(1:T,1);
y  = y(1:T,1);
P  = P(1:T,1);

%% Calibration

pssi = 0.6/0.672;              % calculated from Bar and Leukhina's (2010, JEG) mortality data
alfa = 0.463;                  % Bar and Leukhina (2010, RED) 
gama = 0.418;                  % Bar and Leukhina (2010, RED)


gP_g = gP(tau1,1);                  
sp_g = 0.6;                       

w_g  = y(tau1,1);
wp_g = y(tau1+1,1);

rhho = (sqrt(((gP_g-pssi*sp_g)*gama*w_g)^2+4*pssi*gP_g*(gama^2)*sp_g*w_g*wp_g)-(gP_g-pssi*sp_g)*gama*w_g)/(2*pssi*gP_g);

A_x  = 1;                          % A(1200) = 1
w_x  = y(1,1);
P_x  = P(1,1);

X    = ((w_x^(1/alfa))*pssi*P_x)/(A_x*(pssi+(gama/rhho)*w_x));

%% Identification of A_t and s_t

A = NaN(T-1,1);
s = NaN(T-1,1);

for t = 2:T-1
    A(t,1) = ((y(t,1)^(1/alfa))*pssi*P(t,1))/(X*(pssi+(gama/rhho)*y(t,1)));
    s(t,1) = (gP(t-1,1)*rhho*(rhho*pssi+gama*y(t-1,1)))/(gama*y(t-1,1)*(rhho*pssi+gama*y(t,1))) ;
end

A(1,1) = 1;

d2 = d(1:T-1,1);

%% Figure 8

figure(8)
subplot(2,1,1)
[AX,H1,H2] = plotyy(d,P,d,y,'plot');
set(get(AX(1),'Ylabel'),'String','million') 
set(get(AX(2),'Ylabel'),'String','1860s=1')
set(AX(1),'YColor','blue') 
set(AX(2),'YColor','red') 
set(AX(1),'XLim',[1260 1660])
set(AX(2),'XLim',[1260 1660])
set(AX(1),'YLim',[1 6])
set(AX(2),'YLim',[0.1 0.35])
set(H1,'LineStyle','-','Color','blue','Marker','o','MarkerEdgeColor','blue','MarkerFaceColor',[204/256 229/256 255/256],'LineWidth',1.5,'MarkerSize',5)
set(H2,'LineStyle','-','Color','red','Marker','s','MarkerEdgeColor','red','MarkerFaceColor',[255/256 204/256 204/256],'LineWidth',1.5,'MarkerSize',5)
grid on
legend('Population (left axis)','GDP per capita (right axis)','Location','SouthEast','Orientation','horizontal')
subplot(2,1,2)
[AX,H1,H2] = plotyy(d2,A,d2,s,'plot');
set(get(AX(1),'Ylabel'),'String','1270=1') 
set(get(AX(2),'Ylabel'),'String','1600=0.6')
set(AX(1),'YColor',[0 152/256 0]) 
set(AX(2),'YColor',[193/256 0 193/256]) 
set(AX(1),'XLim',[1260 1660])
set(AX(2),'XLim',[1260 1660])
set(H1,'LineStyle','none','Marker','o','MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
set(H2,'LineStyle','none','Marker','s','MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
grid on
legend('Technology (left axis)','Survival (right axis)','Location','SouthEast','Orientation','horizontal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%