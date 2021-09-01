function run_Benchmark_Decennial(LivStd,alfa_t,gama_t,pssi_t,Tend_t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%  
%%%%%                   Calibration and Identification                %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Guide
%
%  Description:
%   
%   This program calibrates the structural parameters of the model and  
%   identifies technology and survival terms. 
%
%   The program uses Gregory Clark's (2010, REH) data, and it generates 
%   Figures 1 and 2 and Table 1.
%
%  Input Arguments:
%
%   LivStd - This input specifies the measure of living standards.
%
%         LivStd = 0 => NNIpc (NDP)      (Clark)                  *default*
%         LivStd = 1 => NNIpc (DE)       (Clark)
%         LivStd = 2 => real wage        (Clark)
%         LivStd = 3 => real wage        (Allen, unskilled)
%         LivStd = 4 => real wage        (Allen, skilled)
%
%   alfa_t - This input specifies the value of alpha.
%
%         alfa_t = 0 => alfa = 0.463                              *default*
%         alfa_t = 1 => alfa = 1.25*0.463
%
%   gama_t - This input specifies the value of gamma.
%
%         gama_t = 0 => gama = 0.418                              *default*
%         gama_t = 1 => gama = 1.25*0.418
%
%   pssi_t - This input specifies the value of gamma.
%
%         pssi_t = 0 => pssi = 0.6/0.672                          *default*
%         pssi_t = 1 => pssi = 1
%
%   Tend_t - This input specifies the end date of the Malthusian era.
%
%         Tend_t  = 0 => Tend = 1650                              *default*
%         Tend_t  = 1 => Tend = 1800
%
%  Examples:
%
%   Benchmark   :  run_Benchmark_Decennial(0,0,0,0,0)
%   Liv. Std.   :  run_Benchmark_Decennial(2,0,0,0,0)
%   Parameters  :  run_Benchmark_Decennial(0,1,1,1,0)
%   End Date    :  run_Benchmark_Decennial(0,0,0,0,1)
%
%  Copyright 2021 M. Aykut Attar
%
%  Date: April 2, 2021

%% Data

Z = xlsread('Benchmark_Decennial.xls'); % Decennial data

%% Horizon

if Tend_t == 0
    Thor = 48;                     % Tend = 1650
else
    Thor = 63;                     % Tend = 1800
end

D = Z(1:Thor,:);                                                         
clear Z;                                   
                                  
d = D(:,1);                        % decades         
T = size(d,1);                     % horizon

tau1 = 40;                         % auxiliary times indexes for figures
tau2 = 36;                         
tau3 = 39;
tau4 = 26;
tau5 = 8;

%% Observed Variables

P = D(:,2);                        % Population                   (million)

if LivStd == 0
    y = D(:,4)/100;                % NNIpc (NDP)                 (1860 = 1)
elseif LivStd == 1
    y = D(:,3)/100;                % NNIpc (DE)                  (1860 = 1)
elseif LivStd == 2
    y = D(:,5)/100;                % Real Wage                   (1860 = 1)
elseif LivStd == 3
    y = D(:,8)/100;                % Real Wage (Allen, u)        (1860 = 1)
else
    y = D(:,9)/100;                % Real Wage (Allen, s)        (1860 = 1)
end

if (LivStd == 3) || (LivStd == 4)
    P(1:10,:) = [];
    y(1:10,:) = [];
    d(1:10,:) = [];
    T = size(P,1);
    tau1 = 30;
    tau2 = 26;
    tau3 = 29;
    tau4 = 16;
    tau5 = 1;
    D(1:10,:) = [];
end

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

if alfa_t == 0
    alfa = 0.463;                  % Bar and Leukhina (2010, RED)  
else
    alfa = 1.25*0.463;
end

if gama_t == 0
    gama = 0.418;                  % Bar and Leukhina (2010, RED)
else
    gama = 1.25*0.418;
end

if pssi_t == 0
    pssi = 0.6/0.672;              % calculated from Bar and Leukhina's (2010, JEG) mortality data  
else
    pssi = 1; 
end

gP_g = gP(tau1,1);                  
sp_g = 0.67;                       % s(1600) = 0.67

if (LivStd == 3) || (LivStd == 4)
    sp_g = 0.6;
end

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

%% Efficiency and Survival

eff_o = D(1:T-1,6);
eff_e = D(1:T-1,7);
a_la0 = D(1:T-1,10);             % natural logarithm & non-scaled
a_la1 = exp(a_la0);              % levels & non-scaled
a_las = a_la1(tau2,1)/A(tau2,1); % re-scaling factor A(1550)=A_LA(1550)
a_la  = a_la1/a_las;
sr4   = D(1:T-1,11);
sr9   = D(1:T-1,12);
sr14  = D(1:T-1,13);
sr19  = D(1:T-1,14);
sr24  = D(1:T-1,15);
sbl5  = D(1:T-1,16);
sbl10 = D(1:T-1,17);
sbl15 = D(1:T-1,18);
sbl25 = D(1:T-1,19);

%% Figures 1 and 2

figure(1)
[AX,H1,H2] = plotyy(d,P,d,y,'plot');
set(get(AX(1),'Ylabel'),'String','million') 
set(get(AX(2),'Ylabel'),'String','1860s=1')
set(AX(1),'YColor','blue') 
set(AX(2),'YColor','red') 
set(H1,'LineStyle','-','Color','blue','Marker','o','MarkerEdgeColor','blue','MarkerFaceColor',[204/256 229/256 255/256],'LineWidth',1.5,'MarkerSize',5)
set(H2,'LineStyle','-','Color','red','Marker','s','MarkerEdgeColor','red','MarkerFaceColor',[255/256 204/256 204/256],'LineWidth',1.5,'MarkerSize',5)
grid on
legend('Population (left axis)','Income per capita (right axis)','Orientation','vertical','Location','SouthEast')



figure(2)
subplot(2,1,1)
rectangle('Position',[1200 1.1 450 0.25],'FaceColor','white','EdgeColor','white')
hold on
text(1335,1.3,'50-year average temperature in \circC','FontSize',8,'FontWeight','bold')
text(1220,1.2,'10.1','FontSize',8,'FontWeight','bold','Color',[243/255 156/255 18/255])
text(1270,1.2,'10.2','FontSize',8,'FontWeight','bold','Color',[230/255 126/255 34/255])
text(1320,1.2,'9.8','FontSize',8,'FontWeight','bold','Color',[240/255 200/255 111/255])
text(1370,1.2,'9.5','FontSize',8,'FontWeight','bold','Color',[125/255 206/255 160/255])
text(1420,1.2,'9.1','FontSize',8,'FontWeight','bold','Color',[133/255 193/255 233/255])
text(1470,1.2,'9.0','FontSize',8,'FontWeight','bold','Color',[46/255 134/255 193/255])
text(1520,1.2,'9.3','FontSize',8,'FontWeight','bold','Color',[133/255 213/255 200/255])
text(1570,1.2,'8.8','FontSize',8,'FontWeight','bold','Color',[33/255 16/255 172/255])
text(1620,1.2,'8.8','FontSize',8,'FontWeight','bold','Color',[33/255 16/255 172/255])
plot(d2,eff_o,'Color','red','LineWidth',2)
plot(d2,a_la,'-.','Color','black','LineWidth',2)
plot(d2,A,'o','MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
hold off
ylim([0.2 1.35])
yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('1200=1')
grid on
title('Identified Technology Term')
legend('Clark (2010)','Lee & Anderson (2002)','Technology','Location','SouthWest')
subplot(2,1,2)
rectangle('Position',[1200 1 450 0.6],'FaceColor','white','EdgeColor','white')
hold on
rectangle('Position',[1348 1 302 0.6],'FaceColor',[255/256 230/256 255/256],'EdgeColor',[255/256 230/256 255/256])
ciplot(sr24(tau5:tau4,1),sr4(tau5:tau4,1),d2(tau5:tau4,1),[0.9 0.9 0.9])
ciplot(sbl25(tau3:tau3+7,1),sbl5(tau3:tau3+7,1),d2(tau3:tau3+7,1),[0 0 0])
plot(d2,s,'s','MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
text(1315,1.595,'Great Famine','Rotation',90,'HorizontalAlignment','right','FontSize',8,'FontWeight','bold','Color',[153/256 0 153/256]);
text(1348,1.595,'Black Death','Rotation',90,'HorizontalAlignment','right','FontSize',8,'FontWeight','bold','Color',[153/256 0 153/256]);
text(1550,1.3,'Plague Centuries','HorizontalAlignment','right','FontSize',8,'FontWeight','bold','Color','black');
hold off
ylim([0 1.6])
yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('1600=0.67')
grid on
title('Identified Survival Term')
legend('Russell (1948)','Bar & Leukhina (2010a)','Survival','Location','SouthEast')

%% Saving benchmark A and s for diagnostic tests, irfs and the M cycle
% These runs should be executed with LivStd=alfa_t=gama_t=pssi_t=Tend_t=0.
save('dtestsdata.mat','A','s','d2')
save('irfsdata.mat','A','s','d2','X','rhho','pssi','gama','alfa')
save('mcycledata.mat','A','s','d2','P','y','X','rhho','pssi','gama','alfa')

%% Saving A and s for robustness analyses
% These runs should be executed one at a time.
% Abm=A; sbm=s;
% save('bm.mat','Abm','sbm') % LivStd=alfa_t=gama_t=Tend_t=0
% ACy=A; sCy=s;
% save('Cy.mat','ACy','sCy') % LivStd=1
% ACw=A; sCw=s;
% save('Cw.mat','ACw','sCw') % LivStd=2
% AA1=A; sA1=s;
% save('A1.mat','AA1','sA1') % LivStd=3
% AA2=A; sA2=s;
% save('A2.mat','AA2','sA2') % LivStd=4
% Aa=A; sa=s;
% save('alfa.mat','Aa','sa') % LivStd=0 alfa_t=1
% Ag=A; sg=s;
% save('gama.mat','Ag','sg') % alfa_t=0 gama_t=1
% Ap=A; sp=s;
% save('pssi.mat','Ap','sp') % alfa_t=0 gama_t=0 pssi_t=1
% Aed=A; sed=s;
% save('enddate.mat','Aed','sed') % alfa_t=0 gama_t=0 pssi_t=0 Tend_t=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%