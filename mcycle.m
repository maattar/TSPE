%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%  
%%%%%                The Source of the Malthusian Cycle               %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning:

clc; clear; close all;

%% Data

load('mcycledata.mat'); 

T = size(A,1);
P = P(1:T,:); y = y(1:T,:); 

D = [A s d2 P y];
D(1,:) = [];
T = size(D,1);
A=D(:,1); s=D(:,2); d=D(:,3); P=D(:,4); y=D(:,5);


%% First-Differences

dA = zeros(T-1,1);
ds = zeros(T-1,1);
for t=1:T-1
    dA(t,:) = A(t+1,:)-A(t,:);
    ds(t,:) = s(t+1,:)-s(t,:);
end

dd = d(1:T-1,1);

sigA = std(dA);
sigs = std(ds);

%% Simulation Horizon

simT = size(dd,1);

%% The Initial Steady-State

As = mean(A);
ss = mean(s);
ws = rhho/(gama*ss);
bs = 1/ss;
Ls = As*X*(1/pssi)*((gama/rhho)^(1/alfa))*(ss^((1-alfa)/alfa));
Ps = (pssi+(1/ss))*Ls*ss;

%% Without Technology Shocks

wwA  = zeros(simT,1);    wwA(1,1)  = ws;
bwA  = zeros(simT,1);    bwA(1,1)  = bs;
PwA  = zeros(simT,1);    PwA(1,1)  = Ps;
LwA  = zeros(simT,1);    LwA(1,1)  = Ls;
AwA  = zeros(simT,1);    AwA(1,1)  = As;
swA  = zeros(simT,1);    swA(1,1)  = ss;
eAwA = zeros(simT,1);    
eswA = ds;    

for z = 1:simT
    AwA(z+1,1) = AwA(z,1) + eAwA(z,1);
    swA(z+1,1) = swA(z,1) + eswA(z,1);
    wwA(z,1)   = ((AwA(z,1)*X)/(pssi*swA(z,1)*LwA(z,1)))^(alfa); 
    bwA(z,1)   = (gama/rhho)*wwA(z,1);
    PwA(z,1)   = (pssi+bwA(z,1))*LwA(z,1)*swA(z,1);
    LwA(z+1,1) = swA(z,1)*bwA(z,1)*LwA(z,1);
end

%% Without Survival Shocks

wws  = zeros(simT,1);    wws(1,1)  = ws;
bws  = zeros(simT,1);    bws(1,1)  = bs;
Pws  = zeros(simT,1);    Pws(1,1)  = Ps;
Lws  = zeros(simT,1);    Lws(1,1)  = Ls;
Aws  = zeros(simT,1);    Aws(1,1)  = As;
sws  = zeros(simT,1);    sws(1,1)  = ss;
eAws = dA;    
esws = zeros(simT,1);    

for z = 1:simT
    Aws(z+1,1) = Aws(z,1) + eAws(z,1);
    sws(z+1,1) = sws(z,1) + esws(z,1);
    wws(z,1)   = ((Aws(z,1)*X)/(pssi*sws(z,1)*Lws(z,1)))^(alfa); 
    bws(z,1)   = (gama/rhho)*wws(z,1);
    Pws(z,1)   = (pssi+bws(z,1))*Lws(z,1)*sws(z,1);
    Lws(z+1,1) = sws(z,1)*bws(z,1)*Lws(z,1);
end

%% Figure 5

figure(5)
subplot(2,2,1)
plot(dd,P(2:simT+1,1),'Color',[0.375 0.375 0.375],'LineWidth',2)
hold on
plot(dd,PwA,'o','MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
hold off
grid on
ylabel('million')
title('Population')
ylim([1.5 8])
%xlim([d3(2,1) d3(simT,1)])
subplot(2,2,2)
plot(dd,y(2:simT+1,1),'Color',[0.375 0.375 0.375],'LineWidth',2)
hold on
plot(dd,wwA,'o','MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
hold off
grid on
ylabel('1860s=1')
title('Income per capita')
ylim([0.1 0.9])
% xlim([d3(2,1) d3(simT,1)])
legend('data','w / o Technology Shocks','Location','SouthWest','Orientation','vertical')
subplot(2,2,3)
plot(dd,P(2:simT+1,1),'Color',[0.375 0.375 0.375],'LineWidth',2)
hold on
plot(dd,Pws,'s','MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
hold off
grid on
ylabel('million')
title('Population')
ylim([1.5 8])
% xlim([d3(2,1) d3(simT,1)])
subplot(2,2,4)
plot(dd,y(2:simT+1,1),'Color',[0.375 0.375 0.375],'LineWidth',2)
hold on
plot(dd,wws,'s','MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
hold off
grid on
ylabel('1860s=1')
title('Income per capita')
ylim([0.1 0.9])
% xlim([d3(2,1) d3(simT,1)])
legend('data','w / o Survival Shocks','Location','SouthWest','Orientation','vertical')



%% Correlations

[r_data,p_data] = corr(P(2:simT,1),y(2:simT,1));

[r_wA,p_wA]     = corr(PwA,wwA);

[r_ws,p_ws]     = corr(Pws,wws);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%