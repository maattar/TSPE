%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%  
%%%%%                    The Impulse-Response Analysis                %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning:

clc; clear; close all;

%% Data

load('irfsdata.mat'); 

D = [A s d2];
D(1,:) = [];
T = size(D,1);
A=D(:,1); s=D(:,2); d=D(:,3);

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

simT = 9;

%% The Initial Steady-State

As = mean(A);
ss = mean(s);
ws = rhho/(gama*ss);
bs = 1/ss;
Ls = As*X*(1/pssi)*((gama/rhho)^(1/alfa))*(ss^((1-alfa)/alfa));
Ps = (pssi+(1/ss))*Ls*ss;

%% Responses to the Technology Shock

wsim_A  = zeros(simT,1);    wsim_A(1,1) = ws;
bsim_A  = zeros(simT,1);    bsim_A(1,1) = bs;
Psim_A  = zeros(simT,1);    Psim_A(1,1) = Ps;
Lsim_A  = zeros(simT,1);    Lsim_A(1,1) = Ls;
Asim_A  = zeros(simT,1);    Asim_A(1,1) = As;
ssim_A  = zeros(simT,1);    ssim_A(1,1) = ss;
epsA_A  = zeros(simT,1);    epsA_A(1,1) = sigA;
epss_A  = zeros(simT,1);    epss_A(1,1) = 0;

for z = 1:simT-1
    Asim_A(z+1,1) = Asim_A(z,1) + epsA_A(z,1);
    ssim_A(z+1,1) = ssim_A(z,1) + epss_A(z,1);
    wsim_A(z,1)   = ((Asim_A(z,1)*X)/(pssi*ssim_A(z,1)*Lsim_A(z,1)))^(alfa); 
    bsim_A(z,1)   = (gama/rhho)*wsim_A(z,1);
    Psim_A(z,1)   = (pssi+bsim_A(z,1))*Lsim_A(z,1)*ssim_A(z,1);
    Lsim_A(z+1,1) = ssim_A(z,1)*bsim_A(z,1)*Lsim_A(z,1);
end

%% Responses to the Survival Shock

wsim_s  = zeros(simT,1);    wsim_s(1,1) = ws;
bsim_s  = zeros(simT,1);    bsim_s(1,1) = bs;
Psim_s  = zeros(simT,1);    Psim_s(1,1) = Ps;
Lsim_s  = zeros(simT,1);    Lsim_s(1,1) = Ls;
Asim_s  = zeros(simT,1);    Asim_s(1,1) = As;
ssim_s  = zeros(simT,1);    ssim_s(1,1) = ss;
epsA_s  = zeros(simT,1);    epsA_s(1,1) = 0;
epss_s  = zeros(simT,1);    epss_s(1,1) = sigs;

for z = 1:simT-1
    Asim_s(z+1,1) = Asim_s(z,1) + epsA_s(z,1);
    ssim_s(z+1,1) = ssim_s(z,1) + epss_s(z,1);
    wsim_s(z,1)   = ((Asim_s(z,1)*X)/(pssi*ssim_s(z,1)*Lsim_s(z,1)))^(alfa); 
    bsim_s(z,1)   = (gama/rhho)*wsim_s(z,1);
    Psim_s(z,1)   = (pssi+bsim_s(z,1))*Lsim_s(z,1)*ssim_s(z,1);
    Lsim_s(z+1,1) = ssim_s(z,1)*bsim_s(z,1)*Lsim_s(z,1);
end

%% Figure 4

figure(4)
subplot(1,2,1)
plot((1:simT-2)',ws+zeros(simT-2,1),'--','Color','black','LineWidth',2)
hold on
plot((1:simT-2)',wsim_A(2:simT-1,1),'o-','Color',[0 102/256 0],'MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
plot((1:simT-2)',wsim_s(2:simT-1,1),'s-','Color',[153/256 0 153/256],'MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
hold off
grid on
xlim([0 8])
ylim([0.50 0.65])
title('Responses of Income per capita')
xlabel('decades after shock')
ylabel('1860s=1')
legend('Initial Steady-State','Technology Shock','Survival Shock','location','NorthEast');
subplot(1,2,2)
plot((1:simT-2)',Ps+zeros(simT-2,1),'--','Color','black','LineWidth',2)
hold on
plot((1:simT-2)',Psim_A(2:simT-1,1),'o-','Color',[0 102/256 0],'MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
plot((1:simT-2)',Psim_s(2:simT-1,1),'s-','Color',[153/256 0 153/256],'MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
hold off
grid on
xlim([0 8])
ylim([3.5 4.5])
title('Responses of Population')
ylabel('million')
xlabel('decades after shock')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%