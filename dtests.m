%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%  
%%%%%                         Diagnostic Tests                        %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning:

clc; clear; close all;

%% Data

load('dtestsdata.mat');          % Identified Technology and Survival Terms

D = [A s d2];
D(1,:) = [];
T = size(D,1);
A=D(:,1); s=D(:,2); d=D(:,3);

%% Unit Root Tests on A_t and s_t

% Phillips-Perron
[h_A_pp,pValue_A_pp,~,~,reg_A_pp] = pptest(A,'model','AR');
[h_s_pp,pValue_s_pp,~,~,reg_s_pp] = pptest(s,'model','AR');

% Augmented Dickey-Fuller
[h_A_adf,pValue_A_adf,~,~,reg_A_adf] = adftest(A,'model','AR');
[h_s_adf,pValue_s_adf,~,~,reg_s_adf] = adftest(s,'model','AR');

% KPSS
[h_A_kpss,pValue_A_kpss,~,~,reg_A_kpss] = kpsstest(A,'trend',false);
[h_s_kpss,pValue_s_kpss,~,~,reg_s_kpss] = kpsstest(s,'trend',false);

%% First-Differences

dA = zeros(T-1,1);
ds = zeros(T-1,1);
for t=1:T-1
    dA(t,:) = A(t+1,:)-A(t,:);
    ds(t,:) = s(t+1,:)-s(t,:);
end

dd = d(1:T-1,1);

sigma_A = std(dA);
sigma_s = std(ds);

%% Various Tests on d(A_t) and d(s_t)

% Z test
[h_dA_z,p_dA_z] = ztest(dA,0,std(dA));
[h_ds_z,p_ds_z] = ztest(ds,0,std(ds));

% Jarque-Bera test
[h_dA_jb,p_dA_jb] = jbtest(dA);
[h_ds_jb,p_ds_jb] = jbtest(ds);

% Lilliefors-van Soest test
[h_dA_l,p_dA_l] = lillietest(dA);
[h_ds_l,p_ds_l] = lillietest(ds);

% Kolmogorov-Smirnov test
[h_dA_ks,p_dA_ks] = kstest(dA/sigma_A);
[h_ds_ks,p_ds_ks] = kstest(ds/sigma_s);

% LBQ test
[h_dA_lbq,p_dA_lbq] = lbqtest(dA);
[h_ds_lbq,p_ds_lbq] = lbqtest(ds);

% Engle's ARCH test
[h_dA_arch,p_dA_arch] = archtest(dA);
[h_ds_arch,p_ds_arch] = archtest(ds);

%% Figure 3

figure(3)
subplot(2,1,1)
rectangle('Position',[1550,-0.3,100,0.6],'FaceColor',[0.9 0.9 0.9])
hold on
stem(dd,dA,'-o','Color',[0 102/256 0],'MarkerEdgeColor',[0 102/256 0],'MarkerFaceColor',[153/256 255/256 153/256],'LineWidth',1.5,'MarkerSize',5)
hold off
title('Technology Shocks')
grid on
box on
ylim([-0.3 0.3])
set(gca, 'Layer', 'top');
subplot(2,1,2)
rectangle('Position',[1550,-0.3,100,0.6],'FaceColor',[0.9 0.9 0.9])
hold on
stem(dd,ds,'-s','Color',[153/256 0 153/256],'MarkerEdgeColor',[153/256 0 153/256],'MarkerFaceColor',[255/256 153/256 255/256],'LineWidth',1.5,'MarkerSize',5)
hold off
title('Survival Shocks')
grid on
box on
ylim([-0.3 0.3])
set(gca, 'Layer', 'top');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%