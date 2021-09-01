%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%  
%%%%%                   Calibration and Identification                %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cleaning
clc; clear; close all;

%% Data

Z = xlsread('Clark_Quadranscentennial.xls'); % Quadranscentennial data

%% Horizon

%Thor = 21; % Tend = 1650
Thor = 27; % Tend = 1800

D = Z(1:Thor,:);                                                         
clear Z;                                   
                                  
d = D(:,1); % periods         
T = size(d,1); % horizon

% tau1 = 40;                         % auxiliary times indexes for figures
% tau2 = 36;                         
% tau3 = 39;
% tau4 = 26;
% tau5 = 8;

%% Observed Variables

P = D(:,2);      % population (million)
y = D(:,3)/100;  % income per capita (1860s=1)

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

pssi = 0.6/0.672; % calculated from Bar & Leukhina's (2010, JEG) mort. data
alfa = 0.463;     % Bar and Leukhina (2010, RED)
gama = 0.418;     % Bar and Leukhina (2010, RED)

tau1 = 16;

gP_g = gP(tau1,1);                  
sp_g = 0.67;                       % s(1600) = 0.67

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

AQ=A; sQ=s;
save('enddateQ.mat','AQ','sQ')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%