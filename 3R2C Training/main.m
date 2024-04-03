%	Example of parameter estimation and online disturbance observer
% 
%	DC house thermal model
% 
%	Authors: Trevor Bird and Elias N Pergantis
% Please cite: 
% % E. N. Pergantis, Priyadarshan, N. A. Theeb, P. Dhillon,
% J. P. Ore, D. Ziviani, E. A. Groll, K. J. Kircher, Field
% demonstration of predictive heating control for an all-
% electric house in a cold climate, Applied Energy 360
% (2024) 122820
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% load the excel sheet

% x_dot = A*x + B*u + W*(vbar + vhat)
% y = C*x

% d = vbar + vhat - some offset with varying value

% States
% x(1) = T_in % indoor air temperature [C]
% x(2) = T_m % mass temperature [C]

% Inputs
% u(1) = T_out % outdoor air temperature [C] 
% u(2) = Q_hp % Heat Pump supply heat [kW]
% u(3) = I/1000 % Global horizontal radiation [kW/m^2] the /1000 is to have
% all powers in the correct units

% Observables
% y(1) = T_in, 1x1 matrix, T_m is not observed

% d -> Disturbance term
% This term is to be obtained from the fit for ML fitting purposes
% d = f(Hour_of_day, Day_of_week, wind_speed, T_out)

%% Load the data

% I saved the data in a .mat so it would be faster to load..
data = readtimetable('DC_House_Nov_April.xlsx');
data = rmmissing(data);
tdata = 1:length(data.Time);	% time
Y = data.IndoorTemp_C_;
U = [data.OutdoorTemp_C_ data.HeatSupply_kW_ data.SolarIrradiation_W_m_2_/1000];

% all of the data
% tstart = 1;
% tend = size(Y,1);
% trim some of the data off..
tstart = 500;
tend = 1000;

tdata = tdata(tstart:tend);
tdata = tdata - tdata(1);
Y = Y(tstart:tend);
U = U(tstart:tend,:);

% sample time is one hour
Ts = 1;

figure; hold on;
plot(tdata,Y) % Indoor Temperature

figure; hold on;
plot(tdata,U(:,1)) % Outdoor temperature
figure; hold on;
plot(tdata,U(:,2)) % Heat supply
figure; hold on;
plot(tdata,U(:,3)) % Solar irradiation 

%% simulate the system for an initial guess

% can only measure the first state
C = [ 1 0 ];

% Initial guesses..
A_m = 10;		% [m^2]
A_in = 5;		% [m^2]
R_in_out = 5;	% [K/kW]
R_in_m = 2.5;	% [K/kW]
R_m_out = 2.5;	% [K/kW]
C_in = 5;		% [kWh/K]
C_m = 5;		% [kWh/K]

% x0 = [ A_m ; A_in ; R_in_out ; R_in_m ; R_m_out ; C_in ; C_m ];

% these are ones I found from the optimization routine
% a pretty good one found from the genetic alg
x0 = [3.02900539760059
3.95080553134958
7.99693842055921
0.380910424508521
4.47607200917178
9.23811845095962
14.8150486391551
-4.01687858868181
21.4568782170021];
% last two entries are vbar and guess of initial building mass temperature

% get the A and B matrices in continuous time
[Ac,Bc,Wc] = getSSDChouse(x0(1:7));

% assign initial condition
T0 = [ Y(1) x0(9) ];

% can change the options for the ode solver
% ODEoptions = odeset('RelTol',1e-6,'AbsTol',1e-10);
tsim = tdata;
tic
[t, Tsim] = ode15s(@(t,T) contSimDC(t,T,U,x0(8),tsim,Ac,Bc,Wc), tsim, T0);%,ODEoptions);
toc

figure; hold on;
plot(tdata,Tsim(:,1),'b-')
plot(tdata,Y,'b-.')

figure; hold on;
plot(tdata,Tsim(:,1),'b-',tdata,Tsim(:,2),'r-')

%% use optimization to fit the parameters

% bounds on the parameters
lb = zeros(size(x0,1),1)+1e-12;	% physical parameters must be positive
lb(8) = -10;	% disturbance offset could be large and negative
lb(9) = 15;		% initial building mass temperature lb
% ub = [ 20 ; 20 ; 15 ; 15 ; 15 ; 20 ; 20 ];
ub = [ 20 ; 10 ; 10 ; 5 ; 5 ; 15 ; 15 ];
ub(8) = 10;		% disturbance offset could be large and positive
ub(9) = 25;		% initial building mass temperature ub

% potential bounds
% Area for solar radation
% 0 < A_m < 20 % m^2
% 0 < A_in < 10 % m^2

% Resistances
% 0 < R_in_out < 10 K/kW
% 0 < R_in_m < 5    K/kW
% 0 < R_m_out < 5   K/kW

% Capacitances
% 0 < C_in < 10     kWh/C
% 0 < C_m < 10      kWh/C

%% try a genetic algorithm.. this seems to work a lot better

% runs a simulation for each of the candidates and calculates RMSE against
% the measurements Y
options = optimoptions('ga','PlotFcn', {@gaplotdistance,@gaplotrange,@gaplotbestf},'UseParallel',true,'MutationFcn','mutationadaptfeasible');
x = ga(@(x)getFitCost(x,Y,U,tdata),size(x0,1),[],[],[],[],lb,ub,[],options);
save('IDparameters_GA_700_1000','x')

%% try a pso
% options = optimoptions('particleswarm','SwarmSize',10,'PlotFcn', 'pswplotbestf','UseParallel',true);
% [x,fval,exitflag,output]  = particleswarm(@(x)getFitCost(x,Y,U,tdata),size(x0,1), lb,ub,options)


%% simulate the new identified system

% find new A and B matrices
[Ac,Bc,Wc] = getSSDChouse(x(1:7));

% assume they start in equilibrium.. from the data
T0 = [ Y(1) x(9) ];
% simulate the candidate model
[~, ySim] = ode15s(@(t,T) contSimDC(t,T,U,x(8),tsim,Ac,Bc,Wc), tsim, T0);%,ODEoptions);
figure; hold on;
plot(tdata,ySim(:,1),'b-','LineWidth',2)
plot(tdata,ySim(:,2),'r-','LineWidth',0.5)
% plot(tdata,ySim(:,1),'b-')
plot(tdata,Y,'k-.','LineWidth',2)
ylabel('Temperature')
title('3R2C Fit GA RMSE = 1.4 no restart')
legend('T_{air} RC','T_{m} RC','T_{air} exp.')

% a little bit better..
params0 = x';

% plot the error
error = Y-ySim(:,1);
figure; hold on;
plot(tdata,error)
a = 1;
B = 1;
%% now try a pattern search starting from this point..

options = optimoptions('patternsearch','UseParallel',true);
x_PS = patternsearch(@(x)getFitCost(x,Y,U,tdata),params0,[],[],[],[],lb,ub,[],options);
% % just stays at initial local minimum potentially global, can change to
% multi search if needed

%% Use a Kalman filter to identify the disturbance based on these model parameters

% define the covariance matrices.. these are all tuning parameters that will change how the estimator responds
% R = cov(Y)^2;
R = 0.1;	% how much do we trust the measurement?
Qx = eye(2)/100;	% how much do we trust the model?
Qv = eye(1)*2;	% assuming the disturbance is v = vbar + vhat.. what is variance in vhat from step to step if it is driven by gaussian noise
Q = blkdiag(Qx,Qv);

x = x0;
% find discrete time model
% find new A and B matrices
[Ac,Bc,Wc] = getSSDChouse(x(1:7));
sysc = ss(Ac,[Bc,Wc],eye(2),0);
sysd = c2d(sysc,Ts);
Ad = sysd.A;
Bd = sysd.B(:,1:3);
Wd = sysd.B(:,4);

% initial error covaraince.. how much do we trust these initial values?
Pp1 = blkdiag(100*eye(2),1000);

% and for augmented system with vhat as a state - X = [ x ; vhat ]
Akfd = [Ad Wd ; zeros(1,2) 1];
Bkfd = [ Bd ; zeros(1,3) ];
Wkfd = [ Wd ; 0 ];
Ckfd = [ 1 0 0 ];

% run the kalman filter over the data set
xp1 = [Y(1);x(9);0];	% assume vhat is vero at the start

% vectors for memory
save_ekf_state(1,:) = xp1(1:2);
save_ekf_dist(1,:) = xp1(3);

tic
for i = 1:numel(tdata)-1
	
	[xp1 , Pp1] = kf_d(Y(i+1),Pp1,xp1,U(i,:)',Akfd,Bkfd,Wkfd*x(8),Ckfd,Q,R);
	save_ekf_state(i+1,:) = xp1(1:2);
	save_ekf_dist(i+1,:) = xp1(3);

end
toc

figure; hold on;
plot(tdata,save_ekf_state(:,1),'b-',tdata,save_ekf_state(:,2),'r-')
plot(tdata,Y,'b-.')

figure; hold on;
plot(tdata,save_ekf_dist,'k-')

%% try simulation with the observed disturbances

% find new A and B matrices
[Ac,Bc,Wc] = getSSDChouse(x(1:7));

% assume they start in equilibrium.. from the data
T0 = [ Y(1) x(9) ];
% simulate the candidate model
[~, ySim] = ode15s(@(t,T) contSimDC(t,T,U,x(8)+save_ekf_dist,tsim,Ac,Bc,Wc), tsim, T0);%,ODEoptions);

figure; hold on;
plot(tdata,ySim(:,1),'b-',tdata,ySim(:,2),'r-')
% plot(tdata,ySim(:,1),'b-')
plot(tdata,Y,'b-.')


% Code for computing heat supply prediction of model
% w_dis = Wd*(-4.5)
% dis = [w_dis(1) *  ones(length(Y_heat_calc(1:end-1,:)),1) w_dis(2)*ones(length(Y_heat_calc(1:end-1,:)),1)] ;
% Y_heat_calc = [Y ySim(:,2)];
% Q = pinv(Bd)*(Y_heat_calc(2:end,:) - (Ad*Y_heat_calc(1:end-1,:)')' + dis )'
% 
% figure(66)
% plot(Q(2,:)), hold on
% plot(U(1:end,2))
% 
% rmse(U(1:end-1,2)',Q(2,:))