% Authors Elias N. Pergantis, Kevin J. Kircher. If you use please cite:

% E. N. Pergantis, Priyadarshan, N. A. Theeb, P. Dhillon,
% J. P. Ore, D. Ziviani, E. A. Groll, K. J. Kircher, Field
% demonstration of predictive heating control for an all-
% electric house in a cold climate, Applied Energy 360
% (2024) 122820

close all
clear all
clc

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
%% introduction
% This script fits a 2R1C thermal circuit model of the DC house.
% set(groot, 'defaultAxesTickLabelInterpreter','none'); set(groot, 'defaultLegendInterpreter','none');

%% input data
% DC house data import
fileName = 'DC_House_Nov_April.xlsx';
opts = detectImportOptions(fileName);
opts.PreserveVariableNames = 1;
rawData = readtable(fileName,opts);
rawData = rmmissing(rawData); % Remove missing data, experiment with backfill or mean
timeStamp = rawData{:,1};
T = rawData{:,2}; % indoor temperature, C
Qdotc = rawData{:,6}; % heat pump thermal power, kW
Tout = rawData{:,3}; % outdoor temperature, C

% discard the non-heating data
iKeep = timeStamp >= datetime(2022,11,11,0,0,0); % indices of heating data
T = T(iKeep);
Qdotc = Qdotc(iKeep);
Tout = Tout(iKeep);
timeStamp = timeStamp(iKeep);

% retime the data
dt = 1; % time step, h
% TT = retime(timetable(timeStamp,T,Tout,Qdotc),'regular','timestep',hours(dt));

% timing
K = length(timeStamp)-1; % number of time steps
t = (0:dt:K*dt)'; % time span as array, h

%% exploratory data analysis
% plot parameters
lw = 3; % line width
fs = 12; % font size

% load vs. outdoor temperature
figure(1), clf
subplot(4,2,[2 4 6 8]), plot(T-Tout,Qdotc,'ko'), grid on, axis square
xlabel('Indoor-outdoor temperature difference ($^\circ$C)','fontsize',fs)
ylabel('Heat pump thermal power (kW)','fontsize',fs)
ylim([0 6])

% add trend line
p = polyfit(T-Tout,Qdotc,1);
px = [min(T-Tout) max(T-Tout)];
py = polyval(p,px);
hold on, plot(px,py,'m--','linewidth',lw);
text(mean(xlim),max(ylim),sprintf('$y = %.3g x %.3g$',p(1),p(2)),...
    'color','m','fontsize',fs,'verticalalignment','top','horizontalalignment','center')

% indoor temperature vs. time
tLim = [timeStamp(1) timeStamp(end)];
TLim = 21 + 3*[-1 1];
subplot(4,2,1), plot(timeStamp,T,'k','linewidth',lw), grid on
xlim(tLim)
ylim(TLim)
ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)

% outdoor temperature vs. time
subplot(4,2,3), plot(timeStamp,Tout,'k','linewidth',lw), grid on
xlim(tLim)
ylabel({'Outdoor','temperature ($^\circ$C)'},'fontsize',fs)

% load vs. time
subplot(4,2,5), plot(timeStamp,Qdotc,'k','linewidth',lw), grid on
xlim(tLim)
ylabel({'Heat pump','thermal','power (kW)'},'fontsize',fs)
sgtitle('Exploratory data analysis','fontsize',fs)

%% thermal resistance fit
% extract thermal resistance training and validation data
% (here this means night with steady indoor temperature)
tStart = datetime(2022,11,15,20,0,0); % start of unsteady period
tEnd = datetime(2023,1,1,16,0,0); % end of unsteady period
iFit = find(... % indices of steady night hours
    (timeStamp <= tStart | timeStamp >= tEnd) & ...
    (hour(timeStamp) <= 6 | hour(timeStamp) >= 20));

% define training and validation data
it = iFit(1:round(0.75*length(iFit))); % training indices
iv = iFit(length(it)+1:end); % validation indices
Xt = [ones(length(it),1), Qdotc(it)]; % training feature, kW
yt = T(it) - Tout(it); % training target, C
Xv = [ones(length(iv),1), Qdotc(iv)]; % validation feature, kW
yv = T(iv) - Tout(iv); % validation target, C

% joint fit of thermal resistance and exogenous thermal power
cvx_begin quiet
    variable betaHat(2,1)
    minimize(sum_square(yt - Xt*betaHat))
    subject to
        0 <= betaHat(1) <= 9 % temperature disturbance constraints, C
        0 <= betaHat(2) <= 5 % thermal resistance constraints, C/kW
cvx_end
% betaHat = (Xt'*Xt)\(Xt'*yt); % least-squares fit
RoutHat = betaHat(2); % estimated thermal resistance, C/kW
QdoteHat = betaHat(1)/RoutHat; % estimated exogenous thermal power, C/kW

% report fit statistics
ytHat = Xt*betaHat; % training prediction, C
yvHat = Xv*betaHat; % validation prediction, C
et = ytHat - yt; % training prediction error, C
ev = yvHat - yv; % validation prediction error, C
fprintf('---------------------------------------------------------------\n')
fprintf('Thermal resistance training RMSE: %.3g C. ',sqrt(mean(et.^2)))
fprintf('Validation: %.3g C.\n',sqrt(mean(ev.^2)))

% use real thermal load to predict indoor temperature
TtHat = Tout(it) + RoutHat*(Qdotc(it) + QdoteHat);
TvHat = Tout(iv) + RoutHat*(Qdotc(iv) + QdoteHat);

% use real indoor temperature to predict thermal load
QdotctHat = (T(it) - Tout(it))/RoutHat - QdoteHat;
QdotcvHat = (T(iv) - Tout(iv))/RoutHat - QdoteHat;

% plot temperature fit
QdotLim = [0 5.5];
figure(2), clf
subplot(2,2,1), plot(timeStamp(it),T(it),'ko',timeStamp(it),TtHat,'mo','linewidth',lw), grid on
ylim(TLim), ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Training','fontsize',fs)
subplot(2,2,2), plot(timeStamp(iv),T(iv),'ko',timeStamp(iv),TvHat,'mo','linewidth',lw), grid on
ylim(TLim), ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Validation','fontsize',fs)

% plot thermal load fit
subplot(2,2,3), plot(timeStamp(it),Qdotc(it),'ko',timeStamp(it),QdotctHat,'mo','linewidth',lw), grid on
ylim(QdotLim), ylabel({'Thermal','load (kW)'},'fontsize',fs)
legend('Truth','Estimate','fontsize',fs,'location','south','orientation','horizontal')
subplot(2,2,4), plot(timeStamp(iv),Qdotc(iv),'ko',timeStamp(iv),QdotcvHat,'mo','linewidth',lw), grid on
ylim(QdotLim), ylabel({'Thermal','load (kW)'},'fontsize',fs)
sgtitle('Thermal resistance fit','fontsize',fs)

% refit using all data
cvx_begin quiet
    variable betaHat(2,1)
    minimize(sum_square([yt;yv] - [Xt;Xv]*betaHat))
    subject to
        0 <= betaHat(1) <= 9 % temperature disturbance constraints, C
        0 <= betaHat(2) <= 5 % thermal resistance constraints, C/kW
cvx_end
RoutHat = betaHat(2); % estimated thermal resistance, C/kW
% QdoteHat = betaHat(1)/RoutHat; % estimated exogenous thermal power, C/kW

%% thermal capacitance fit
% define thermal mass temperature and air-mass resistance
Tm0 = mean(T(timeStamp >= datetime(2022,11,11,19,0,0) & ...
    timeStamp <= datetime(2022,11,15,10,0,0))); % thermal mass temperature
Tm = Tm0*ones(K+1,1);

% add mass temperature to input data plot
figure(1), subplot(4,2,1), hold on
plot(timeStamp,Tm,'m--','linewidth',lw)
legend('Air','Mass','fontsize',fs,'location','north','orientation','horizontal')

% define indoor air-mass thermal resistance values
n = 100; % number of Rm values
RmValues = linspace(0.01,10,n); % indoor air-mass thermal resistance, C/kW

% extract data with indoor temperature excitation;
it = find(timeStamp >= datetime(2022,11,25,19,0,0) & ...  % second overnight step response
    timeStamp <= datetime(2022,12,3,2,0,0));
iv = find(timeStamp >= datetime(2022,12,3,2,0,0) & ... % first overnight step response
    timeStamp <= datetime(2022,12,5,2,0,0));

% data storage
aHatValues = zeros(n,1); % discrete-time dynamics parameter
QdoteHatValues = zeros(n,1); % exogenous thermal power, kW
tRMSE = zeros(n,1); % training RMSE, C
vRMSE = zeros(n,1); % validation RMSE, C

% loop over Rm values
for i=1:n
    % define air-mass thermal resistance
    Rm = RmValues(i);
    
    % define effective temperature and resistance
    RHat = RoutHat*Rm/(RoutHat + Rm);
    theta = (Rm*Tout + RoutHat*Tm)/(RoutHat + Rm);

    % define training and validation data
    Xt = [ones(length(it),1), T(it) - theta(it) - RHat*Qdotc(it)]; % training feature, C
    yt = T(it+1) - theta(it) - RHat*Qdotc(it); % training target, C
    Xv = [ones(length(iv),1), T(iv) - theta(iv) - RHat*Qdotc(iv)]; % validation feature, C
    yv = T(iv+1) - theta(iv) - RHat*Qdotc(iv); % validation target, C

    % joint fit of thermal resistance and exogenous thermal power
    betaHat = (Xt'*Xt)\(Xt'*yt); % least-squares fit
    aHatValues(i) = betaHat(2); % estimated discrete-time dynamics parameter
    QdoteHatValues(i) = betaHat(1)/(RHat*(1-betaHat(2))); % estimated exogenous thermal power, C/kW

    % compute training and validation RMSEs
    ytHat = Xt*betaHat; % training prediction, C
    yvHat = Xv*betaHat; % validation prediction, C
    et = ytHat - yt; % training prediction error, C
    ev = yvHat - yv; % validation prediction error, C
    tRMSE(i) = sqrt(mean(et.^2)); % training RMSE, C
    vRMSE(i) = sqrt(mean(ev.^2)); % validation RMSE, C
end

% set the weights to balance training and validation RMSE in cost function
lambda = 0.2; % weight on normalized training RMSE
tWeight = lambda/mean(tRMSE);
vWeight = (1-lambda)/mean(vRMSE);

% plot the RMSEs and cost function
figure(3), clf
subplot(3,1,1), plot(RmValues,tRMSE/mean(tRMSE),'k','linewidth',lw), grid on
ylabel({'Normalized','training','RMSE ($^\circ$C)'},'fontsize',fs)
title('Air-mass resistance tuning','fontsize',fs)
subplot(3,1,2), plot(RmValues,vRMSE/mean(vRMSE),'m','linewidth',lw), grid on
ylabel({'Normalized','validation','RMSE ($^\circ$C)'},'fontsize',fs)
subplot(3,1,3), plot(RmValues,tWeight*tRMSE + vWeight*vRMSE,'c','linewidth',lw), grid on
ylabel({'Cost','function','($^\circ$C)'},'fontsize',fs)
xlabel('$R_m$ ($^\circ$C/kW)','fontsize',fs)

% pick an Rm value that minimizes a mixture of training and validation RMSE
[~,iStar] = min(tWeight*tRMSE + vWeight*vRMSE);
RmHat = RmValues(iStar); % C/kW

% redefine effective temperature and resistance
RHat = RoutHat*RmHat/(RoutHat + RmHat);
theta = (RmHat*Tout + RoutHat*Tm)/(RoutHat + RmHat);

% redefine training and validation data
Xt = [ones(length(it),1), T(it) - theta(it) - RHat*Qdotc(it)]; % training feature, C
yt = T(it+1) - theta(it) - RHat*Qdotc(it); % training target, C
Xv = [ones(length(iv),1), T(iv) - theta(iv) - RHat*Qdotc(iv)]; % validation feature, C
yv = T(iv+1) - theta(iv) - RHat*Qdotc(iv); % validation target, C

% joint fit of thermal resistance and exogenous thermal power
betaHat = (Xt'*Xt)\(Xt'*yt); % least-squares fit
aHat = betaHat(2); % estimated discrete-time dynamics parameter
QdoteHat = betaHat(1)/(RHat*(1-betaHat(2))); % estimated exogenous thermal power, C/kW

% extract associated parameter estimates
% aHat = aHatValues(iStar);
% QdoteHat = QdoteHatValues(iStar); % kW
% betaHat = [QdoteHat*RHat*(1-aHat);aHat];

% report fit statistics
ytHat = Xt*betaHat; % training prediction, C
yvHat = Xv*betaHat; % validation prediction, C
et = ytHat - yt; % training prediction error, C
ev = yvHat - yv; % validation prediction error, C
fprintf('Discrete-time dynamics parameter training RMSE: %.3g C. ',sqrt(mean(et.^2)))
fprintf('Validation: %.3g C.\n',sqrt(mean(ev.^2)))

% use real thermal load to predict indoor temperature
TtHat = aHat*T(it) + (1-aHat)*(theta(it) + RHat*(Qdotc(it) + QdoteHat)); % one-step training temperature predictions
TvHat = aHat*T(iv) + (1-aHat)*(theta(iv) + RHat*(Qdotc(iv) + QdoteHat)); % one-step validation temperature predictions

% use real indoor temperature to predict thermal load
QdotctHat = ((T(it+1) - aHat*T(it))/(1-aHat) - theta(it))/RHat - QdoteHat; % one-step training load predictions
QdotcvHat = ((T(iv+1) - aHat*T(iv))/(1-aHat) - theta(iv))/RHat - QdoteHat; % one-step validation load predictions

% plot temperature fit
QdotLim = ceil(max(abs([QdotctHat;QdotcvHat]))/5)*5*[-1 1];
figure(4), clf
subplot(2,2,1), plot(timeStamp(it),T(it),'ko',timeStamp(it),TtHat,'mo','linewidth',lw), grid on
ylim(TLim), ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Training','fontsize',fs)
subplot(2,2,2), plot(timeStamp(iv),T(iv),'ko',timeStamp(iv),TvHat,'mo','linewidth',lw), grid on
ylim(TLim), ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Validation','fontsize',fs)

% plot thermal load fit
subplot(2,2,3), plot(timeStamp(it),Qdotc(it),'ko',timeStamp(it),QdotctHat,'mo','linewidth',lw), grid on
ylim(QdotLim), ylabel({'Thermal','load (kW)'},'fontsize',fs)
legend('Truth','Estimate','fontsize',fs,'location','south','orientation','horizontal')
subplot(2,2,4), plot(timeStamp(iv),Qdotc(iv),'ko',timeStamp(iv),QdotcvHat,'mo','linewidth',lw), grid on
ylim(QdotLim), ylabel({'Thermal','load (kW)'},'fontsize',fs)
sgtitle('Thermal capacitance fit','fontsize',fs)

% retrain over the full excitation data set
y = [yt;yv];
X = [Xt;Xv];
betaHat = (X'*X)\(X'*y); % least-squares fit
aHat = betaHat(2); % estimated discrete-time dynamics parameter
% QdoteHat = betaHat(1)/(RHat*(1-betaHat(2))); % kW

% report estimates
fprintf('Exogenous thermal power estimate: QdoteHat = %.5g kW.\n',QdoteHat)
CHat = -dt/(RHat*log(aHat)); % air thermal capacitance estimate
fprintf('Air thermal capacitance estimate: CHat = %.5g kWh/C.\n',CHat)
cp = 2.792e-4; % specific heat of air at constant volume, kWh/kg/C
rho = 1.2754; % density of air, kg/m^3
V = 208*6; % enclosed air volume, m^3
fprintf('For the enclosed air, rho*cv*V ~= %.5g kWh/C.\n',rho*cp*V)

%% exogenous thermal power fit
% define training and validation targets
% (this is the exogenous thermal power consistent with measurements and
% parameter fits RHat and aHat)
Qdote = ((T(2:K+1) - aHat*T(1:K))/(1-aHat) - theta(1:K))/RHat - Qdotc(1:K); % target, kW

% add exogenous thermal power to input data plot
figure(1)
subplot(4,2,7), plot(timeStamp(1:K),Qdote,'k','linewidth',lw), grid on
xlim(tLim)
ylabel({'Exogenous','thermal','power (kW)'},'fontsize',fs)

% this is now a supervised learning (time series regression) problem
% .
% .
% .
% insert model (ANN, SVM, SARIMAX, ...) training and validation here
% .
% .
% .
% what follows is just a simple autoregressive model

% define candidate memories
mValues = 0:24/dt;
nm = length(mValues);

% data storage
tRMSE = zeros(nm,1); % training RMSEs
vRMSE = zeros(nm,1); % validation RMSEs

% loop over memories
for i=1:nm
    % generate target and features
    m = mValues(i); % memory
    y = Qdote(m+1:K); % target vector
    X = zeros(K-m,m+1); % feature matrix
    X(:,1) = 1; % constant feature
    for j=1:m
        X(:,j+1) = Qdote(j:K-m+j-1);
    end

    % separate training and validation data
    n = length(y); % total sample size
    nt = round(2*n/3); % training sample size
    nv = n - nt; % validation sample size
    it = nv+1:n; % training indices
    iv = 1:nv; % validation indices
    Xt = X(it,:); % training feature matrix
    Xv = X(iv,:); % validation feature matrix
    yt = y(it); % training target vector
    yv = y(iv); % training target vector

    % fit
    betaHat = (Xt'*Xt)\(Xt'*yt); % least-squares training
    
    % RMSEs
    ytHat = Xt*betaHat; % training prediction
    yvHat = Xv*betaHat; % validation prediction
    tRMSE(i) = sqrt(mean((yt-ytHat).^2));
    vRMSE(i) = sqrt(mean((yv-yvHat).^2));
end

% plot RMSE
figure(5), clf
subplot(2,1,1), plot(mValues,tRMSE,'ko-','linewidth',lw), grid on
ylabel('Training RMSE','fontsize',fs)
title('Exogenous thermal power model memory tuning','fontsize',fs)
subplot(2,1,2), plot(mValues,vRMSE,'mo-','linewidth',lw), grid on
ylabel('Validation RMSE','fontsize',fs)
xlabel('Number of time steps of memory','fontsize',fs)

% pick a decent m and retrain
m = 2;
y = Qdote(m+1:K); % target vector
X = zeros(K-m,m+1); % feature matrix
X(:,1) = 1; % constant feature
for j=1:m
    X(:,j+1) = Qdote(j:K-m+j-1);
end

% fit on the training data
n = length(y); % total sample size
nt = round(2*n/3); % training sample size
nv = n - nt; % validation sample size
it = nv+1:n; % training indices
Xt = X(it,:); % training feature matrix
yt = y(it); % target vector
betaHat = (Xt'*Xt)\(Xt'*yt); % least-squares training

% compute training and validation predictions
yHat = X*betaHat; % predictions
QdoteHat = [Qdote(1:m);yHat];

%% full model assessment
% predict temperature using real thermal load
THat = aHat*T(1:K) + (1-aHat)*(theta(1:K) + RHat*(Qdotc(1:K) + QdoteHat(1:K))); % one-step training temperature predictions

% predict thermal load using real temperature
QdotcHat = ((T(2:K+1) - aHat*T(1:K))/(1-aHat) - theta(1:K))/RHat - QdoteHat(1:K); % predicted heat pump thermal power (used in fit), kW

% plot the temperature fit
iPlot = 1:K; % indices to plot
figure(6), clf
subplot(2,1,1), plot(timeStamp(iPlot),T(iPlot+1),'ko-',timeStamp(iPlot),THat(iPlot),'mo-','linewidth',lw), grid on
ylabel('Indoor temperature ($^\circ$C)','fontsize',fs)
xline(timeStamp(nv),'b--','linewidth',lw)
legend('True','Predicted','fontsize',fs,'location','north','orientation','horizontal')
title('Full model fit','fontsize',fs)
text(timeStamp(nv),min(ylim),'$\leftarrow$ Validate \quad','fontsize',fs,'color','b',...
    'horizontalalignment','right','verticalalignment','bottom')
text(timeStamp(nv),min(ylim),'\quad Train $\rightarrow$','fontsize',fs,'color','b',...
    'horizontalalignment','left','verticalalignment','bottom')

% plot the thermal load fit
subplot(2,1,2), plot(timeStamp(iPlot),Qdotc(iPlot),'ko-',timeStamp(iPlot),QdotcHat(iPlot),'mo-','linewidth',lw), grid on
ylabel('Heat pump thermal power (kW)','fontsize',fs)
xline(timeStamp(nv),'b--','linewidth',lw)

%% report the final parameter values
% retrain disturbance model on full dataset
betaHat = (X'*X)\(X'*y); % least-squares training
QdoteBaseHat = betaHat(1); % baseline exogenous thermal power
alphaHat = betaHat(2:m+1); % autoregresion parameter vector

% report final results
fprintf('---------------------------------------------------------------\n')
fprintf('Indoor air-outdoor air resistance estimate: RoutHat = %.5g C/kW.\n',RoutHat)
fprintf('Indoor air-thermal mass resistance estimate: RmHat = %.5g C/kW.\n',RmHat)
fprintf('Effective 1R1C resistance estimate: RHat = %.5g C/kW.\n',RHat)
CHat = -dt/(RHat*log(aHat)); % air thermal capacitance estimate
fprintf('Indoor air capacitance estimate: CHat = %.5g kWh/C.\n',CHat)
fprintf('Discrete-time dynamics parameter estimate with dt = %.3g h: aHat = %.5g.\n',dt,aHat)
fprintf('Baseline exogenous thermal power estimate: QdoteBaseHat = %.5g kW.\n',QdoteBaseHat)
fprintf('Autoregressive model parameters with memory m = %i:\n',m)
for i=1:m
    fprintf('\talphaHat_%i = %.5g\n',i,alphaHat(i))
end











