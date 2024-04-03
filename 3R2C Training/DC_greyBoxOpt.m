% this function calls fmincon to try and fit the model with parameters x to
% the data Y for inputs U
%
% iteratively calls getFitCost(x,Y,U,tdata) to run simulation and calculate
% RMSE
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function x = DC_greyBoxOpt(x0,ub,lb,Y,U,tdata)

% I think it is fun to watch.. lol
fig = figure; hold on;
h2 = plot(tdata,Y,'b-');
h1 = plot(tdata,Y,'b-');

% for nested objective and constraint functions
xLast = 0;	% remember the parameters for the last simulation call
yLast = 0;	% remember the simulation results for the last call
JLast = 0;	% remember the cost function for the last call

% UseParallel
% opts = optimoptions('fmincon','Algorithm','sqp','UseParallel',true,'Display','iter','PlotFcn','optimplotfval');
opts = optimoptions('fmincon','Algorithm','interior-point','UseParallel',false,'Display','iter','PlotFcn','optimplotfval','StepTolerance',1e-14);

problem = createOptimProblem('fmincon','objective',...
    @(x) cost(x,Y,U,tdata),'x0',x0,'lb',lb,'ub',ub,'options',opts);
ms = MultiStart;
[x,f] = run(ms,problem,5);
	
	% cost function
	function J = cost(x,Y,U,tdata)
	
		% check if we need to resimulate the system
		if ~isequal(x,xLast)
			[J,ySim] = getFitCost(x,Y,U,tdata);
			xLast = x;
			yLast = ySim;
			JLast = J;
		else
			ySim = yLast;
			J = JLast;
		end
		
		if numel(ySim) < numel(Y)
			% sometimes the simulation fails.. this value is way above the initial guess x0
			J = 1000;
		else
			figure(fig)
			delete(h1)
			h1 = plot(tdata,ySim(:,1),'r--');
		end
		
	end

end

