function [J,ySim] = getFitCost(x,Y,U,tdata)

[Ac,Bc,Wc] = getSSDChouse(x(1:7));
		
% assign initial condition
T0 = [ Y(1) x(9) ];
% simulate the candidate model
[~, ySim] = ode15s(@(t,T) contSimDC(t,T,U,x(8),tdata,Ac,Bc,Wc), tdata, T0);%,ODEoptions);

if size(ySim,1) ~= size(Y,1)
	stop = 1;
end

% minimize the root mean square error
err = sqrt(sum((ySim(:,1) - Y).^2)/size(Y,1));
J = err;

end