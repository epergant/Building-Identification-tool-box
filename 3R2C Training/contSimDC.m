% continuous time simulation file for DC house model
%
%	To be called by ode solver
% 
% [~, ySim] = ode15s(@(t,T) contSimDC(t,T,U,v,tsim,Ac,Bc,Wc), tsim, T0);
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function dT = contSimDC(t,x,U,D,td,Ac,Bc,Wc)

if size(U,1) > 1
	u(1,1) = interp1(td,U(:,1),t,'linear');
	u(2,1) = interp1(td,U(:,2),t,'linear');
	u(3,1) = interp1(td,U(:,3),t,'linear');
else
	u = U;
end

if size(D,1) > 1
	d = interp1(td,D,t,'linear');
else
	d = D;
end

dT = Ac*x + Bc*u + Wc*d;

end