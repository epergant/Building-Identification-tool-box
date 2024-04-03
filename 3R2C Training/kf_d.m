% kalman filter to observe the disturbance heat load occuring DC house
% measurements
%
%	[xp1 , Pp1] = kf_d(y,P0,x0,u0,A,B,c,C,Q,R)
%
%	Inputs:
%		y : [1x1] measurement at time k+1
%		P0 : [3x3] covariance of estimation error at time k
%		x0 : [3x1] estimate of state and disturbance at time k
%		u0 : [3x1] inputs applied at time k
%		A : [3x3] discrete time dynamics for [x;d]
%		B : [3x3] discrete time input matrix for u
%		c : [3x1] discrete time offset from constant disturbance
%		C : [1x3] output matrix s.t. y = Cx
%		Q : [3x3] process noise variance
%		R : [1x1] measurement variance
% 
%	Outputs:
%		xp1 : [3x1] update estimate at time k + 1
%		Pp1 : [3x3] estimated error covariance at time k + 1
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [xp1 , Pp1] = kf_d(y,P0,x0,u0,A,B,c,C,Q,R)

% update prediction
xm = A*x0 + B*u0 + c;
Pm = A*P0*A' + Q;

% Kalman gain
L = Pm*C'/(C*Pm*C' + R);

% measurement update
xp1 = xm + L*(y-C*xm);
Pp1 = Pm - L*C*Pm;

end