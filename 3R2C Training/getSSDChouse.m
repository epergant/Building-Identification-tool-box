% returns the state space matrices for the DC house model
%
%	Inputs:
%		x : [7x1] list of paramters
%	Outputs:
%		Dynamic matrics s.t. xdot = Ax + Bu + Wv
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


function [ A , B , W ] = getSSDChouse(x)

% extract parameters from x
A_m = x(1);
A_in = x(2);
R_in_out = x(3);
R_in_m = x(4);
R_m_out = x(5);
C_in = x(6);
C_m = x(7);

% build matrices such that xdot = Ax + Bu + Wv
A(1,1) = -1/(R_in_out * C_in) - 1/(R_in_m * C_in);
A(1,2) = 1/(R_in_m * C_in);
A(2,1) = 1/(R_in_m * C_m);
A(2,2) = - 1/(R_in_m * C_m) - 1/(R_m_out * C_m);

B(1,1) = 1/(R_in_out * C_in);
B(1,2) = 1/C_in;
B(1,3) = A_in/C_in;
B(2,1) = 1/(R_m_out * C_m);
B(2,2) = 0;
B(2,3) = A_m/C_m;

W(1,1) = (1/C_in);
W(2,1) = (1/C_m);

end