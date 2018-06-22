function u = innovm(Y,X,theta,N,T,L,c);

% function u = innovm(Y,X,theta,N,T,L);

% this file computes recursive residuals for an unrestricted VAR
% N = number of equations
% T = number of time periods
% L = number of lags
% Y = matrix of dependent variables; NxT
% X = matrix of right hand variables; identical in all equations; Txk
% theta = matrix of time-varying coefficients; N*(1+NL) x T

for i = 1:N,
   u(i,1:T) = Y(i,1:T)- sum((X(1:T,:).*theta((i-1)*(N*L+c)+1:i*(N*L+c),1:T)')');
end
