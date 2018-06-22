function [TQ,DF] = IWPQ(theta,N,T,L,TQ0,df);

%function [TQ,DF] = IWPQ(theta,N,T,L,TQ0,df);

% This file computes posterior estimates of TQ
% for an informative prior, in which Q is inverse
% wishart with degrees of freedom df and scale matrix
% TQ0.  The posterior density is also inverse 
% wishart, with scale matrix TQ and degrees
% of freedom DF

% N = number of equations
% T = number of time periods
% L = number of lags

v(:,2:T) = theta(:,2:T) - theta(:,1:T-1);
TQ = TQ0 + v*v';
DF = df + T;