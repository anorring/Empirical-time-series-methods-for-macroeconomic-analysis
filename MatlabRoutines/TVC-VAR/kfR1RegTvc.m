function [S0,P0,P1] = kfR1RegTvc(Y,X,Q,R,SI,PI,T);

% function [S0,P0,P1] = KFR(Y,X,Q,B,H,SI,PI,T,N,L);

% This file performs the forward kalman filter recursions
% for the random coefficients VAR.  R is time varying, depening
% on the stochastic volatilities: R = inv(B)*H(t)*B

% C is zero

% SI,PI are the initial values for the recursions, S(1|0) and P(1|0)

P1(:,:,1) = PI; % P(1|0)
K = (P1(:,:,1)*X(:,1))*inv( X(:,1)'*P1(:,:,1)*X(:,1) + R(1)); % K(1)
P0(:,:,1) = P1(:,:,1) - K*(X(:,1)'*P1(:,:,1)); % P(1|1)
S0(:,1) = SI + K*( Y(:,1) - X(:,1)'*SI ); % S(1|1)

% Iterating through the rest of the sample
for i = 2:T,
   P1(:,:,i) = P0(:,:,i-1) + Q; % P(t|t-1)
   K = (P1(:,:,i)*X(:,i))*inv( X(:,i)'*P1(:,:,i)*X(:,i) + R(i)); % K(t)
   P0(:,:,i) = P1(:,:,i) - K*(X(:,i)'*P1(:,:,i)); % P(t|t)
   S0(:,i) = S0(:,i-1) + K*( Y(:,i) - X(:,i)'*S0(:,i-1) ); % S(t|t)
end
