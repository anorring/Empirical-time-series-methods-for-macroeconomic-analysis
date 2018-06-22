function SA = GIBBS1Reg(Y,X,S0,P0,P1,T,NP);

% function SA = GIBBS1(Y,X,S0,P0,P1,T,N,L);

% This file performs the first step of the Gibbs sampler,
% drawing from p(theta | Q,R,Y) 

wa = randn(NP,T); % draws for state innovations

% Backward recursions and sampling
% Terminal state
%SA(:,T) = S0(:,T) + real(sqrtm(P0(:,:,T)))*wa(:,T);
SA(:,T) = S0(:,T) + chol(P0(:,:,T))'*wa(:,T);

% iterating back through the rest of the sample
for i = 1:T-1,
   PM = P0(:,:,T-i)*inv(P1(:,:,T-i+1));
   P = P0(:,:,T-i) - PM*P0(:,:,T-i);
   SM = S0(:,T-i) + PM*(SA(:,T-i+1) - S0(:,T-i));
   %SA(:,T-i) = SM + real(sqrtm(P))*wa(:,T-i);
   SA(:,T-i) = SM + chol(P)'*wa(:,T-i);
end
   
% Updating SI,PI
%PM = P0(:,:,1)*inv(P1(:,:,1));
%PI = P0(:,:,1) - PM*P0(:,:,1);
%SI = S0(:,1) + PM*(SA(:,1) - S0(:,1));
   
   