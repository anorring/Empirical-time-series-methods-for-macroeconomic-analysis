function SA = GIBBS1(Y,X,S0,P0,P1,T,N,L,c);

SA = zeros(N*(c+N*L),T); % artificial states
SM = zeros(N*(c+N*L),1); % backward update for conditional mean of state vector
PM = zeros(N*(c+N*L),N*(c+N*L)); % backward update for projection matrix
P = PM; % backward update for covariance matrix
wa = randn(N*(c+N*L),T); % draws for state innovations

% Backward recursions and sampling
% Terminal state
SA(:,T) = S0(:,T) + real(sqrtm(P0(:,:,T)))*wa(:,T);

% iterating back through the rest of the sample
for i = 1:T-1,
   PM = P0(:,:,T-i)*inv(P1(:,:,T-i+1));
   P = P0(:,:,T-i) - PM*P0(:,:,T-i);
   SM = S0(:,T-i) + PM*(SA(:,T-i+1) - S0(:,T-i));
%    SA(:,T-i) = SM + real(sqrtm(P))*wa(:,T-i);
   SA(:,T-i) = SM + chol(P)'*wa(:,T-i);
end
   
   