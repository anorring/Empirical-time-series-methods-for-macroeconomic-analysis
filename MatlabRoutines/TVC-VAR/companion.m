% Construct the companion representation of a state space model
% with state vector theta p lags in the VAR representation and 
% n variables and c=1 constant c=0 no constant.

function C = companion(theta,p,n,c)
if c == 1
    np = (n*p+1);
elseif c == 0
    np = (n*p);
end   
theta = theta';
for i = 1:n
    M(i,1:np) = theta((i-1)*np+1:i*np);
end
if c == 1;
    M = M(:,2:end);
end
C = [M ; eye(n*(p-1)) zeros(n*(p-1),n)];
