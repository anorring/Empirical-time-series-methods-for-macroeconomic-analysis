function [PC, IC] = baingcriterion(X, rmax, rvar)
if nargin == 2
    rvar = rmax;
end
N = size(X, 2);
T = size(X, 1);
W = diag(std(X));
x = center(X)*(W^-1);

Gamma0 = cov(x);
opts.disp = 0;
[H, D] = eigs(Gamma0, rmax,'LM',opts);

for r = 1:rmax
I = x*(eye(N) - H(:,1:r)*H(:,1:r)');
V(r) = sum(sum(I.^2))/(N*T);
end
e = (N + T)/(N*T);

penalty1 = (1:rmax)*e*log(1/e);
penalty2 = (1:rmax)*e*log(min(N,T));
penalty3 = (1:rmax)*(log(min(N,T))/min(N,T));

PC = [V + V(rvar)*penalty1; V + V(rvar)*penalty2; V + V(rvar)*penalty3]';
IC = [log(V) + penalty1; log(V) + penalty2; log(V) + penalty3]';
