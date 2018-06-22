function qhat = swdftest(X,rhat,nlagsvar)
[T, N] = size(X);
[Fhat, R] = principalcomponents(X,rhat);
[beta,e] = myvar(Fhat,nlagsvar);
eta = standardize(X(nlagsvar+1:end,:) - X(nlagsvar+1:end,:)*R*R' + e*R');
[PC, IC] = baingcriterion(eta, rhat);
[minimum, qhat] = min(IC(:,1:2));

