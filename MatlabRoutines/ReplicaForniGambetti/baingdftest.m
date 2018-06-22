function [stat, e] = baingdftest(X,rhat,pmax)

[T, N] = size(X);
Fhat = principalcomponents(X,rhat);
[beta,e] = myvar(Fhat,pmax);
xepsilon = [1 1];
eps = xepsilon/min([N^(.4);T^(.4)]);
[stat1,stat2,evals_d]=dfactest(e,N,T,eps,1);
xepsilon=[1.25 2.25];
eps = xepsilon/min([N^(.4);T^(.4)]);
[stat3,stat4,evals_d]=dfactest(e,N,T,eps,2);
stat = [stat1 stat2; stat3 stat4];

