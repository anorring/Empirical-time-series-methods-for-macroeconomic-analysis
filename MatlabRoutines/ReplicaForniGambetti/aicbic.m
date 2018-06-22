function [aic, bic, e] = aicbic(y , pmax)
[T,N] = size(y);
for p = 1:pmax
[beta, e] = myvar(y, p);
a = log(det((e'*e)/T));
b = (p*N^2 + N)/T;
aic(p) = a + 2*b;
bic(p) = a + log(T)*b;
end
[minaic,aic] = min(aic);
[minbic,bic] = min(bic);