function XC = center(X)
[T n] = size(X);
XC = X - ones(T,1)*(sum(X)/T); 