function [y, M, s,MM,ss] = standardize(x);
s = std(x);
M = mean(x);
ss = ones(size(x , 1) , 1)*s;
MM = ones(size(x , 1) , 1)*M;
y = (x - MM)./ss;