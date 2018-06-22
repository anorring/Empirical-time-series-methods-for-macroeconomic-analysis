function [y, M, s,MM,ss] = standardize(x,w);
if nargin == 1
    w = ones(1,size(x,2));
end
s = std(x)./w;
M = mean(x);
ss = ones(size(x , 1) , 1)*s;
MM = ones(size(x , 1) , 1)*M;
y = (x - MM)./ss;