%Creates matrices for VAR(k)
function [Y,X] = VarStr(y,c,k)
[T N]=size(y);

Y = y(k+1:T,:);

X=[];

for j=1:k,
    X=[X y(k+1-j:T-j,:)];
end

if c == 1
    X=[ones(size(X,1),1) X];
end
