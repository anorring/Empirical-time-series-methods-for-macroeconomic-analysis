% [reg,u]=myvar(y,k,c)
% var multivariato con k lags;
% identificazione Wold; reg=parametri con st.err.; u=residui.
% c=0: ne' costante ne' trend; c=1: costante; c=2: costante e trend.
function [reg,u]=myvar(y,k,c);
if nargin==2, c=1; end
s=size(y);
T=s(1); N=s(2);
for i=1:N,
yy(:,i)=y(k+1:T,i);
for j=1:k,
xx(:,k*(i-1)+j)=y(k+1-j:T-j,i);
end; end;
if c==1, xx=[ones(T-k,1) xx];
elseif c==2,
xx=[ones(T-k,1) (1:T-k)' xx];
end
for i=1:N,
[reg(2*(i-1)+1:2*(i-1)+2,:),u(:,i)]=myols(yy(:,i),xx);
end;
