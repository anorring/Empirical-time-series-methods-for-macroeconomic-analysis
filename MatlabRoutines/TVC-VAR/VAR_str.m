%Creates matrices for VAR(k)
function [yy,x] = VAR_str(y,c,k)
s=size(y);
T=s(1); N=s(2);
for i=1:N,
yy(:,i)=y(k+1:T,i);
for j=1:k,
xx(:,k*(i-1)+j)=y(k+1-j:T-j,i);
end; end;
if c==0,xx=xx;
elseif c==1,xx=[ones(T-k,1) xx];
elseif c==2,xx=[ones(T-k,1) (1:T-k)' xx];
end
z(:,1)=xx(:,1);
for ij = 1:k
    zz(:,(ij-1)*N+1:N*ij) = xx(:,ij+c:k:size(xx,2));
end
if c==1,x = [z zz];
elseif c==0 x=zz;
end