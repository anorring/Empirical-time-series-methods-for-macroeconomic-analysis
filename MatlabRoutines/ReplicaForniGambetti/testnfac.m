clear;

r=3;
N=100;
T=50;
randn('state',999);

e=randn(T,N);
f=randn(T,r);
lambda=randn(N,r);
x=f*lambda'+e;

rmax=10;
DEMEAN=1;
disp(sprintf('Demean %d',DEMEAN));
disp('detrmining number of factors');
disp(sprintf('T= %d N= %d',size(x)));
for i=1:2;
  disp([nbpiid(x,rmax,i,DEMEAN)   nbplog(x,rmax,i,DEMEAN)]);
end;  
