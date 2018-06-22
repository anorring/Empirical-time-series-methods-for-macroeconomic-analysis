function [test1,test2,s]=dfactest(e,  N,T, xepsilon,method);
% program to compute the number of dynamic factors
% e is the residual matrix from estimating a VAR in Fhat
% Fhat is estimated using the L'L/N=I normalizationy


  if method==1; sigF=cov(e); end;
  if method==2; sigF=corrcoef(e); end;
r=cols(sigF);
stat1=zeros(1,r-1);
stat2=zeros(1,r-1);


CT=sqrt(min([N;T]));
[u,s,v]=svd(sigF);
s2=s.*s;
nf=sum(diag(s2));
W1=diag(s2);

for k=1:r-1;
  stat1(1,k)=sqrt(W1(k+1))/sqrt(nf);
  stat2(1,k)=sqrt(sum(W1(k+1:r)))/sqrt(nf);
end; % end k

test1=zeros(1,rows(xepsilon));
test2=zeros(1,rows(xepsilon));

for iepsilon=1:rows(xepsilon);
epsilon=xepsilon(iepsilon,:);
dum1=stat1 < epsilon(1);
dum2=stat2 < epsilon(2);



 test1(iepsilon)=r;
 test2(iepsilon)=r; 

if sum(dum1 )>0;
test1(iepsilon)=min(find(dum1==1));
end;
if sum(dum2 )>0;
test2(iepsilon)=min(find(dum2==1));
end;
end; % end iepsilon;

s=diag(s);





