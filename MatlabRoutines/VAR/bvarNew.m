function [Bet,Sigma,B,S,N,nu]=bvarNew(data,nlagsvar,c,ndraws,B0,N0,S0,nu0)
if nargin==4
    B0=0;
    N0=0;
    nu0=0;
    S0=0;
end
[T n]=size(data);

% Y lhs variables, X regreessors
[Y,X] = VAR_str(data,c,nlagsvar);
np=size(X,2)*n;

%Ols estimates
Bols=inv(X'*X)*X'*Y;
res=Y-X*Bols;
Sols=cov(res);

% Posterior parameters
nu=T+nu0;
N=N0+X'*X;
B=inv(N)*(N0*B0+X'*X*Bols);
S=nu0/nu+T/nu*Sols+(1/nu)*(Bols-B0)'*N0*inv(N)*X'*X*(Bols-B0);

for i=1:ndraws
     PS = real(sqrtm(inv(nu*S))); 
     u = randn(n,nu);
     Sigma(:,:,i) = inv(PS*u*u'*PS');
%         Sigma(:,:,i)=wishrnd(inv(S)/nu,nu);
%         Sigma(:,:,i)=inv(Sigma(:,:,i));
         %Bet(:,i)=mvnrnd(vec(B,2),kron(Sigma(:,:,i),inv(N)),1)';
     [uu s vv]=svd((kron(Sigma(:,:,i),inv(N))));
     Bet(:,i)=reshape(B,np,1)+real(uu*sqrt(s)*vv')*randn(np,1);
end


