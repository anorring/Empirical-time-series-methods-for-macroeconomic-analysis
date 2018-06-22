function [bqirf bqsirfBoot] = BQBoot(y,p,H,c,MaxBoot,cl,ind,opt)
n=size(y,2); % number of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OLS estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Y X] = VAR_str(y,c,p);       % yy and XX all the sample
T=size(Y,1);
Bols=inv(X'*X)*X'*Y;
VecB=reshape(Bols,n*(n*p+1),1);
B=companion(VecB,p,n,c);
C=Bols(1,:)';
u=Y-X*Bols;
CovU=cov(u);
u=u';
%%%% impulse response functions
im=inv(eye(size(B,1))-B);
SSigma=zeros(n*p,n*p);
SSigma(1:n,1:n)=CovU;
CS=im*SSigma*im';
CC=chol(CS(1:n,1:n))';
S=inv(im(1:n,1:n))*CC;

for h=1:H
    irf(:,:,h)=B^(h-1);
    bqirf(:,:,h)=irf(1:n,1:n,h)*S;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Bootstrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:MaxBoot
    %i
    %%%% generate new series
    for t=1:T
        if t==1
            YB(:,1)=X(1,2:end)';
        else
            a=fix(rand(1,1)*size(u,2))+1;
            uu=u(:,a);
            YB(:,t)=[C;zeros(n*p-n,1)]+B*YB(:,t-1)+[uu;zeros(n*p-n,1)];
        end
    end

    %%%% estimate the new VAR
    yb=YB(1:n,:)';
    [Yn Xn] = VAR_str(yb,c,p);     % yy and XX all the sample

    T=size(Yn,1);
    Bolsn=inv(Xn'*Xn)*Xn'*Yn;
    CovUBoot=cov(Yn-Xn*Bolsn);
    VecBn=reshape(Bolsn,n*(n*p+1),1);
    Bn=companion(VecBn,p,n,c);
   
    imBoot=inv(eye(size(Bn,1))-Bn);
    SSigmaBoot=zeros(n*p,n*p);
    SSigmaBoot(1:n,1:n)=CovUBoot;
    CSBoot=imBoot*SSigmaBoot*imBoot';
    CCBoot=chol(CSBoot(1:n,1:n))';
    SBoot=inv(imBoot(1:n,1:n))*CCBoot;

    %%%% impulse response functions
    for h=1:H
        irfBoot(:,:,h,i)=Bn^(h-1);
        %%%% cholesky
        BQBoot(:,:,h,i)=irfBoot(1:n,1:n,h,i)*SBoot;
    end
end
lb=round(MaxBoot*(1-cl)/2);
ub=round(MaxBoot*(cl+(1-cl)/2));
me=round(MaxBoot*.5);

bqirf(ind,:,:)=cumsum(bqirf(ind,:,:),3);
BQBoot(ind,:,:,:)=cumsum(BQBoot(ind,:,:,:),3);
bqsirfBoot=sort(BQBoot,4);

if opt==1
    k=0;
    figure(2)
    for ii=1:n
        for jj=1:n
            k=k+1;
            subplot(n,n,k),plot(1:H,squeeze(bqirf(ii,jj,:)),'k',...
            1:H,squeeze(bqsirfBoot(ii,jj,:,[lb ub])),':k'),axis tight%,1:H,squeeze(mean(csirfBoot(ii,jj,:,:),4)),'--k'),axis tight
        end
    end
end
