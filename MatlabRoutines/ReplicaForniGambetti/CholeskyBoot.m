function [cirf csirfBoot] = CholeskyBoot(y,p,H,c,MaxBoot,cl,ind,opt)
if nargin<8
    opt=0;
end
n=size(y,2); % number of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OLS estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Y X] = VAR_str(y,c,p);       % yy and XX all the sample
T=size(Y,1);
Bols=inv(X'*X)*X'*Y;
VecB=reshape(Bols,n*(n*p+1),1);
B=companion(VecB,p,n,c);
%C=VecB(1:n*p+c:end);
C=Bols(1,:)';
%%%% residuals
% for t=1:T
%     XX(:,:,t)=kron(eye(n),X(t,:)');
%     u(:,t)=Y(t,:)'-XX(:,:,t)'*VecB;
% end
u=Y-X*Bols;
CovU=cov(u);
u=u';
%%%% impulse response functions
for h=1:H
    irf(:,:,h)=B^(h-1);
    cirf(:,:,h)=irf(1:n,1:n,h)*chol(CovU)';
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
    CovUBoot(:,:,i)=cov(Yn-Xn*Bolsn);
    VecBn=reshape(Bolsn,n*(n*p+1),1);
    Bn=companion(VecBn,p,n,c);
   
    %%%% impulse response functions
    for h=1:H
        irfBoot(:,:,h,i)=Bn^(h-1);
        %%%% cholesky
        CholBoot(:,:,h,i)=irfBoot(1:n,1:n,h,i)*chol(CovUBoot(:,:,i))';
    end
end
lb=round(MaxBoot*(1-cl)/2);
ub=round(MaxBoot*(cl+(1-cl)/2));
me=round(MaxBoot*.5);

cirf(ind,:,:)=cumsum(cirf(ind,:,:),3);
CholBoot(ind,:,:,:)=cumsum(CholBoot(ind,:,:,:),3);
csirfBoot=sort(CholBoot,4);

if opt==1
    k=0;
    figure(2)
    for ii=1:n
        for jj=1:n
            k=k+1;
            subplot(n,n,k),plot(1:H,squeeze(cirf(ii,jj,:)),'k',...
            1:H,squeeze(csirfBoot(ii,jj,:,[lb ub])),':k'),axis tight%,1:H,squeeze(mean(csirfBoot(ii,jj,:,:),4)),'--k'),axis tight
        end
    end
end
