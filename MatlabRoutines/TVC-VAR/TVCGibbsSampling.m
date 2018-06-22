function [PB, PR, PQ, Q0]=TVCGibbsSampling(Y,L,maxgbit,burn,jump,para,years,init)

%%%% Load parameters of KSC(98)
MixingKSC

%%%% Set parameters
gbit = 1;
post_index = 1;
trs=1;
trsMax=10;
MaxEig=1.1;
MaxRoot=MaxEig+1;
c=1; % constant term

T0=4*years; % period in the initial sample 4*8-L
[yy xx] = VAR_str(Y(1:end,:),c,L);     % yy and XX all the sample
[TT r]=size(yy);

for t=1:TT
    XX(:,:,t)=kron(eye(r),xx(t,:)');
end

% y=yy(T0/2+1:end,:);
% x=xx(T0/2+1:end,:);
% X=XX(:,:,T0/2+1:end);
y=yy(init:end,:);
x=xx(init:end,:);
X=XX(:,:,init:end);
r=size(y,2);
[T rp] = size(x);

%%%% OLS in the pre sample
[B0,VB0,R0] = sur(yy',XX,T0);
resid=innovm(yy(1:T0,:)',xx(1:T0,:),B0*ones(1,T0),r,T0,L,c)'; % residuals

%%%% Variance and mean of theta0
VB0=1*VB0;

%%%% Variance and mean of log(sigma0)
H0=log(diag(sqrt(R0)));
VH0=eye(r);%diag(diag(R0))*2;


%%%% Mean and variance of alpha0

co1=[];
co2=[];
for i=1:r-1
    dd=size(co1,1)+1:size(co1,1)+size(co2,1)+1;
    A0(dd)=inv(resid(:,1:i)'*resid(:,1:i))*resid(:,1:i)'*resid(:,i+1);
    VA0(dd,dd) = inv(resid(:,1:i)'*resid(:,1:i))*cov(resid(:,i+1)-resid(:,1:i)*A0(dd)');
    co2=ones(i,1);
    co1=[co1;co2];
end
A0=A0';
% covariance matrix of innovation in the off-diagonal elements
U0=para(2)*VA0;
dfU0=[2:r];


%%%% Scale matrix for the variance in coefficients innovations
Q0=para(1)*VB0;
% Q0(1,:)=0.5*Q0(1,:);
% Q0(:,1)=0.5*Q0(:,1);
% Q0(9,:)=0.5*Q0(9,:);
% Q0(:,9)=0.5*Q0(:,9);

df0=size(B0,1)+1;
TQ0=df0*Q0;
df = T+df0;

%%%% Scale matrix for the variance in volatilities innovations
W0=para(3)*eye(r);

dfW0=r+1;

%%%% Initialize Yhat, ystar,y2star
dy = diff(y(:,1:r));
e(:,1:r) = dy - ones(T-1,1)*mean(dy);
Yhat=[zeros(r,1) e']';
InitS=cov(Yhat);
InitA=chol(InitS)';
ystar=inv(InitA*inv(diag(sqrt(diag(InitS)))))*Yhat';
y2star=[log(exp(H0).^2'); log(ystar(:,2:end)'.^2)];

%%%% Initialize states
for i=1:r
    for t=1:T
        for j=1:7
            probs(j)=ksc(j,2)*normpdf(y2star(t,i),y2star(t,i)+ksc(j,3)-1.2704,ksc(j,4));
        end
        s0(i,t)=DiscreteDraw(probs./sum(probs));
    end
end
s=s0;

%%%% Initialize variances
U=U0;
W=W0;
Q=Q0;
%R=R0;

%%%% Coefficitent matrix in the state space for volatilities
h=eye(r)*2;

%%%%%%%% Initial draws: this is useful for the last part of the program
P1(:,:,1) = VH0; % P(1|0)
K = (P1(:,:,1))*h'*inv(h* P1(:,:,1)*h' + diag(ksc(s(:,1),4))); % K(1)
P0(:,:,1) = P1(:,:,1) - K*h*P1(:,:,1); % P(1|1)
S0(:,1) = H0 + K*( y2star(1,:)' - h*H0 - (ksc(s(:,1),3)-1.2074)); % S(1|1)
for i = 2:T,
    P1(:,:,i) = P0(:,:,i-1) + W; % P(t|t-1)
    K = (P1(:,:,i))*h'*inv(h*P1(:,:,i)*h' + diag(ksc(s(:,i),4))); % K(t)
    P0(:,:,i) = P1(:,:,i) - K*h*(P1(:,:,i)); % P(t|t)
    S0(:,i) = S0(:,i-1) + K*( y2star(i,:)' - h*S0(:,i-1) - (ksc(s(:,i),3)-1.2704) ); % S(t|t)
end
wa = randn(r,T);
% SA(:,T) = S0(:,T) + real(sqrtm(P0(:,:,T)))*wa(:,T);
SA(:,T) = S0(:,T) + chol(P0(:,:,T))'*wa(:,T);
for i = 1:T-1,
    PM = P0(:,:,T-i)*inv(P1(:,:,T-i+1));
    P = P0(:,:,T-i) - PM*P0(:,:,T-i);
    SM = S0(:,T-i) + PM*(SA(:,T-i+1) - S0(:,T-i));
%     SA(:,T-i) = SM + real(sqrtm(P))*wa(:,T-i);
    SA(:,T-i) = SM + chol(P)'*wa(:,T-i);
end
lh=SA';
H=exp(2*SA'); %% compute variance fr

%%%% 2) Draw covariances
TA=zeros(r*(r-1)/2,r*(r-1)/2);
A=[];
AA=[];
for i=1:r-1
%     dd=i*(i-1)/2+1:i*(i+1)/2;
    dd=size(A,1)+1:size(A,1)+size(AA,1)+1;
    [AS0,AP0,AP1] = kfR1RegTvc(Yhat(:,i+1)',-Yhat(:,1:i)',U(dd,dd),H(:,i+1),A0(dd),VA0(dd,dd),T);
    AA = gibbs1Reg(Yhat(:,i+1),-Yhat(:,1:i),AS0,AP0,AP1,T,length(dd));
    [TU(dd,dd),dfU] = iwpQ(AA,r,T,L,dfU0(i)*U0(dd,dd),dfU0(i));
    A=[A;AA];
    g=real(sqrtm(inv(TU(dd,dd))));
    ug=randn(i,dfU);
    U(dd,dd) = inv(g*ug*ug'*g');
    CH=A;
    clear AS0 AP0 AP1;
end

for t=1:T
    CF(:,:,t)= chofac(r,A(:,t));
    R0(:,:,t) = inv(CF(:,:,t))*diag(H(t,:))*inv(CF(:,:,t))';
end
R=R0;

%%%% initial draw of states for AR coefficients
[SS0,PP0,PP1] = kfR(y',X,Q0,R,B0,VB0,T,r,L,c);
while MaxRoot>=MaxEig
    %[S0,P0,P1] = kfR(y',X,Q0,CF,H,B0,VB0,T,r,L);
    th_backward0 = gibbs1(y',X,SS0,PP0,PP1,T,r,L,c);
    %[th_backward0,P_backward0] = tvc_backward_filter(y,X,B0,VB0,R0,Q0,0);
    for j =1:T
        comp_mat = companion(th_backward0(:,j),L,r,c);
        AutVal(j) =  max(abs(eig(comp_mat)));
    end
    MaxRoot=max(AutVal)
end
B=th_backward0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Initialize matrices of draws
PQ=zeros(rp*r,rp*r,(maxgbit-burn)/jump);
PR=zeros(r,r,T,(maxgbit-burn)/jump);
PB=zeros(rp*r,T,(maxgbit-burn)/jump);
PW=zeros(r,r,(maxgbit-burn)/jump);
PU=zeros(r*(r-1)/2,r*(r-1)/2,(maxgbit-burn)/jump);
Ps=zeros(r,T,(maxgbit-burn)/jump);
PA=zeros(r*(r-1)/2,T,(maxgbit-burn)/jump);

PQ(:,:,post_index)=Q;
PR(:,:,:,post_index)=R;
PB(:,:,post_index)=B;
PW(:,:,post_index)=W;
PU(:,:,post_index)=U;
Ps(:,:,post_index)=s;
PA(:,:,post_index)=A;


%%%% Start the chain
disp('begin MCMC')

while gbit<=maxgbit
    gbit
    itnum=num2str(gbit);
    if itnum(size(itnum,2))==num2str(0) & itnum(size(itnum,2)-1)==num2str(0)
        disp(['Iteration   ' num2str(gbit) '    Collected Draws   ' num2str(post_index)...
            '   Done   ' num2str(gbit/maxgbit*100) '%    MaxEig   ' num2str(max(AutVal)) ])
    end

    %%%% 1) Draw stochastic volatilities
    P1(:,:,1) = VH0; % P(1|0)
    K = (P1(:,:,1))*h'*inv(h* P1(:,:,1)*h' + diag(ksc(s(:,1),4))); % K(1)
    P0(:,:,1) = P1(:,:,1) - K*h*P1(:,:,1); % P(1|1)
    S0(:,1) = H0 + K*( y2star(1,:)' - h*H0 - (ksc(s(:,1),3)-1.2074)); % S(1|1)
    for i = 2:T,
        P1(:,:,i) = P0(:,:,i-1) + W; % P(t|t-1)
        K = (P1(:,:,i))*h'*inv(h*P1(:,:,i)*h' + diag(ksc(s(:,i),4))); % K(t)
        P0(:,:,i) = P1(:,:,i) - K*h*(P1(:,:,i)); % P(t|t)
        S0(:,i) = S0(:,i-1) + K*( y2star(i,:)' - h*S0(:,i-1) - (ksc(s(:,i),3)-1.2704) ); % S(t|t)
    end
    wa = randn(r,T);
%     SA(:,T) = S0(:,T) + real(sqrtm(P0(:,:,T)))*wa(:,T);
    SA(:,T) = S0(:,T) + chol(P0(:,:,T))'*wa(:,T);

    for i = 1:T-1,
        PM = P0(:,:,T-i)*inv(P1(:,:,T-i+1));
        P = P0(:,:,T-i) - PM*P0(:,:,T-i);
        SM = S0(:,T-i) + PM*(SA(:,T-i+1) - S0(:,T-i));
%         SA(:,T-i) = SM + real(sqrtm(P))*wa(:,T-i);
        SA(:,T-i) = SM + chol(P)'*wa(:,T-i);
    end
    lh=SA';
    H=exp(2*SA'); %% compute variance fr

    %%%% 2) Draw covariances
    TA=zeros(r*(r-1)/2,r*(r-1)/2);
    A=[];
    AA=[];
    for i=1:r-1
%         dd=i*(i-1)/2+1:i*(i+1)/2;
        dd=size(A,1)+1:size(A,1)+size(AA,1)+1;
        [AS0,AP0,AP1] = kfR1RegTvc(Yhat(:,i+1)',-Yhat(:,1:i)',U(dd,dd),H(:,i+1),A0(dd),VA0(dd,dd),T);
        AA = gibbs1Reg(Yhat(:,i+1),-Yhat(:,1:i),AS0,AP0,AP1,T,length(dd));
        [TU(dd,dd),dfU] = iwpQ(AA,r,T,L,dfU0(i)*U0(dd,dd),dfU0(i));
        A=[A;AA];
        g=real(sqrtm(inv(TU(dd,dd))));
        ug=randn(i,dfU);
        U(dd,dd) = inv(g*ug*ug'*g');
        CH=A;
        clear AS0 AP0 AP1;
    end
    for t=1:T
        CF(:,:,t)= chofac(r,A(:,t));
        R(:,:,t) = inv(CF(:,:,t))*diag(H(t,:))*inv(CF(:,:,t))';
    end


    %%%% Draw VAR coefficients
    [SS0,PP0,PP1] = kfR(y',X,Q,R,B0,VB0,T,r,L,c);
    B = gibbs1(y',X,SS0,PP0,PP1,T,r,L,c);


    %%%% compute new Yhat
    Yhat=innovm(y',x,B,r,T,L,c)';

    %%%% Compute new ystar
    for t=1:T
        ystar(:,t)=CF(:,:,t)*Yhat(t,:)';
    end

    %%%% Compute new y2star
    y2star=log(ystar'.^2+0.001);

    %%%% Compute new probabilities
    for i=1:r
        for t=1:T
            for j=1:7
                probs(j)=ksc(j,2)*normpdf(y2star(t,i),2*SA(i,t)+ksc(j,3)-1.2704,ksc(j,4));
            end
            s(i,t)=DiscreteDraw(probs./sum(probs));
        end
    end

    %%%% Draw SV innovation variance
    [TW,dfW] = iwpQ(SA,r,T,L,dfW0*W0,dfW0);
    W = gibbs2Q(TW,dfW,r,0,1);

    %%%% Draw coefficients innovations variance
    [TQ,df] = iwpQ(B,r,T,L,TQ0,df0);
    Q = gibbs2Q(TQ,df,r,L,c);

    %%%% Check roots
    rootInd=1;
    timInd=T+1;
    while rootInd==1 & timInd<=T
        comp_mat = companion(B(:,timInd),L,r,c);
        AutVal =  max(abs(eig(comp_mat)));
        if max(AutVal)>MaxEig & trs<trsMax
            disp('stuck')
            rootInd=0;
            trs=trs+1;
        elseif max(AutVal)>MaxEig & trs==trsMax
            disp('go back')
            Q=PQ(:,:,post_index);
            R=PR(:,:,:,post_index);
            B=PB(:,:,post_index);
            W=PW(:,:,post_index);
            U=PU(:,:,post_index);
            s=Ps(:,:,post_index);
            A=PA(:,:,post_index);
            Yhat=innovm(y',x,B,r,T,L,c)';
            for t=1:T
                ystar(:,t)=chofac(r,A(:,t))*Yhat(t,:)';
            end
            y2star=log(ystar'.^2+0.001);
            rootInd=0;
            trs=0;
            %disp('eccomi')
        elseif max(AutVal)<MaxEig
            timInd=timInd+1;
        end
    end

    % sample from markov chain
    if timInd==T+1
        %B=th_b;
        if gbit <= burn
            post_index=1;
            PQ(:,:,post_index)=Q;
            PR(:,:,:,post_index)=R;
            PB(:,:,post_index)=B;
            PW(:,:,post_index)=W;
            PU(:,:,post_index)=U;
            Ps(:,:,post_index)=s;
            PA(:,:,post_index)=A;
        elseif gbit > burn
            check_index = (gbit-burn)/jump- fix((gbit-burn)/jump);
            if check_index == 0
                post_index = (gbit-burn)/jump+1;
            end
            PQ(:,:,post_index)=Q;
            PR(:,:,:,post_index)=R;
            PB(:,:,post_index)=B;
            PW(:,:,post_index)=W;
            PU(:,:,post_index)=U;
            Ps(:,:,post_index)=s;
            PA(:,:,post_index)=A;
        end
        gbit = gbit+1;
    end
end
