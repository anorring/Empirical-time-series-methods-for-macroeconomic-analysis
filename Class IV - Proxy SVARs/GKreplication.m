%% Replication Trial:

clc,clear;

load 'GKdata.mat'

% load('trialdata.mat')

%% Estimate VAR and implement Cholesky :

finaldata=[log(CPI)*100, log(INDPRO)*100, DGS1, ebpm];

%% 2) Estimate a VAR(12) with finaldata:

% First of all, we have to create the T x n matrices of the SUR
% representation, i.e. Y and X:

n=size(finaldata,2); % # of variables
p=12; % 12 lags
c=1;

[pi_hat,~,~,~,~,err]=VAR(finaldata,p,1);
BigA=[pi_hat(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix

%% 3) Wold representation impulse responses:

C=zeros(n,n,50);

for j=1:50
    BigC=BigA^(j-1);
    C(:,:,j)=BigC(1:n,1:n); % Impulse response functions of the Wold representation
end

%% Cholesky:

T=length(finaldata)-n*p-n; % -p lags for n variables and -n constants
omega=(err'*err)./T; %estimate of omega
S=chol(omega,'lower'); %cholesky factorization, lower triangular matrix

D=zeros(n,n,50);
for i=1:50
    D(:,:,i)=C(:,:,i)*S; %cholesky wold respesentation
end

D_wold=(reshape(permute(D,[3 2 1]),50,n*n,[]));

% Bootstrap:

hor=50;
iter=1000;
conf=95;

[HighD,LowD]=bootstrapchol(finaldata,hor,c,iter,conf,p,n); % function that performs bootstrap

% Plot:

colorBNDS=[0.7 0.7 0.7];
VARnames={'One-Year Rate'; 'CPI'; 'IP';  'Excess Bond Premium'};

figure2=figure(2);

for j=[11,3,7,15;1:n]
    subplot(n,1,j(2))
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD(:,j(1)) ; flipud(LowD(:,j(1)))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,D_wold(:,j(1)),'LineWidth',3.5,'Color','k'); axis tight
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames(j(2)),'FontWeight','bold') 
    legend({'95% confidence bands','IRF'},'FontSize',15)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
end


%% Estimate VAR and implement GK :

clc; clear;

load 'GKdata.mat' 

finaldata=[log(CPI)*100, log(INDPRO)*100, ebpm, FEDFUNDS];

%% Estimate a VAR(12) with finaldata:

% First of all, we have to create the T x n matrices of the SUR
% representation, i.e. Y and X:

n=size(finaldata,2); % # of variables
p=12; % 12 lags
c=1; % constant

% First Step - Estimate VAR and store residuals:

[pi_hat,Y,X,Y_initial,Yfit,err]=VAR(finaldata,p,c);
BigA=[pi_hat(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix


%% 3) Wold representation impulse responses:

C=zeros(n,n,50);

for j=1:50
    BigC=BigA^(j-1);
    C(:,:,j)=BigC(1:n,1:n); % Impulse response functions of the Wold representation
end

%% External instrument:

% Second Step - Two stage regression:

eps_p=err(length(finaldata)-p-length(FF4)+1:end,end);
eps_q=err(length(finaldata)-p-length(FF4)+1:end,1:end-1);

Z = FF4;

[ols_instr,stderr_ols,tstat_ols]=OLS(eps_p,Z,1);

u_hat_p=[ones(length(Z),1),Z]*ols_instr;

sq_sp=(u_hat_p'*u_hat_p)\u_hat_p'*eps_q;

% Step 3 - Normalise s_p=1 and compute IRFs

s=[sq_sp,1];

H=zeros(n,1,50);
for i=1:50
    H(:,:,i)=C(:,:,i)*s'; % IV wold respesentation
end

H_wold=(reshape(permute(H,[3 2 1]),50,n,[]));

c=1;
cc=1;

Z = [NaN(length(finaldata)-p-length(FF4),1);FF4];

% Bootstrap:

hor=50;
iter=1000;
conf=95;

[HighH,LowH]=bootstrapVARIVnewmonthly(finaldata,hor,c,cc,iter,conf,p,n,Z);

% Plot:

colorBNDS=[0.7 0.7 0.7];
VARnames={'One-Year Rate'; 'CPI'; 'IP';  'Excess Bond Premium'};

figure3=figure(3);

for j=[4,1,2,3;1:n]
    subplot(n,1,j(2))
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighH(:,j(1)) ; flipud(LowH(:,j(1)))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,H_wold(:,j(1)),'LineWidth',3.5,'Color','k'); axis tight
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames(j(2)),'FontWeight','bold') 
    legend({'95% confidence bands','IRF'},'FontSize',15)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
end

