%%% Empirical Time Series Methods for Macroeconomic Analysis
%%% Structural VARs
%%% Author: Nicolo' Maffei Faccioli

%% Class Example:

clc; clear;

%% Import the data and create spread and GDP's growth rate:

load 'datastruc.mat' 

finaldata=[TFPgr(2:end),diff(stockprices)]; % stockprices is already in logs

%% Estimate a VAR(4) with finaldata:

% First of all, we have to create the T x n matrices of the SUR
% representation, i.e. Y and X:

n=2; % # of variables
p=4; % 4 lags
c=1; % including constant

[pi_hat,Y,X,Y_initial,Yfit,err]=VAR(finaldata,p,c); % VAR estimation

BigA=[pi_hat(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix

%% Wold representation impulse responses + Cholesky:

C=zeros(n,n,20);

for j=1:20
    BigC=BigA^(j-1);
    C(:,:,j)=BigC(1:n,1:n); % Impulse response functions of the Wold representation
end

C_wold=reshape(permute(C,[3 2 1]),20,n*n,[]);

% Bootstrap:

hor=20;
iter=1000;
conf=68;

[HighC,LowC]=bootstrapVAR(finaldata,hor,c,iter,conf,p,n); % function that performs bootstrap

colorBNDS=[0.7 0.7 0.7];
VARnames={'TFP growth'; 'SP500 growth'};

figure2=figure(2);

for k=1:n
    for j=1:n
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighC(:,j+n*k-n); flipud(LowC(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,C_wold(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames{k})
    legend({'68% confidence bands','IRF'},'FontSize',12)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
    end
end

% Cholesky:

T=length(Y)-n*p-n; % -p lags for n variables and -n constants
omega=(err'*err)./T; %estimate of omega
S=chol(omega,'lower'); %cholesky factorization, lower triangular matrix

D=zeros(n,n,20);
for i=1:20
    D(:,:,i)=C(:,:,i)*S; %cholesky wold respesentation
end

D_wold=(reshape(permute(D,[3 2 1]),20,n*n,[]));

[HighD,LowD]=bootstrapchol(finaldata,hor,c,iter,conf,p,n); % function that performs bootstrap for cholesky

% Plot:

figure3=figure(3);

for k=1:n
    for j=1:n
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD(:,j+n*k-n); flipud(LowD(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,D_wold(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames{k})
    legend({'68% confidence bands','IRF'},'FontSize',12)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
    end
end

%% News Shock:

eta=S\err';

% Plot:

DATES(1)=[];

figure4=figure(4);
plot([eta(2,:)' err(:,2)],'LineWidth',1), axis('tight')
ylabel('News shock and Error','FontWeight','Bold')
xlabel('Year','FontWeight','Bold')
legend('Cholesky','Error','italic')
title('Plot of News Shock (Cholesky) vs. Error ','FontWeight','Bold')
set(gca,'Xtick',1:60:length(eta),'XtickLabel',DATES(1:60:end)) % Changing what's on the x axis with the "dates" in form dd-mm-yy

%% Blanchard & Quah - Long-run restriction:

A1=sum(C,3);
S=chol(A1*omega*A1')';
K=A1\S;
F=zeros(2,2,20);

for i=1:20
    F(:,:,i)=C(:,:,i)*K;   
end

F_wold=reshape(permute(F, [3 2 1]),[],n*n);

[HighF,LowF]=bootstrapBQ(finaldata,hor,1,iter,conf,p,n);

figure5=figure(5);

for k=1:n
    for j=1:n
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighF(:,j+n*k-n); flipud(LowF(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,F_wold(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames{k})
    legend({'68% confidence bands','IRF'},'FontSize',12)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
    end
end

w=K\err';

figure6=figure(6);
plot([w(2,:)' err(:,2)],'LineWidth',1), axis('tight')
ylabel('News shock and Error','FontWeight','Bold')
xlabel('Year','FontWeight','Bold')
legend('B&Q','Error','italic')
title('Plot of News Shock (B&Q) vs. Error ','FontWeight','Bold')
set(gca,'Xtick',1:60:length(eta),'XtickLabel',DATES(1:60:end)) % Changing what's on the x axis with the "dates" in form dd-mm-yy

%% Cumulative Impulse Responses:

% Cholesky:

D_wold_cum=cumsum(D_wold);
HighD_cum=cumsum(HighD);
LowD_cum=cumsum(LowD);

VARnames={'TFP'; 'SP500'};

figure7=figure(7);

for k=1:n
    for j=1:n
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD_cum(:,j+n*k-n); flipud(LowD_cum(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,D_wold_cum(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames{k})
    legend({'68% confidence bands','IRF'},'FontSize',12)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
    end
end

% Blanchard and Quah:

F_wold_cum=cumsum(F_wold);
HighF_cum=cumsum(HighF);
LowF_cum=cumsum(LowF);

figure8=figure(8);

for k=1:n
    for j=1:n
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighF_cum(:,j+n*k-n); flipud(LowF_cum(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,F_wold_cum(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames{k})
    legend({'68% confidence bands','IRF'},'FontSize',12)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
    end
end

%% Variance Decomposition Analysis:

% Cholesky:

F_TFP=sum((D_wold(:,1)).*(D_wold(:,1)))+sum((D_wold(:,2)).*(D_wold(:,2)));
F_SP=sum((D_wold(:,3)).*(D_wold(:,3)))+sum((D_wold(:,4)).*(D_wold(:,4)));
decomposion_TFP1=sum((D_wold(:,1)).*(D_wold(:,1)))/F_TFP;
decomposion_TFP2=sum((D_wold(:,2)).*(D_wold(:,2)))/F_TFP;
decomposion_SP1=sum((D_wold(:,3)).*(D_wold(:,3)))/F_SP;
decomposion_SP2=sum((D_wold(:,4)).*(D_wold(:,4)))/F_SP;

table(decomposion_TFP1,decomposion_TFP2)
table(decomposion_SP1, decomposion_SP2)

% Blanchard & Quah:

F_TFPBQ=sum((F_wold(:,1)).*(F_wold(:,1)))+sum((F_wold(:,2)).*(F_wold(:,2)));
F_SPBQ=sum((F_wold(:,3)).*(F_wold(:,3)))+sum((F_wold(:,4)).*(F_wold(:,4)));
decomposion_TFP1BQ=sum((F_wold(:,1)).*(F_wold(:,1)))/F_TFPBQ;
decomposion_TFP2BQ=sum((F_wold(:,2)).*(F_wold(:,2)))/F_TFPBQ;
decomposion_SP1BQ=sum((F_wold(:,3)).*(F_wold(:,3)))/F_SPBQ;
decomposion_SP2BQ=sum((F_wold(:,4)).*(F_wold(:,4)))/F_SPBQ;

table(decomposion_TFP1BQ,decomposion_TFP2BQ)
table(decomposion_SP1BQ, decomposion_SP2BQ)

%% Variance Decomposition:

%MiddleDsquare=D_wold.^2; % For Cholesky, growth rates
MiddleDsquare=F_wold_cum.^2; % For B&Q, growth rates

%MiddleDsquare=D_wold_cum.^2; % For Cholesky, levels
%MiddleDsquare=F_wold_cum.^2; % For B&Q, levels

denom=zeros(hor,n);

for k=[1:n;1:n:n*n]
    
    denom(:,k(1))=cumsum(sum(MiddleDsquare(:,k(2):k(2)+n-1),2));
    
end

denomtot=zeros(hor,n*n);

for k=[1:n;1:n:n*n]
    
    denomtot(:,k(2):k(2)+n-1)=denom(:,k(1)).*ones(hor,n);
    
end

vardec=zeros(hor,n);

for j=1:n*n
    
    vardec(:,j)=cumsum(MiddleDsquare(:,j))./denomtot(:,j);
    
end

figure2=figure(2);

% VARnames={'TFP growth'; 'SP500 growth'}; % growth rates
VARnames={'TFP'; 'SP500'}; % levels



for k=[1:n;1:n:n*n]
    subplot(1,n,k(1))
    area(0:hor-1,vardec(:,k(2):k(2)+n-1),'LineWidth',1); hold on; 
    legend({'TFP shock','SP500 shock'},'FontSize',13)
    set(gca,'FontSize',15)
    title(VARnames{k(1)})
    xlim([0 hor-1]);
    ylim([0 1])
end

