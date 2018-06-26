%%% Empirical Time Series Methods for Macroeconomic Analysis
%%% Estimating VARs
%%% Author: Nicolo' Maffei Faccioli

%% Class Example:

clc; clear;

%% Import the data and create spread and GDP's growth rate:

load 'data.mat' 

GDPC1gr=diff(log(GDPC1)).*100; % real GDP's growth rate
spread=(GS10-FEDFUNDS); % spread between 10 yrs and FF

finaldata=[GDPC1gr,spread(2:end)];

% Plot:

DATE(1)=[]; % exclude the first date as we lost one observation by taking diff(log(.))

figure1=figure(1);
plot([GDPC1gr spread(2:end)],'LineWidth',2), axis('tight')
ylabel('GDPC1 growth rate and spread','FontWeight','Bold')
xlabel('Year','FontWeight','Bold')
legend('GDPC1 growth rate','spread GS10-FF','italic')
title('Plot of US GDP growth rate and GS10-FF spread ','FontWeight','Bold')
set(gca,'FontSize',15)
set(gca,'Xtick',1:20:length(GDPC1gr),'XtickLabel',DATE(1:20:end)) % Changing what's on the x axis with the "dates" in form dd-mm-yy
xtickangle(90)

%% Estimate a VAR(4) with finaldata:

% First of all, we have to create the T x n matrices of the SUR
% representation, i.e. Y and X:

n=2; % # of variables
p=4; % 4 lags
c=1; % including constant

[pi_hat,Y,X,Y_initial,Yfit,err]=VAR(finaldata,p,c); % VAR estimation, c - constant

pi_hat % the first row corresponds to the constant

BigA=[pi_hat(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, npxnp matrix
% start from 2 to exclude the constant, mistake --> error



%% Wold representation impulse responses:

C=zeros(n,n,20);

for j=1:20
    BigC=BigA^(j-1);
    C(:,:,j)=BigC(1:n,1:n); % Impulse response functions of the Wold representation
end

C_wold=reshape(permute(C,[3 2 1]),20,n*n,[]);  % "will make our life easier"

% Bootstrap:

hor=20;
iter=1000;
conf=68;

[HighC,LowC]=bootstrapVAR(finaldata,hor,c,iter,conf,p,n); % function that performs bootstrap

% Plot:

colorBNDS=[0.7 0.7 0.7];
VARnames={'Real GDP growth'; 'Spread'};

figure2=figure(2);

for k=1:n
    for j=1:n
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighC(:,j+n*k-n); flipud(LowC(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,C_wold(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;  % the confidence bands
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames{k})
    legend({'68% confidence bands','IRF'},'FontSize',12)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
    end
end


%% One period ahead forecasts:

% First of all, we have to estimate the parameters using an initial sample
% of T_0=123:

T_0=123;
ldiff=length(finaldata)-T_0;
forecasts1=zeros(ldiff,n);

% Notice that here I make a loop that repeats the steps until the end of
% the sample:

for k=1:ldiff
    
    [pi_hat,Y,X,Y_initial,Yfit,err]=VAR(finaldata(1:123+k-1,:),p,1);
    Ylast= reshape(flipud(Y(end-p+1:end,:))',1,[]); 
    forecasts1(k,:)=pi_hat(1,:)+Ylast*pi_hat(2:end,:); % 1 period ahead
    
end

y_hat11=[NaN(123,1); forecasts1(:,1)]; % here I include vectors of NaN in order to include both actual and forecasts in the same graph
y_hat21=[NaN(123,1); forecasts1(:,2)];


figure3=figure(3);
plot([GDPC1gr y_hat11],'LineWidth',2), axis('tight')
ylabel('GDP growth rate','FontWeight','Bold')
xlabel('Year','FontWeight','Bold')
legend('GDP growth rate','VAR(4) forecast 1 period ahead','italic')
title('Plot of US GDP growth rate and VAR(4) forecasts one periods ahead','FontWeight','Bold')
set(gca,'FontSize',15)
set(gca,'Xtick',1:20:length(GDPC1gr),'XtickLabel',DATE(1:20:end)) % Changing what's on the x axis with the "dates" in form dd-mm-yy
xtickangle(90)

figure4=figure(4);
plot([spread(2:end) y_hat21],'LineWidth',2), axis('tight')
ylabel('spread GS10-FF','FontWeight','Bold')
xlabel('Year','FontWeight','Bold')
legend('spread GS10-FF','VAR(4) forecast 1 period ahead','italic')
title('Plot of spread GS10-FF and VAR(4) forecasts one periods ahead ','FontWeight','Bold')
set(gca,'FontSize',15)
set(gca,'Xtick',1:20:length(spread(2:end)),'XtickLabel',DATE(1:20:end)) % Changing what's on the x axis with the "dates" in form dd-mm-yy
xtickangle(90)






