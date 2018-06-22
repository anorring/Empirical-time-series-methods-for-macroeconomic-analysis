%%% Bayesian VAR - Sign Restrictions
%%% Author: Nicol� Maffei Faccioli

%% Importing the data:

clc; clear;

load 'newtrialAP.mat'
load 'APWW20081.mat'
load 'GDPMA.mat'
load 'GS10.mat' 
load 'tenyrs.mat'

% Asset Purchase Announcements:

DATES{83, 1}(1) = "01/11/2014";
DATES{84, 1}(1) = "01/12/2014";

figure1=figure(1);
plot(APWW20081(1:end-5),'LineWidth',3.5), axis tight
ylabel('US AP Announcements','FontWeight','Bold')
xlabel('Year','FontWeight','Bold')
legend('AP Announcements','Error','italic')
title('Plot of US Asset Purchase Announcements ($bn) ','FontWeight','Bold')
set(gca,'Xtick',1:5:length(APWW20081(1:end-5)),'XtickLabel',DATES(1:5:end)) % Changing what's on the x axis with the "dates" in form dd-mm-yy
set(gca,'FontSize',16)
xtickangle(90)
grid on

%% Estimating a Bayesian VAR(2) with Gibbs Sampler:

y=[log(GDPMA(1:end-24))*100, log(CPI(1:end-24))*100, APWW20081(1:end)./143.839, GS10(1:end-24), log(SP500(1:end-24))*100];

[T,n]=size(y);
p=2; % # of lags
c=1; % include constant | if c=0, no constant

% Draws set up:

sel=2;
maxdraws=3000;
b_sel=1001:sel:maxdraws; % exclude the first 1000 draws
drawfin=size(b_sel,2); % number of stored draws

% Set up the loop for each draw :

PI=zeros(n*p+c,n,maxdraws+1);
BigA=zeros(n*p,n*p,maxdraws+1);
Q=zeros(n,n,maxdraws+1);
errornorm=zeros(T-p,n,maxdraws+1);
fittednorm=zeros(T-p,n,maxdraws+1);

for i=1:maxdraws+1
    
    [PI(:,:,i),BigA(:,:,i),Q(:,:,i),errornorm(:,:,i),fittednorm(:,:,i)]=BVAR(y,p,c);
    
end

%% Sign Restricted IRFs:

hor=48;
restr=1; % # of restricted periods (1 = on impact) after the shock
candidateirf=zeros(n,n,hor,drawfin); % candidate impulse response

% Discard the first draws:

BigA=BigA(:,:,b_sel);
Q=Q(:,:,b_sel);

% Set up 4-D matrices for IRFs to be filled in the loop:

C=zeros(n,n,hor,drawfin); 
D=zeros(n,n,hor,drawfin);


h = waitbar(0,'Please wait, some fantastic results are coming ...');

for k=1:drawfin

for j=1:48
    BigC=BigA(:,:,k)^(j-1);
    C(:,:,j,k)=BigC(1:n,1:n); % Impulse response functions of the Wold representation
end

% Cholesky factorization:

S=chol(Q(:,:,k),'lower'); % lower triangular matrix

for i=1:48
    D(:,:,i,k)=C(:,:,i,k)*S; % Cholesky Wold respesentation
end

% Sign Restrictions Loop:

control=0;

while control==0
    
    % STEP 1:
   
    W=mvnrnd(zeros(n),eye(n)); % draw an independent standard normal nxn matrix
    [Qr,~]=qr_dec(W); % QR decomposition with the with the diagonal or R normalized to be positive
    
    % STEP 2 - Compute candidate IRFs:
    
    for i=1:hor
            candidateirf(:,:,i,k)=D(:,:,i,k)*Qr';
    end
    
    % STEP 3 - Check if the candidates satisfy all the sign restrictions:
    
    % N.B. : In what follows, following RWZ (2010), I use, to gain
    % efficiency, the fact that changing the sign of any of the columns of
    % matrix Q results in another orthogonal matrix. If all of the
    % candidate IRFs have wrong signs, then, by changing the sign of, say,
    % column j of Q, we obtain a new rotation matrix that, multiplied by D
    % will give us a candidate that satisfies the sign restrictions.

   %=== Responses to a Supply Shock:
   
   a = (candidateirf(1,1,1:restr,k) > 0) .* (candidateirf(2,1,1:restr,k) < 0) .* (candidateirf(4,1,1:restr,k) > 0) .* (candidateirf(5,1,1:restr,k) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock.
      am = (candidateirf(1,1,1:restr,k) < 0) .* (candidateirf(2,1,1:restr,k) > 0) .* (candidateirf(4,1,1:restr,k) < 0) .* (candidateirf(5,1,1:restr,k) < 0);
      
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Qr(1,:) = -Qr(1,:); 
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %=== Responses to a Demand shock:
   
   a = (candidateirf(1,2,1:restr,k) > 0) .* (candidateirf(2,2,1:restr,k) > 0) .* (candidateirf(4,2,1:restr,k) > 0) .* (candidateirf(5,2,1:restr,k) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      am = (candidateirf(1,2,1:restr,k) < 0) .* (candidateirf(2,2,1:restr,k) < 0) .* (candidateirf(4,2,1:restr,k) < 0) .* (candidateirf(5,2,1:restr,k) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Qr(2,:) = -Qr(2,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end

   %=== Responses to an AP shock:
 
   a = (candidateirf(1,3,1:restr,k) > 0) .* (candidateirf(2,3,1:restr,k) > 0) .* (candidateirf(3,3,1:restr,k) > 0).* (candidateirf(4,3,1:restr,k) < 0).* (candidateirf(5,3,1:restr,k) > 0);
   if (max(a)==0)
      %--- Swiching the sign of the shock and normalize.
      am = (candidateirf(1,3,1:restr,k) < 0) .* (candidateirf(2,3,1:restr,k) < 0) .* (candidateirf(3,3,1:restr,k) < 0).* (candidateirf(4,3,1:restr,k) > 0).* (candidateirf(5,3,1:restr,k) < 0);
      if (min(am)==0)
         continue;   %The restrictions are not satisfied.  Go the beginning to redraw.
      else
         %--- Normalizing according to the switched sign.
         Qr(3,:) = -Qr(3,:);
      end
   elseif (min(a)==0)
      continue;  %The restrictions are not satisfied.  Go the beginning to redraw.
   end
   
   %--- Terminating condition: all restrictions are satisfied.
   
   control=1;
   
end

   %--- Given the properly selected Q matrix, compute the responses that
   % satisfy all the sign restrictions and store these:
   
    for i=1:hor
            candidateirf(:,:,i,k)=D(:,:,i,k)*Qr';
    end
    
    waitbar(k/drawfin,h,sprintf('Please wait, cool results are on the way! Percentage completed %2.2f',(k/drawfin)*100))

end

%% Reshape the matrices into a 3D object:

% For each draw, compute a matrix with the IRFs for each variable and each
% shock for the entire horizon considered, i.e. hor x n*n for the # of
% draws:

candidateirf_wold=zeros(hor,n*n,drawfin); 

for k=1:drawfin
    
candidateirf_wold(:,:,k)=(reshape(permute(candidateirf(:,:,:,k),[3 2 1]),hor,n*n,[]));

end

% Create Probability Bands:

conf=68;

LowD=zeros(hor,n*n);
MiddleD=zeros(hor,n*n);
HighD=zeros(hor,n*n);

for k=1:n
 for j=1:n
        Dmin = prctile(candidateirf_wold(:,j+n*k-n,:),(100-conf)/2,3); %16th percentile
        LowD(:,j+n*k-n) = Dmin; %lower band
        Dmiddle=prctile(candidateirf_wold(:,j+n*k-n,:),50,3); %50th percentile
        MiddleD(:,j+n*k-n) = Dmiddle; %lower band
        Dmax = prctile(candidateirf_wold(:,j+n*k-n,:),(100+conf)/2,3); %84th percentile
        HighD(:,j+n*k-n) = Dmax; %upper band
 end
end

% Plot:

colorBNDS=[0.7 0.7 0.7];
VARnames={ 'GDP'; 'CPI';'AP Announcements'; '10 yrs'; 'Real Equity Prices'};

figure2=figure(2);

for k=1:n
    for j=1:n
    subplot(n,n,j+n*k-n)
    fill([0:hor-1 fliplr(0:hor-1)]' ,[HighD(:,j+n*k-n); flipud(LowD(:,j+n*k-n))],...
        colorBNDS,'EdgeColor','None'); hold on;
    plot(0:hor-1,MiddleD(:,j+n*k-n),'LineWidth',3.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[0 0 0],'LineStyle','-'); hold off;
    title(VARnames{k})
    %legend({'68% confidence bands','IRF'},'FontSize',12)
    set(gca,'FontSize',15)
    xlim([0 hor-1]);
    end
end



