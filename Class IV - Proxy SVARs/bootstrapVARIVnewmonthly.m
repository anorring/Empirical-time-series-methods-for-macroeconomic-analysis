function [HighH,LowH]=bootstrapVARIVnewmonthly(data,hor,c,cc,iter,conf,p,n,Z)

% Function to perform bootstrap method
% Author: Nicolo' Maffei Faccioli

% Step 1: 

[pi_hat,~,~,Y_initial,~,err]=VAR(data,p,c);

T=length(data)-size(Y_initial,1); % Data after losing lags

if c==1
    cons=pi_hat(1,:)';
    PI=pi_hat(2:end,:)';
else, cons=0;
    PI=pi_hat';
    
end

ExtraBigC=zeros(hor,n,iter);

for j=1:iter
    
% STEP 2: generate a new sample with bootstrap

% Wild bootstrap based on simple distribution (~Rademacher)
rr = 1-2*(rand(T,1)>0.5);
u=zeros(length(err),size(err,2));

for k=1:n
    u(:,k) = err(:,k).*rr;
end

Znew = Z.*rr;
Znew(isnan(Znew))=[];

Y_new=zeros(T,n);
Y_in= reshape(flipud(Y_initial)',1,[]);                  % 1*(n*p) row vector of [y_0 ... y_-p+1]
  for i=1:T
    Y_new(i,:)= cons' + Y_in*PI' + u(i,:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
    Data_new=[Y_initial ; Y_new];                        % Add the p initial lags to recreate the sample
  
% STEP 3: Esimate VAR(p) with new sample and obtain IRF:

[pi_hat_c,~,~,~,~,err_c]=VAR(Data_new,p,c);

eps_p=err_c(length(data)-p-length(Z)+1:end,end);
eps_q=err_c(length(data)-p-length(Z)+1:end,1:end-1);

ols_instr=OLS(eps_p,Znew,cc);

if cc==1
u_hat_p=[ones(length(Znew),1),Znew]*ols_instr;
else, u_hat_p=Znew*ols_instr;
end

sq_sp=(u_hat_p'*u_hat_p)\u_hat_p'*eps_q;

s=[sq_sp,1];

if c==1
BigA_c=[pi_hat_c(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np+1 x np+1 matrix
else, BigA_c=[pi_hat_c'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np x np matrix
end

C_c=zeros(n,n,hor);
H=zeros(n,1,hor);

for l=1:hor
    BigC_c=BigA_c^(l-1);
    C_c(:,:,l)=BigC_c(1:n,1:n); 
    H(:,:,l)=C_c(:,:,l)*s'; % Impulse response functions of the IV Wold representation (Bootstrap)
end

C_c_1= reshape(permute(H,[3 2 1]),hor,n,[]);
ExtraBigC(:,:,j)=C_c_1;

end

% Create bands:

LowH=zeros(hor,n);
HighH=zeros(hor,n);

for k=1:n
        Cmin = prctile(ExtraBigC(:,k,:),(100-conf)/2,3); %16th percentile
        LowH(:,k) = Cmin; %lower band
        Cmax = prctile(ExtraBigC(:,k,:),(100+conf)/2,3); %84th percentile
        HighH(:,k) = Cmax; %upper band
end

end

