function [HighD,LowD]=bootstrapchol(data,hor,c,iter,conf,p,n)

% Function to perform bootstrap method
% Author: Nicolo' Maffei Faccioli

% Step 1: 

[pi_hat,~,~,Y_initial,~,err]=VAR(data,p,c);

T=length(data)-length(Y_initial); % Data after losing lags

if c==1
    cons=pi_hat(1,:)';
    PI=pi_hat(2:end,:)';
else cons=0;
    PI=pi_hat';
    
end

ExtraBigD=zeros(hor,n*n,iter);

for j=1:iter
    
% STEP 2: generate a new sample with bootstrap

Y_new=zeros(T,n);
Y_in= reshape(flipud(Y_initial)',1,[]);                  % 1*(n*p) row vector of [y_0 ... y_-p+1]
  for i=1:T;
    Y_new(i,:)= cons' + Y_in*PI' + err(randi(T),:);  
    Y_in=[Y_new(i,:) Y_in(1:n*(p-1))];
  end
    Data_new=[Y_initial ; Y_new];                        % Add the p initial lags to recreate the sample
  
% STEP 3: Esimate VAR(p) with new sample and obtain IRF:

[pi_hat_c,~,~,~,~,err_c]=VAR(Data_new,p,c);

if c==1
BigA_c=[pi_hat_c(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np+1 x np+1 matrix
else BigA_c=[pi_hat_c'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np x np matrix
end

t=length(Data_new)-n*p-n; % -p lags for n variables and -n constants
omega=(err_c'*err_c)./t; %estimate of omega
S=chol(omega,'lower')'; %cholesky factorization, lower triangular matrix

C_c=zeros(n,n,hor);

for l=1:hor
    BigC_c=BigA_c^(l-1);
    C_c(:,:,l)=BigC_c(1:n,1:n)*S; % Impulse response functions of the Wold representation (Bootstrap)
end

C_c_1= reshape(permute(C_c,[3 2 1]),hor,n*n,[]);
%ExtraBigC(:,:,j)=reshape(C_c_1,[hor p]);
ExtraBigD(:,:,j)=C_c_1;
end

% Create bands:

LowD=zeros(hor,n*n);
HighD=zeros(hor,n*n);

for k=1:n
 for j=1:n
        Dmin = prctile(ExtraBigD(:,j+n*k-n,:),(100-conf)/2,3); %16th percentile
        LowD(:,j+n*k-n) = Dmin; %lower band
        Dmax = prctile(ExtraBigD(:,j+n*k-n,:),(100+conf)/2,3); %84th percentile
        HighD(:,j+n*k-n) = Dmax; %upper band
 end
end

end

