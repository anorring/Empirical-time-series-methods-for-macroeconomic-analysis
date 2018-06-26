function [HighC,LowC]=bootstrapVAR(data,hor,c,iter,conf,p,n)

% Function to perform bootstrap method
% Author: Nicolo' Maffei Faccioli

% Step 1: estimate and store

[pi_hat,~,~,Y_initial,~,err]=VAR(data,p,c);

T=length(data)-length(Y_initial); % Data after losing lags

if c==1
    cons=pi_hat(1,:)';
    PI=pi_hat(2:end,:)';
else cons=0;
    PI=pi_hat';
end

ExtraBigC=zeros(hor,n*n,iter);

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

[pi_hat_c,~,~,~,~,~]=VAR(Data_new,p,c);

if c==1
BigA_c=[pi_hat_c(2:end,:)'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np+1 x np+1 matrix
else BigA_c=[pi_hat_c'; eye(n*p-n) zeros(n*p-n,n)]; % BigA companion form, np x np matrix
end

C_c=zeros(n,n,hor);

for l=1:hor
    BigC_c=BigA_c^(l-1);
    C_c(:,:,l)=BigC_c(1:n,1:n); % Impulse response functions of the Wold representation (Bootstrap)
end

C_c_1= reshape(permute(C_c,[3 2 1]),hor,n*n,[]);
ExtraBigC(:,:,j)=C_c_1;

end

% Create bands:

LowC=zeros(hor,n*n);
HighC=zeros(hor,n*n);

for k=1:n
 for j=1:n
        Cmin = prctile(ExtraBigC(:,j+n*k-n,:),(100-conf)/2,3); %16th percentile
        LowC(:,j+n*k-n) = Cmin; %lower band
        Cmax = prctile(ExtraBigC(:,j+n*k-n,:),(100+conf)/2,3); %84th percentile
        HighC(:,j+n*k-n) = Cmax; %upper band
 end
end

end

