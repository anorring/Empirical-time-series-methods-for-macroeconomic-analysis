% Function to perform OLS
% Author: Nicol? Maffei Faccioli

function [beta_hat,stderr_ols,tstat_ols]=OLS(y,reg,c)

N=length(reg);

if c==1 % including the constant
    X=[ones(N,1) reg];
    beta_hat=(X'*X)\X'*y;
    
    % Standard errors:
    K=size(X,2);
    u=y-X*beta_hat;
    sigma2=(u'*u)/(N-K);
    stderr_ols=diag(sqrt(sigma2.*inv(X'*X)));
    tstat_ols=beta_hat./stderr_ols;     
    
else X=reg; % not including the constant
     beta_hat=(X'*X)\X'*y;
     
     % Standard errors and :
     K=size(X,2);
     u=y-X*beta_hat;
     sigma2=(u'*u)/(N-K);
     stderr_ols=diag(sqrt(sigma2.*inv(X'*X)));
     tstat_ols=beta_hat./stderr_ols;
        
end


end
