%%% BVAR - Diffuse/Uninformative Prior (Uniform -inf, +inf)
%%% Author: Nicolo' Maffei Faccioli

function [PI,BigA,Q,errornorm,fittednorm]=BVAR(y,p,c)

    [Traw,K]=size(y);
    
    T=Traw-p;

    [pi_hat,Y,X,~,~,err]=VAR(y,p,c); % ML = OLS estimate 
    
    sigma=err'*err; % SSE
    
    % Draw Q from IW ~ (S,v), where v = T-(K*p+1)-K-1 :
    
    Q=iwishrnd(sigma,T-size(X,2)-K-1); 
    
    % Compute the Kronecker product Q x inv(X'X) and vectorize the matrix
    % pi_hat in order to obtain vec(pi_hat):
    
    XX=kron(Q,inv(X'*X)); 
    s=size(pi_hat)';
    vec_pi_hat=reshape(pi_hat,s(1)*s(2),1);
    
    % Draw PI from a multivariate normal distribution with mean vec(pi_hat)
    % and variance Q x inv(X'X):
    
    PI=mvnrnd1(vec_pi_hat,XX,1);
    PI=PI';
    PI=reshape(PI,[K*p+c,K]); % reshape PI such that Y=X*PI+e, i.e. PI is (K*p+c)x(K).
    
    % Create the companion form representation matrix A:
    
    BigA=[PI(1+c:end,:)'; eye(K*p-K) zeros(K*p-K,K)]; % (K*p)x(K*p) matrix
    
    % Store errors and fitted values:
    
    errornorm=Y-X*PI;
    fittednorm=X*PI;
    
end
    
    
    
    
    


