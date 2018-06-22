function [pi_hat,Y,X,Y_initial,Yfit,err]=VAR(data,p,c)

% Function to estimate VAR by OLS, using the companion form representation
% Author: Nicolo' Maffei Faccioli
% INPUTS: data = data at hand, p = # of lags, c=1 if constant, c=0 if not.
% OUTPUTS: pi_hat = estimated coefficients organized in a n*pxn matrix
% (including the constant if c==1, not including it if c==0) ; Y = matrix
% of dependent variables ; X = matrix of regressors (lags + constant if
% c=1) ; Y_initial = initial p elements discarded from data.


[Y,X,Y_initial]=SUR(data,p,c);
pi_hat=(X'*X)\X'*Y;
Yfit=X*pi_hat; %Fitted value of Y
err=Y-Yfit; %Residuals

end