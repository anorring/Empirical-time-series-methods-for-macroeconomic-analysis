%
% B = DfmRawImp(X, q, r,k, h)
% 
% produce la stima come descritto in Forni Lippi Reichlin 2003 
% delle impulse (non cumulate) relative al
% panel X con rango statico r e dinamico q e h ritardi
% B e' tridimensionale con  la cross-sectional unit sulla prima, 
% il common shock sulla seconda,
% il ritardo sulla terza dimensione. Le impulse in B corrispondono ad una 
% identificazione ortonormale, quale risulta determinata 
% dalla scelta degli autovettori operata dalla routine  eigs.
%
% 


function [B, Chi, rsh]  = DfmRawImp(X, q, r, k, h,weights)
N = size(X, 2);
if nargin == 5
    weights = ones(1,N);
end
T = size(X, 1);
WW = diag(std(X)./weights);
x = center(X)*(WW^-1);
Gamma0 = cov(x);
opt.disp = 0;
[W, Lambda] = eigs(Gamma0, r,'LM',opt);
F = x*W;
chi = F*W';
Chi = chi*WW + ones(T,1)*mean(X);
[BB, epsilon, coeff] = woldimpulse(F, k, h + 1);
Sigma = cov(epsilon);
[ K, MM ] = eigs(Sigma, q, 'LM',opt);
M = diag(sqrt(diag(MM)));
for lag = 1 : h + 1
    B(:, :, lag) = WW*W*BB(:, :, lag)*K*M;
end
rsh = epsilon*K*inv(M);


