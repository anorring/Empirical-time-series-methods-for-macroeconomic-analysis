function [B, chi, rsh]  = FAVARRaw(X, Z, k, h)
N = size(X, 2);

T = size(X, 1);
W = [ones(T,1) Z];
AA = (W'*W)\W'*X;
chi = W*AA;
A = AA(2:end,:);
[BB, epsilon, coeff] = woldimpulse(Z, k, h + 1);
Sigma = cov(epsilon);
C = chol(Sigma)';
for lag = 1 : h + 1
    B(:, :, lag) = A'*BB(:, :, lag)*C;
end
rsh = epsilon/C;


