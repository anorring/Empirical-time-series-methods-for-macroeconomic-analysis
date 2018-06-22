function [B, chi, rsh, LR,sh]  = FAVARLR(X, Z, k, h)
T = size(X, 1);
W = [ones(T,1) Z];
AA = inv(W'*W)*W'*X;
chi = W*AA;
A = AA(2:end,:);
[BB, epsilon, LR] = woldimpulseLR(Z, k, h + 1);
Sigma = cov(epsilon);
% C = chol(Sigma)';
%LR=A'*LR;
%LR=LR(intv,:);
H = inv(LR)*chol(LR*Sigma*LR')';
check_sign=LR*H;
if check_sign(1,1)<0
    H(:,1)=-H(:,1);
end
for lag = 1 : h + 1
    B(:, :, lag) = A'*BB(:, :, lag)*H;
end
rsh = epsilon/H;
sh = inv(H)*rsh';


