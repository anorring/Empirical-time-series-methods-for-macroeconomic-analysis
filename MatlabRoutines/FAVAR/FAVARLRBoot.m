%
%
%
%
function B = FAVARLRBoot(Data, Z,k,h,nrepli,block,L)
[T N] = size(Data);
r = size(Z, 2);
B = zeros( N, r, h + 1,nrepli);

[VarPa C X u] = VarParameters(Z,k,1);

W = [ones(T,1) Z];
AA = inv(W'*W)*W'*Data;
chi = W*AA;

Idio = Data - chi;

for j=1:nrepli
    Z_boot = GenerateNewSeries(VarPa,C,X,u,k);
    W_boot = [ones(T-k,1) Z_boot];
    if block==1
        t1=size(W_boot,1)-floor(size(Idio(k+1:end,:),1)/L)*L+1;
        chi_boot=W_boot*AA;
        X_boot=chi_boot(t1:end,:)+ boot_block(Idio(k+1:end,:),L);
    else
        X_boot = W_boot*AA +Idio(k+1:end,:);
    end
    B( :, :, :, j)  = FAVARLR(X_boot, Z_boot, k, h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_boot = boot_block(X,L);

OPTS.disp=0;

[T,N] = size(X);

K = floor(T/L);

blks = ceil(rand(K,1)*K);
for i = 1:K
    X_boot(((i-1)*L+1):(i*L),:) = X(((blks(i)-1)*L+1):(blks(i)*L),:);
end;


