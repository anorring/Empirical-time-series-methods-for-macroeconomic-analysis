%
%
% 
%
function B = FAVARLRBlockBoot(X, Z,k,h, L, nrepli)
N = size(X, 2);
r = size(Z, 2);
B = zeros( N, r, h + 1,nrepli);

for j = 1 : nrepli
  XZ_boot = boot_block([X Z],L);
X_boot = XZ_boot(:,1:N);
Z_boot =  XZ_boot(:,N+1:end);

B( :, :, :, j)  = FAVARLR(X_boot, Z_boot, k, h);
    
%B( :, :, :,j)  = DfmLRIdent(BB( :, :, :, j), variables,LR,rsh);

% j
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


   