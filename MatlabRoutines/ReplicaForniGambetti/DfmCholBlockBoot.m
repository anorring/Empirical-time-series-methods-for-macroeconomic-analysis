%
% B  = DfmCholBlockBoot(X, q, r,  h, L, nrepli)
%
%
% Le nouve realizzazioni di X sono ottenute con block bootstrapping, L
% lunghezza dei blocchi in no. periodi
% 
%
function B = DfmCholBlockBoot(X, variables, r, k,h, L, nrepli)
q = length(variables);
N = size(X, 2);
B = zeros( N, q, h + 1,nrepli);
for j = 1 : nrepli
  X_boot = boot_block(X,L);


BB( :, :, :, j) = DfmRawImp(X_boot, q, r,k, h);
B( :, :, :,j) = DfmCholIdent(BB( :, :, :,j),variables);
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


   