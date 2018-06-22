function [w, u, coeffmatrix3D] = woldimpulse(y,nlagsvar,nlagsimp,c)
[T,N] = size(y);
if nargin == 3,
   c = 1;
end
if nargin == 2,
   c = 1;
   nlagsimp = 30;
end
[varoutput, u] = myvar(y,nlagsvar,c);
coeffmatrix = varoutput(1:2:2*N,c+1:N*nlagsvar+c);
coeffmatrix3D = zeros(N,N,nlagsvar+1);
coeffmatrix3D(:,:,1) = eye(N);
for k = 1:nlagsvar,
   coeffmatrix3D(:,:,k+1) = -coeffmatrix(:,k:nlagsvar:N*nlagsvar);
end;
w = invertepolynomialmatrix(coeffmatrix3D,nlagsimp);

