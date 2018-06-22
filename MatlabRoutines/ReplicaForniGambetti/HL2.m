  
function [q, o_log, o,rr,rs] = HL2(panel,  kmax, cmax, nmin, penalty)
% 

%      cmax : c = [0:cmax], e.g cmax = 3
%
%   penalty :
%     p1 = ((M/T)^0.5 + M^(-2) + n^(-1))*log(min([(T/M)^0.5;  M^2; n]))  
%     p2 = (min([(T/M)^0.5;  M^2; n])).^(-1/2)  
%     p3 = (min([(T/M)^0.5;  M^2; n])).^(-1)*log(min([(T/M)^0.5;  M^2; n]))
%
%--------------------------------------------------------------------------
%
%
%  
%


[T,n] = size(panel);
x = standardize(panel);
[a Ns]=sort(rand(n,1));
x = x(:,Ns);
s=0;


for N = nmin:round((n-nmin)/10):n
s=s+1;
tau = round(N*T/n);
M = 2*round(sqrt(tau));
    eigv = subr_1(x(end-tau+1:end,1:N), M); 
    eigv = (eigv >0).*eigv;
    IC1 = flipud(cumsum(flipud(eigv)));
    IC1 = IC1(1:kmax+1,:);
    [a1,a2] = size(IC1);
    

    p = ((M/tau)^0.5 + M^(-2) + N^(-1))*log(min([(tau/M)^0.5;  M^2; N]))*ones(a1,a2);  

    if penalty == 'p2'
        p = (min([(tau/M)^0.5;  M^2; N])).^(-1/2)*ones(a1,a2);  
    elseif penalty == 'p3'
        p = (min([(tau/M)^0.5;  M^2; N])).^(-1)*log(min([(tau/M)^0.5;  M^2; N]))*ones(a1,a2);  
    end
    
    
    for c = 1:floor(cmax*100)
        cc = c/100;
        
        IC_log = log(IC1./N) + kron([0:kmax]',ones(1,a2)).*p*cc;
        IC = (IC1./N) + kron([0:kmax]',ones(1,a2)).*p*cc;
        
        [rr,rs]=find((IC_log == ones(kmax+1,1)*min(IC_log))==1);
        o_log(s,c)=rr-1;
    
        [rr,rs]=find((IC == ones(kmax+1,1)*min(IC))==1);
        o(s,c)=rr-1;
    end %c
end %n


g1 = find(std(o) == 0);
f1 = find(o(end, g1) < kmax);
if numel(f1) > 0;
q1 = o(end, g1(f1(1)) );
else
    q1=999;
end
g2 = find(std(o_log) == 0);
f2 = find(o_log(end, g2) < kmax);
if numel(f2) > 0
q2 = o_log(end, g2(f2(1)) );
else
    q2=999;
end
q = [q1 q2];

set(0,'DefaultLineLineWidth',2);
cr=[1:floor(cmax*100)]'/100;

% %figure
% %subplot(2,1,1)
% %plot(cr,o(end,:),'r-')
% axis tight
% xlabel('c')
% legend('q^{*T}_{c;n}')
% title('estimated number of factors,  IC_1 -  nolog criterion')
% 
% subplot(2,1,2)
% plot(cr,std(o),'b-')
% xlabel('c')
% axis tight
% legend('S_c')
% title('S_c, IC_1 - nolog criterion')
% 
% 
% figure 
% 
% subplot(2,1,1)
% 
% plot(cr,o_log(end,:),'r-')
% axis tight
% xlabel('c')
% legend('q^{*T}_{c;n}')
% title('estimated number of factors, IC_2 - log criterion')
% 
% subplot(2,1,2)
% plot(cr,std(o_log),'b-')
% xlabel('c')
% axis tight
% legend('S_c')
% title('S_c,  IC_2 - log criterion')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  E = subr_1(x, M)

%x = center(x);
% define some useful quantities
[T,N] = size(x);
w = M;
W = 2*w+1;
% compute covariances
SS = zeros(W,N*N);
for k = 1:w+1,
     Sk = center(x(k:T,:))'*center(x(1:T+1-k,:))/(T-k);
     S_k = Sk';
     % the x(:) creates a column vector from the matrix x columnwise
     SS(w-k+2,:) = S_k(:)';  
     SS(w+k,:) = Sk(:)';
end


% compute the spectral matrix in w points (S)
% [0:2*pi/W:4*pi*w/W]
%Factor = exp(-sqrt(-1)*(-w:w)'*(0:2*pi/W:4*pi*w/W));


freq = [0:2*pi/W:pi];
Factor = exp(-sqrt(-1)*(-w:w)'*freq);

ww = 1 - abs([-w:w])/(w+1);

% multiply first line of SS by ww(1)
% multiply second line of SS by ww(2)

S = diag(ww)*SS(1:W,:); 

% inverse to the (:) operation
% S = reshape(S'*Factor,N,N*W);

S = reshape(S'*Factor,N,N*(w+1));

% compute the eigenvalues 
eigenvalues = zeros(N,w+1);
D = eig(S(:,1:N));
eigenvalues(:,1) = flipud(sort(real(D)));

for j = 1:w,
   D = eig(S(:,j*N+1:(j+1)*N));
   eigenvalues(:,1+j) = flipud(sort(real(D)));
end
%eigenvalues
%[eigenvalues(:,1)  eigenvalues(:,2:jj+1)*2]
E = [eigenvalues(:,1)  eigenvalues(:,2:w+1)*2]*ones(w+1,1)/(2*w+1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XC = center(X)

[T n] = size(X);
XC = X - ones(T,1)*mean(X); 
XC = XC./kron(ones(T,1),std(XC));





