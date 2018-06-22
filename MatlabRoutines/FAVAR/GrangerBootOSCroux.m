function  [pvalMSFE pvalreg ER E] = GrangerBootOSCroux(x, y, nvarint ,R, nboot, k)
if nargin <= 5
      k = aicbic([x y], 4); 
end
if nargin == 4
    nboot = 50;
end
[T, n] = size(x);
m = size(y, 2);
l = n + m;
T = T-k;
ER = zeros(T-R,nvarint);
E = zeros(T-R,nvarint);
for s = R+1:T
    [YY , XX] = VAR_str(x(1:s,:), 1, k);
    Y = YY(1:end-1,:);
    X =  XX(1:end-1,:);
    CR = (X'*X)\X'*Y;
    ERR = YY(end,:) - XX(end,:)*CR;
    ER(s-R,:) = ERR(1,1:nvarint);
    [YY , XX] = VAR_str([x(1:s,:) y(1:s,:)], 1, k);
    Y = YY(1:end-1,:);
    X =  XX(1:end-1,:);
    C = (X'*X)\X'*Y;
    EE = YY(end,:) - XX(end,:)*C;
    E(s-R,:) = EE(1,1:nvarint);
end
MSFEpoint =  (log(det(ER'*ER))-log(det(E'*E)));
D = ER - E;
for ss = 1:nvarint
    eps(:,ss) = ER(:,ss) - D(:,ss)*((D(:,ss)'*D(:,ss))\D(:,ss)'*ER(:,ss));
end
regpoint = (T-R)*(log(det(ER'*ER))-log(det(eps'*eps)));



[Y , X] = VAR_str(x, 1, k);
T = size(X,1);
CR = (X'*X)\X'*Y;
UR = Y - X*CR;
%[Y , X] = VAR_str([x y], 1, k);
%C = (X'*X)\X'*Y;
%U = Y - X*C;
%CC = C;
%CC(1:end, 1:n) = zeros(size(C, 1), n);
%CC(1, 1:n) = CR(1,:);
%for j = 1:k
 %   CC(l*(j-1)+2:l*(j-1)+n+1,1:n) = CR((j-1)*n+2:j*n+1,:);
%end
%U(:,1:n) = UR;
%comp  = [CC(2:end,:)'; [eye(l*(k-1)) zeros(l*(k-1),l)]];

CC = CR;
U = UR;
l = n;
comp  = [CR(2:end,:)'; [eye(n*(k-1)) zeros(n*(k-1),l)]];
W(:,1) = X(1,2:end)';




MSFE = zeros(nboot,1);
reg = MSFE;
for i = 1:nboot
       for t = 2:T+1
        W(:,t) = [CC(1,:)'; zeros(l*(k-1),1)] + comp*W(:,t-1) + [U(ceil(rand*T),:)'; zeros(l*(k-1),1)];
    end
    x = W(1:n,2:T+1)';
   % y = W(n+1:l,2:T+1)';
    
  
    BER = zeros(T-R,nvarint);
    BE = zeros(T-R,nvarint);
    for s = R+1:T
        [YY , XX] = VAR_str(x(1:s,:), 1, k);
        Y = YY(1:end-1,:);
        X =  XX(1:end-1,:);
        CR = (X'*X)\X'*Y;
        BERR = YY(end,:) - XX(end,:)*CR;
        BER(s-R,:) = BERR(1,1:nvarint);
        [YY , XX] = VAR_str([x(1:s,:) y(1:s,:)], 1, k);
        Y = YY(1:end-1,:);
        X =  XX(1:end-1,:);
        C = (X'*X)\X'*Y;
        BEE = YY(end,:) - XX(end,:)*C;
        BE(s-R,:) = BEE(1,1:nvarint);
    end
    MSFE(i) =  (log(det(BER'*BER))-log(det(BE'*BE)));
    D = BER - BE;
    for ss = 1:nvarint
        eps(:,ss) = BER(:,ss) - D(:,ss)*((D(:,ss)'*D(:,ss))\D(:,ss)'*BER(:,ss));
    end
    reg(i) = (T-R)*(log(det(BER'*BER))-log(det(eps'*eps)));
    
    
    
end
pvalMSFE = sum(MSFE>MSFEpoint)/nboot;
pvalreg = sum(reg>regpoint)/nboot;
