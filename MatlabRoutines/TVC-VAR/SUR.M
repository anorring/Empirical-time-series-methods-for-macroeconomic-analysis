% sur.m
function [theta,Vtheta,Vu] = sur(Y,X,T);

% function [theta,Vtheta,Vu] = sur(Y,X,T);

% offline 2-step SUR 

% Y is a N by T matrix, with the row indicating equations, and column 
% indicating time periods, i.e, column j is the observation of Y 
% at time j.

% X is a tensor object, with dimension(K_1+k_2+...K_N) by N by T, 
% where K_i is the number of independent variables for equation i,
% and T indicates time as above. For page t, we have (cap)X_t' as 
% in equation (1.2) in Notes. So page t records observations of right
% hand side variables for N equations at time t.

%  initialize moment matrices
[ry,cy] = size(Y);
[rx,cx,px] = size(X);

Mxx = zeros(rx,rx);
Mxy = zeros(rx,1);

% 1st stage estimates; weighting matrix = I
for t = 1:T,
   Mxx = Mxx + X(:,:,t)*X(:,:,t)';
   Mxy = Mxy + X(:,:,t)*Y(:,t);
end
theta = inv(Mxx)*Mxy;

% 1st stage residuals 
e = zeros(ry,T);
for t = 1:T,
   e(:,t) = Y(:,t) - X(:,:,t)'*theta;
end

% 2nd stage estimates; weighting matrix = inv(cov(e))
W = inv(cov(e'));
Mxx = zeros(rx,rx);
Mxy = zeros(rx,1);
for t = 1:T,
   Mxx = Mxx + X(:,:,t)*W*X(:,:,t)';
   Mxy = Mxy + X(:,:,t)*W*Y(:,t);
end
theta = inv(Mxx)*Mxy; 

% 2nd stage residuals
for t = 1:T,
   e(:,t) = Y(:,t) - X(:,:,t)'*theta;
end

Vu = cov(e');
W = inv(Vu);

Mxx = zeros(rx,rx);
for t = 1:T,
   Mxx = Mxx + X(:,:,t)*W*X(:,:,t)';
end
   
Vtheta = inv(Mxx);