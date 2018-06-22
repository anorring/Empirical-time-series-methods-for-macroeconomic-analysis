%%% QR factorization with non-negative diagonals in R:

function [Q,R] = qr_dec(drawW)

[Q,R] = qr(drawW);

for i = 1:size(drawW,1)
    
    if R(i,i)<0
        R(i,:)=-R(i,:);
        Q(:,i)=-Q(:,i);
        
    end
end