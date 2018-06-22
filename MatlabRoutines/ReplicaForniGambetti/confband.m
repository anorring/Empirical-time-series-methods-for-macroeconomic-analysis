function C = confband(BU, B, level)
%BB = cumsum(B,3);
[ N, q, h, nrepli] = size(BU);
m = round((1-level)*nrepli/2);
for k = 1:q
    for j=1:N
CU = sort(squeeze(BU(j,k,:,:)),2);
Clow(j,k,:) = CU(:,m  );
Cup(j,k,:) = CU(:, nrepli-m);
Cmed(j,k,:) = median(CU,2);
    end
end
C(:,:,:,1) = Clow;
C(:,:,:,2) = B;
C(:,:,:,3) = Cup;
C(:,:,:,4) = Cmed;



