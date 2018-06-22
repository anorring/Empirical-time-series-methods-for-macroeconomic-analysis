function  p = ortotest(u, reg, lags)
z = size(u,1);
rog = [ones(z-lags,1)];
for j=1:lags
    rog = [rog reg(lags+1-j:end-j,:) ];
end
[a b]=myols(u(lags+1:end),rog);
R=a(1,end);
dn = lags*size(reg,2);  
dd=z-dn-1;  
p =1-fcdf((R/dn)/((1-R)/dd),dn,dd);
