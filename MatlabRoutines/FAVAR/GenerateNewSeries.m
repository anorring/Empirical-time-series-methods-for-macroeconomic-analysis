function databoot = GenerateNewSeries(B,C,X,u,p);
T=size(X,1);
n=size(u,1);
for t=1:T
    if t==1
        YB(:,1)=X(1,2:end)';
    else
        a=fix(rand(1,1)*size(u,2))+1;
        uu=u(:,a);
        YB(:,t)=[C;zeros(n*p-n,1)]+B*YB(:,t-1)+[uu;zeros(n*p-n,1)];
    end
      
end
databoot=YB(1:n,:)';

