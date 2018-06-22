function y=DiscreteDraw(p)
n=length(p);
a=rand(1,1);
[b c]=sort(p);
pr=cumsum(p(c));
j=1;
while j<=n
    r=(a<pr(j));
    if r==1
        y=c(j);
        break
    else
        j=j+1;
    end
end
    