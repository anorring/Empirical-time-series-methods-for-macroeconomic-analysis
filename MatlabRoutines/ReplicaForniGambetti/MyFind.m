function y=MyFind(x);
y=[];
for i=1:size(x)
    if x(i)~=0
        y=[y; x(i)];
    end
end