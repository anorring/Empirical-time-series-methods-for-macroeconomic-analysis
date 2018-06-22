function [v, stdv, meanv] = ExplainedVariances(bootcumimpulse,cumimpulse,  shock)
a = cumsum(bootcumimpulse.^2,3);
v = squeeze(a(:,shock,:,:))./squeeze(sum(a,2));
clear a;
meanv = mean(v,3);
stdv = std(v,0,3);
a = cumsum(cumimpulse.^2,3);
v = squeeze(a(:,shock,:))./squeeze(sum(a,2));