function [canada,japan,uk] = uhligcumulatedforwardpremium(irf,opt)
canada = cumsum(squeeze(irf(84,3,2:end,1:3))*1200+squeeze(irf(110,3,1:end-1,1:3)));
japan = cumsum(squeeze(irf(85,3,2:end,1:3))*1200+squeeze(irf(111,3,1:end-1,1:3)));
uk = cumsum(-squeeze(irf(87,3,2:end,1:3))*1200+squeeze(irf(112,3,1:end-1,1:3)));
if opt == 1
figure
plot(1:size(irf,3)-1,canada(:,[1 3]),'k:',  'LineWidth',.5);
       hold on;
    plot(1:size(irf,3)-1,canada(:,2),'k', 'LineWidth',2);
    figure
plot(1:size(irf,3)-1,japan(:,[1 3]),'k:',  'LineWidth',.5);
       hold on;
    plot(1:size(irf,3)-1,japan(:,2),'k', 'LineWidth',2);
    figure
plot(1:size(irf,3)-1,uk(:,[1 3]),'k:',  'LineWidth',.5);
       hold on;
    plot(1:size(irf,3)-1,uk(:,2),'k', 'LineWidth',2);
end