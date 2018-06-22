function DoSubPlots(imp,d1,d2,ind,sho,titoli)
ii=0;
for i=1:d1
    for j=1:d2
        ii=ii+1;
        subplot(d1,d2,ii),plot(1:size(imp,3),squeeze(imp(ind(ii),sho,:,[1 3])),'k:',  'LineWidth',1);
        if nargin==6
            title(titoli(ind(ii),:));
        end
        hold on;
        plot(1:size(imp,3),squeeze(imp(ind(ii),sho,:,2)),'k', 'LineWidth',1.5);
        axis tight;hold off
        
%         subplot(d1,d2,ii),plot(1:size(imp,3),squeeze(imp(ind(ii),sho,:,[1 3])),'k:', ...
%             1:size(imp,3),squeeze(imp(ind(ii),sho,:,2)),'k', 'LineWidth',1);
%         axis tight
        
    end
end
