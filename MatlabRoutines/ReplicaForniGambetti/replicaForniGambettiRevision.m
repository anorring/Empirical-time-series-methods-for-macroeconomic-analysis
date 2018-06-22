load GF
load codeGF
load transGF
codeGF(83:87) = 0;
codeGF(71) = 1;
x = standardize(GF);
% the number of static and dynamic factors
% static factors (Bai Ng 2002)
[PC,IC] = baingcriterion(x,25);
[minimum,rhat] = min(IC(:,2));
Fhat = principalcomponents(x,rhat);
% number of lags in the VAR
%[AIC, BIC, e] = aicbic(Fhat,9);
% dynamic factors (Bai Ng 2007, SW 2005)
k = 2;
nrepli=500;
baingqhat = baingdftest(x,rhat,k);
swqhat = swdftest(x,rhat,k);
% dynamic factors (Hallin Liska 2007)
%q1 = HL2(x,  16,  1.5, 30, 'p1');
%q2 = HL2(x,  16,  1.5, 30, 'p2');
%q3 = HL2(x,  16,  1.5, 30, 'p3');
%HLqhat = [q1;q2;q3];
% no of dyn fractors with r static factors according to Bai Ng 2007
rhat
%HLqhat
baingqhat
swqhat
intv = [5 96 75 106];
iv = [107:109 105 89 67 71 3 72 65 54 62 19 44 47 22 20 25 30 ];

%%%% construct labels
for i=1:size(GF,2)
    titoli(i,:)=['                            '];
end
titoli(intv,:)= ...
   ['      Ind. production       ';
    '           CPI              ';
    '    Federal funds rate      ';
    '       Swi/US real ER       '];

titoli(iv,:)= ...
   ['       Jap/US real ER       ';
    '        UK/US real ER       ';
    '       Can/US real ER       ';
    '          Earnings          ';
    '            PPI             ';
    '             M2             ';
    '           Loans            ';
    '   Real pers. consumption   ';
    '       Consumer credit      ';
    '          Orders            ';
    '       Housing starts       ';
    '        Inventories         ';
    '    Capacity utilization    ';
    '           Hours            ';
    '     Hours manufacturing    ';
    '         Employment         ';
    '         Vacancies          ';
    '        Unemployment        ';
    '     Unemployment rate      '];

[virf vbootirf] = CholeskyBoot(GF(:,intv),9,49,1,nrepli,.8,[1 2 ],0);
for i=1:size(vbootirf,4),
    vbootirf(:,:,:,i)=vbootirf(:,:,:,i)/vbootirf(3,3,1,i)*0.5;
end

% vbootirf=vbootirf/virf(3,3,1)*0.5;

virf=virf/virf(3,3,1)*0.5;

varirf = confband(vbootirf,virf,.8);

varirf([1 2 4],:,:,:)=varirf([1 2 4],:,:,:)*100;

[irf, v, stdv, meanv,chi,sh,ci,bootci,imp,bootimp] =DoEverything(GF, intv, rhat, k, 52, codeGF,nrepli);

for i=1:size(bootci,4),
    bootci(:,:,:,i)=bootci(:,:,:,i)/bootci(75,3,1,i)*0.5;
end

% bootci=bootci/ci(75,3,1)*0.5;

ci=ci/ci(75,3,1)*0.5;

irf = confband(bootci, ci, .8);


[canada,japan,uk] = uhligcumulatedforwardpremium(irf,0);
commonvar = var(chi)./var(GF);

index5=MyFind([1:112]'.*(transGF==5)')
index4=MyFind([1:112]'.*(transGF==4)')
index=sort([index5;index4]);

irf(index,:,:,:)=irf(index,:,:,:)*100; 


% figure(1);
% DoSubPlots(varirf,2,2,[1 2 3 4],3,titoli(intv,:));
 
[vvar, stdvar] = ExplainedVariances(vbootirf,virf, 3);
expvariancevar = [vvar(:,1) stdvar(:,1) vvar(:,7) stdvar(:,7) vvar(:,13) stdvar(:,13) vvar(:,49) stdvar(:,49)]

% figure(2);
% DoSubPlots(irf,2,2,intv,3,titoli);

figure(1);
DoSubPlots([varirf;irf],4,2,[1 9 2 100 3 79 4 110],3);


expvariancedfm = [v(intv,1) stdv(intv,1) v(intv,6) stdv(intv,6) v(intv,12) stdv(intv,12) v(intv,48) stdv(intv,48)]
expvariancedfm2 = [v(iv,1) stdv(iv,1) v(iv,6) stdv(iv,6) v(iv,12) stdv(iv,12) v(iv,48) stdv(iv,48)];
%expvariance = [expvariancevar; expvariancedfm]
figure(2);
subplot(3,2,2),plot(1:size(irf,3)-1,canada(:,[1 3]),'k:',  'LineWidth',.5);
hold on;
plot(1:size(irf,3)-1,canada(:,2),'k', 'LineWidth',2);
title(['Canada/US UIP']);
axis tight; hold off;
subplot(3,2,4), plot(1:size(irf,3)-1,japan(:,[1 3]),'k:',  'LineWidth',.5);
hold on;
plot(1:size(irf,3)-1,japan(:,2),'k', 'LineWidth',2);
title(['Japan/US UIP']);
axis tight; hold off;
subplot(3,2,6), plot(1:size(irf,3)-1,uk(:,[1 3]),'k:',  'LineWidth',.5);
hold on;
plot(1:size(irf,3)-1,uk(:,2),'k', 'LineWidth',2);
title(['UK/US UIP']);
axis tight; hold off
subplot(3,2,1),plot(1:size(irf,3),squeeze(irf(109,3,:,[1 3])),'k:',  'LineWidth',1);
title(['Canada/US real exchange rate']);
hold on;
plot(1:size(irf,3),squeeze(irf(109,3,:,2)),'k', 'LineWidth',1.5);
axis tight; hold off
subplot(3,2,3),plot(1:size(irf,3),squeeze(irf(107,3,:,[1 3])),'k:',  'LineWidth',1);
title(['Japan/US real exchange rate']);
hold on;
plot(1:size(irf,3),squeeze(irf(107,3,:,2)),'k', 'LineWidth',1.5);
axis tight; hold off
subplot(3,2,5),plot(1:size(irf,3),-squeeze(irf(108,3,:,[1 3])),'k:',  'LineWidth',1);
title(['UK/US real exchange rate']);
hold on;
plot(1:size(irf,3),-squeeze(irf(108,3,:,2)),'k', 'LineWidth',1.5);
axis tight; hold off


%%%%%% old figures
% figure(4);
% DoSubPlots(irf,2,2,[105 89 67 71],3,titoli);
% figure(5);
% DoSubPlots(irf,3,2,[3 72 65 54 62 19],3,titoli);
% figure(6);
% DoSubPlots(irf,3,2,[44 47 22 20 25 30  ],3,titoli);

figure(3);
DoSubPlots(irf,5,3,[105 89 67 71 3 72 65 54 62 44 47 22 20 25 30],3,titoli);

nsfactors = [4 10 16 ];
nlags = [8 3 2];
for j = 1:3
    [imp3, chi3] = DfmCholImp(GF, intv, nsfactors(j),nlags(j), 48);
    cimp3 = CumImp(imp3, codeGF);
    cimp3(108,3,:) = -cimp3(108,3,:);
    cimp3cr(:,j,:) = squeeze(cimp3([96 106:109],3,:))*100;
    cimp3cr(:,j,:) = cimp3cr(:,j,:)/cimp3(75,3,1)*0.5;
end

figure(4)
for j=1:4
    for k=1:3
        subplot(4,3,3*(j-1)+k)
        plot(squeeze(cimp3cr(j+1,k,:))','k'),
        axis tight
        if j==1 & k==1
            ylabel(titoli(106,:))
            title('4 static factors')
        end
        if j==1 & k==2
            title('10 static factors')
        end
        if j==1 & k==3
            title('16 static factors')
        end
        if j==2 & k==1
            ylabel(titoli(107,:))
        end
        if j==3 & k==1
            ylabel(titoli(108,:))
        end
        if j==4 & k==1
            ylabel(titoli(109,:))
        end
    end
end
    
% figures for appendix
[irf, v, stdv, meanv,chi,sh,ci,bootci,imp,bootimp] =DoEverything(GF, intv, rhat, k, 52, codeGF,10);

ci=ci/ci(75,3,1)*0.5;
ci(index,:,:)=ci(index,:,:)*100;
ci(108,:,:)=-ci(108,:,:);

[irf2, v2, stdv2, meanv2, chi2, sh2, ci2, bootci2] =DoEverything(GF, [intv 107:109], rhat, 2, 52, codeGF,10);

ci2=ci2/ci2(75,3,1)*0.5;
ci2(index,:,:)=ci2(index,:,:)*100;
ci2(108,:,:)=-ci2(108,:,:);

[irf3, v3, stdv3, meanv3, chi3, sh3, ci3, bootci3] =DoEverything(GF, [intv ], 10, 2, 52, codeGF,10);

ci3=ci3/ci3(75,3,1)*0.5;
ci3(index,:,:)=ci3(index,:,:)*100;
ci3(108,:,:)=-ci3(108,:,:);

figure(5);
idx=[intv 107:109 54 65 19  22 30];
for i=1:12
    subplot(4,3,i),plot(1:49,squeeze(ci(idx(i),3,:)),'k',...
                        1:49,squeeze(ci2(idx(i),3,:)),':k',...
                        1:49,squeeze(ci3(idx(i),3,:)),'--k'),title(titoli(idx(i),:))
end
%     ,'LineWidth',1.5
% for j = 1: 16
%     [imp4, chi4] = DfmCholImp(GF, intv,j+3,ceil(30/(j+3)), 48);
%     cimp4 = CumImp(imp4, codeGF);
%     a = cumsum(cimp4.^2,3);
%     v4 = squeeze(a(:,3,:))./squeeze(sum(a,2));
%     t4(:,j) = v4([intv 107:109],[ 48]);
% end
