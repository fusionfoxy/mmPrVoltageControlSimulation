clear all; close all; clc;

filePrefix = '5sec_0985lim';%NewMM
folder = '5sec_plots';
delays = [0:0.05:0.2 0.3:0.1:0.8 1.3 1.6]; %[0:0.1:1 1.5:0.5:10]; %  [0:0.2:1 2:1:6]; %% 
seeds = [1:17 31:41 46]; %1:60; %

alpha = 0.05;
reactive = false;
epsilon = 2;
maxSeed = length(seeds);
%%
% losses = 0:10:90;
% data = zeros(1,length(losses));
res1 = zeros(maxSeed,length(delays));
mm1 = res1;
res0 = zeros(maxSeed,length(delays));
mm0 = res0;
% data = zeros(1,length(delays));
for id =1:length(delays)
    disp(['delay ' num2str(id) ' of ' num2str(length(delays)) ' with ' num2str(maxSeed) ' seeds']);
    for seed = seeds
%         disp([id seed])
        data = load(['results/' filePrefix '_delay' num2str(delays(id)) '_reactive1_seed' num2str(seed) '.mat']);%res/loss/
        tmp = abs(data.vOut(1:data.N,data.lvgc_assetBusses))/data.LV_vBase;
        res1(seed,id) = length(tmp(tmp<0.98));
        tmp = max(data.mm');
        mm1(seed,id) = length(tmp(tmp>epsilon))/length(tmp);
        
        data = load(['results/' filePrefix '_delay' num2str(delays(id)) '_reactive0_seed' num2str(seed) '.mat']);%res/loss/
        tmp = abs(data.vOut(1:data.N,data.lvgc_assetBusses))/data.LV_vBase;
        res0(seed,id) = length(tmp(tmp<0.98));
        tmp = max(data.mm');
        mm0(seed,id) = length(tmp(tmp>epsilon))/length(tmp);
    end

end
data = load(['results/' filePrefix '_delay' num2str(delays(1)) '_reactive' num2str(reactive) '_seed' num2str(1) '.mat']);%res/loss/
for id = length(delays):-1:1
    confRes1(1,id)=nanmean(res1(:,id))+ 2.576*nanstd(res1(:,id))/sqrt(length(res1(~isnan(res1(:,id)),id)));
    confRes1(2,id)=nanmean(res1(:,id))- 2.576*nanstd(res1(:,id))/sqrt(length(res1(~isnan(res1(:,id)),id)));
    confmm1(1,id)=nanmean(mm1(:,id))+ 2.576*nanstd(mm1(:,id))/sqrt(length(mm1(~isnan(mm1(:,id)),id)));
    confmm1(2,id)=nanmean(mm1(:,id))- 2.576*nanstd(mm1(:,id))/sqrt(length(mm1(~isnan(mm1(:,id)),id)));
    confRes0(1,id)=nanmean(res0(:,id))+ 2.576*nanstd(res0(:,id))/sqrt(length(res0(~isnan(res0(:,id)),id)));
    confRes0(2,id)=nanmean(res0(:,id))- 2.576*nanstd(res0(:,id))/sqrt(length(res0(~isnan(res0(:,id)),id)));
    confmm0(1,id)=nanmean(mm0(:,id))+ 2.576*nanstd(mm0(:,id))/sqrt(length(mm0(~isnan(mm0(:,id)),id)));
    confmm0(2,id)=nanmean(mm0(:,id))- 2.576*nanstd(mm0(:,id))/sqrt(length(mm0(~isnan(mm0(:,id)),id)));


%     confRes1(1,id)=sqrt(maxSeed/chi2inv(1-alpha/2,maxSeed))*mean(res1(:,id));
%     confRes1(2,id)=sqrt(maxSeed/chi2inv(alpha/2,maxSeed))*mean(res1(:,id));
%     confmm1(1,id)=sqrt(maxSeed/chi2inv(1-alpha/2,maxSeed))*nanmean(mm1(:,id));
%     confmm1(2,id)=sqrt(maxSeed/chi2inv(alpha/2,maxSeed))*nanmean(mm1(:,id)); 
%     confRes0(1,id)=sqrt(maxSeed/chi2inv(1-alpha/2,maxSeed))*mean(res0(:,id));
%     confRes0(2,id)=sqrt(maxSeed/chi2inv(alpha/2,maxSeed))*mean(res0(:,id));
%     confmm0(1,id)=sqrt(maxSeed/chi2inv(1-alpha/2,maxSeed))*nanmean(mm0(:,id));
%     confmm0(2,id)=sqrt(maxSeed/chi2inv(alpha/2,maxSeed))*nanmean(mm0(:,id)); 
end
filePrefix = ['NewMM_0985lim_eps' num2str(epsilon)];
%%
% Voltages
fig = figure(1);
plot(data(end).tvec,abs(data(end).vOut(1:data(end).N,data(end).lvgc_assetBusses))/data(end).LV_vBase)
ylabel('Voltage [p.u.]')
title('LV Grid')
xlabel('Time [hrs]')
print(figure(1),'-depsc',[folder '/voltages_' filePrefix '_reactive0.eps']);
savefig(figure(1),[folder '/voltages_' filePrefix '_reactive' num2str(reactive) '.fig'])
fig.Position = [0 0+570+30 570 510];

fig = figure(2);
plot(delays*data(end).Ts_AssetData,res1)
xlabel('Delay [s]')
ylabel('number of voltage violations pr day [.]');
title('Control Quality with reactive access')
print(figure(2),'-depsc',[folder '/delayVsCtrl_' filePrefix '_reactive1.eps']);
savefig(figure(2),[folder '/delayVsCtrl_' filePrefix '_reactive1.fig'])
fig.Position = [0+510+50 0+570+30 570 510];


fig = figure(3);
plot(res1,mm1,'*')
ylabel('mmPr [.]')
xlabel('number of voltage violations pr day [.]');
title('mmPr vs Control Quality with reactive access')
print(figure(3),'-depsc',[folder '/mmPrVsCtrl_' filePrefix '_reactive1.eps']);
savefig(figure(3),[folder '/mmPrVsCtrl_' filePrefix '_reactive1.fig'])
fig.Position = [0 0 570 510];

fig = figure(4);
plot(delays*data(end).Ts_AssetData,mm1,'*')
ylabel('mmPr [.]')
xlabel('Delay [s]')
title('mmPr with reactive access')
print(figure(4),'-depsc',[folder '/mmPrVsDelay_' filePrefix '_reactive1.eps']);
savefig(figure(4),[folder '/mmPrVsDelay_' filePrefix '_reactive1.fig'])
fig.Position = [0+510+50 0 570 510];


fig = figure(5);
plot(delays*data(end).Ts_AssetData,mean(res1),delays*data(end).Ts_AssetData,confRes1(1,:),'--',delays*data(end).Ts_AssetData,confRes1(2,:),'--')
xlabel('Delay [s]')
ylabel('mean number of voltage violations pr day [.]');
title('Control Quality with reactive access')
print(figure(5),'-depsc',[folder '/delayVsMeanCtrl_' filePrefix '_reactive1.eps']);
savefig(figure(5),[folder '/delayVsMeanCtrl_' filePrefix '_reactive1.fig'])
fig.Position = [0+510+50 0+570+30 570 510];


fig = figure(6);
plot(mean(res1),nanmean(mm1),'*')
ylabel('mean mmPr [.]')
xlabel('mean number of voltage violations pr day [.]');
title('mmPr vs Control Quality with reactive access')
print(figure(6),'-depsc',[folder '/meanMmPrVsMeanCtrl_' filePrefix '_reactive1.eps']);
savefig(figure(6),[folder '/meanMmPrVsMeanCtrl_' filePrefix '_reactive1.fig'])
fig.Position = [0 0 570 510];


fig = figure(7);
plot(delays*data(end).Ts_AssetData,nanmean(mm1),'-',delays*data(end).Ts_AssetData,confmm1(1,:),'--',delays*data(end).Ts_AssetData,confmm1(2,:),'--')
ylabel('mean mmPr [.]')
xlabel('Delay [s]')
title('mmPr with reactive access')
print(figure(7),'-depsc',[folder '/meanMmPrVsDelay_' filePrefix '_reactive1.eps']);
savefig(figure(7),[folder '/meanMmPrVsDelay_' filePrefix '_reactive1.fig'])
fig.Position = [0+510+50 0 570 510];



%%
fig = figure(8);
plot(delays*data(end).Ts_AssetData,res0)
xlabel('Delay [s]')
ylabel('number of voltage violations pr day [.]');
title('Control Quality with periodical access')
print(figure(8),'-depsc',[folder '/delayVsCtrl_' filePrefix '_reactive0.eps']);
savefig(figure(8),[folder '/delayVsCtrl_' filePrefix '_reactive0.fig'])
fig.Position = [0+510+50 0+570+30 570 510];


fig = figure(9);
plot(res0,mm0,'*')
ylabel('mmPr [.]')
xlabel('number of voltage violations pr day [.]');
title('mmPr vs Control Quality with periodical access')
print(figure(9),'-depsc',[folder '/mmPrVsCtrl_' filePrefix '_reactive0.eps']);
savefig(figure(9),[folder '/mmPrVsCtrl_' filePrefix '_reactive0.fig'])
fig.Position = [0 0 570 510];

fig = figure(10);
plot(delays*data(end).Ts_AssetData,mm0,'*')
ylabel('mmPr [.]')
xlabel('Delay [s]')
title('mmPr with periodical access')
print(figure(10),'-depsc',[folder '/mmPrVsDelay_' filePrefix '_reactive0.eps']);
savefig(figure(10),[folder '/mmPrVsDelay_' filePrefix '_reactive0.fig'])
fig.Position = [0+510+50 0 570 510];


fig = figure(11);
plot(delays*data(end).Ts_AssetData,mean(res0),delays*data(end).Ts_AssetData,confRes0(1,:),'--',delays*data(end).Ts_AssetData,confRes0(2,:),'--')
xlabel('Delay [s]')
ylabel('mean number of voltage violations pr day [.]');
title('Control Quality with periodical access')
print(figure(11),'-depsc',[folder '/delayVsMeanCtrl_' filePrefix '_reactive0.eps']);
savefig(figure(11),[folder '/delayVsMeanCtrl_' filePrefix '_reactive0.fig'])
fig.Position = [0+510+50 0+570+30 570 510];


fig = figure(12);
plot(mean(res0),nanmean(mm0),'*')
ylabel('mean mmPr [.]')
xlabel('mean number of voltage violations pr day [.]');
title('mmPr vs Control Quality with periodical access')
print(figure(12),'-depsc',[folder '/meanMmPrVsMeanCtrl_' filePrefix '_reactive0.eps']);
savefig(figure(12),[folder '/meanMmPrVsMeanCtrl_' filePrefix '_reactive0.fig'])
fig.Position = [0 0 570 510];


fig = figure(13);
plot(delays*data(end).Ts_AssetData,nanmean(mm0),'-',delays*data(end).Ts_AssetData,confmm0(1,:),'--',delays*data(end).Ts_AssetData,confmm0(2,:),'--')
ylabel('mean mmPr [.]')
xlabel('Delay [s]')
title('mmPr with periodical access')
print(figure(13),'-depsc',[folder '/meanMmPrVsDelay_' filePrefix '_reactive0.eps']);
savefig(figure(13),[folder '/meanMmPrVsDelay_' filePrefix '_reactive0.fig'])
fig.Position = [0+510+50 0 570 510];


fig = figure(14);
plot(mean(res0),nanmean(mm0),'*',mean(res1),nanmean(mm1),'*')
ylabel('mean mmPr [.]')
xlabel('mean number of voltage violations pr day [.]');
title('mmPr vs Control Quality with periodical access')
legend('Periodical','Reactive');
print(figure(14),'-depsc',[folder '/meanMmPrVsMeanCtrl_' filePrefix '_reactive01.eps']);
savefig(figure(14),[folder '/meanMmPrVsMeanCtrl_' filePrefix '_reactive01.fig'])
fig.Position = [0 0 570 510];

















