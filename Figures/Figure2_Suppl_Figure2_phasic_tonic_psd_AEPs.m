clear all;
close all;

addpath(genpath('\\surrey.ac.uk\personal\hs301\m17462\matlab\eeglab')); % eeglab toolbox, see README on where to find this
addpath(genpath('\\surrey.ac.uk\personal\hs301\m17462\matlab\Henry\useful_functions')); % contains linspecer function, circular statistics toolbox functions, echt function, shadedErrorBar function, see README on where to find this
addpath(genpath('\\surrey.ac.uk\personal\hs301\m17462\matlab\colorGradient'));  % colorGradient function, see README on where to find this

Savefolder = 'D:\Valeria\RSN\data\for_sharing\data_to_make_figures\Figures';

load('D:\Valeria\RSN\data\for_sharing\data_to_make_figures\psd_allsub_mICA_avref_12-Mar-2023.mat');
load('D:\Valeria\RSN\data\for_sharing\data_to_make_figures\power_allsub_mICA_avref_09-Mar-2023.mat');

load('D:\Valeria\RSN\data\for_sharing\data_to_make_figures\ERP_allsub_REM_mICA_avref04-Jun-2023.mat');
ERP_REM = ERP_all;
clear ERP_all

load('D:\Valeria\RSN\data\for_sharing\data_to_make_figures\ERP_allsub_wake_mICA_avref02-Jun-2023.mat');
ERP_wake = ERP_all;
clear ERP_all

incl_sub = setdiff(1:19,12);

conditions = {'Alpha Phase 0'; 'Alpha Phase 90'; 'Alpha Phase 180'; 'Alpha Phase 270';...
            'Theta Phase 0'; 'Theta Phase 90'; 'Theta Phase 180'; 'Theta Phase 270';}


%% Average across on and off blocks and calculate change
ch = [2 34 65 94]; % Fz, AFz, AFF1h, AFF2h (4 closest electrodes to phase-locking)

mpsd_con_ch_stage = mpsd_con_ch; % mpsd_con_ch: sub x ch x con x bins x ep

mpsd_con_cluster = squeeze(nanmean(mpsd_con_ch(:,ch,:,:,:),2));
mpsd_con_cluster_phasic = squeeze(nanmean(mpsd_con_ch_phasic(:,ch,:,:,:),2));
mpsd_con_cluster_tonic = squeeze(nanmean(mpsd_con_ch_tonic(:,ch,:,:,:),2));

on_block = 7:12;
off_block = 1:6;

psd_ON = nanmean(mpsd_con_cluster(:,:,:,on_block),4);
psd_OFF = nanmean(mpsd_con_cluster(:,:,:,off_block),4);

psd_ON_phasic = nanmean(mpsd_con_cluster_phasic(:,:,:,on_block),4);
psd_OFF_phasic = nanmean(mpsd_con_cluster_phasic(:,:,:,off_block),4);

psd_ON_tonic = nanmean(mpsd_con_cluster_tonic(:,:,:,on_block),4);
psd_OFF_tonic = nanmean(mpsd_con_cluster_tonic(:,:,:,off_block),4);

% psd_ONOFF_change = (psd_ON - psd_OFF)./psd_OFF *100;

psd_OFF_phasic_mallcon  = squeeze(nanmean(psd_OFF_phasic,2));
mpsd_OFF_phasic_mallcon = nanmean(psd_OFF_phasic_mallcon(incl_sub,:),1);
sem_psd_OFF_phasic_mallcon = nanstd(psd_OFF_phasic_mallcon(incl_sub,:),1)./sqrt(length(incl_sub));

psd_OFF_tonic_mallcon  = squeeze(nanmean(psd_OFF_tonic,2));
mpsd_OFF_tonic_mallcon = nanmean(psd_OFF_tonic_mallcon(incl_sub,:),1);
sem_psd_OFF_tonic_mallcon = nanstd(psd_OFF_tonic_mallcon(incl_sub,:),1)./sqrt(length(incl_sub));

psd_wake_EO_e_cluster = squeeze(nanmean(mfft_wake_EO_e(:,ch,:),2));
psd_wake_EO_m_cluster = squeeze(nanmean(mfft_wake_EO_m(:,ch,:),2));

psd_wake_EC_e_cluster = squeeze(nanmean(mfft_wake_EC_e(:,ch,:),2));
psd_wake_EC_m_cluster = squeeze(nanmean(mfft_wake_EC_m(:,ch,:),2));

for s = 1:size(psd_wake_EO_e_cluster,1)
    
    psd_wake_EO_cluster(s,:) = nanmean(vertcat(psd_wake_EO_e_cluster(s,:),psd_wake_EO_m_cluster(s,:)));
    psd_wake_EC_cluster(s,:) = nanmean(vertcat(psd_wake_EC_e_cluster(s,:),psd_wake_EC_m_cluster(s,:)));
    
end

mpsd_wake_EO_cluster = nanmean(psd_wake_EO_cluster(incl_sub,:),1);
sem_psd_wake_EO_cluster = nanstd(psd_wake_EO_cluster(incl_sub,:),1)./sqrt(length(incl_sub));

mpsd_wake_EC_cluster = nanmean(psd_wake_EC_cluster(incl_sub,:),1);
sem_psd_wake_EC_cluster = nanstd(psd_wake_EC_cluster(incl_sub,:),1)./sqrt(length(incl_sub));


%%

m_ERP.tonic = nanmean(ERP_REM.REM_ERP{1}(incl_sub,:),1);
sem_ERP.tonic = nanstd(ERP_REM.REM_ERP{1}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_ERP{1}(incl_sub,:),1));

m_ERP.phasic = nanmean(ERP_REM.REM_ERP{2}(incl_sub,:),1);
sem_ERP.phasic = nanstd(ERP_REM.REM_ERP{2}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_ERP{2}(incl_sub,:),1));

m_ERP.wake_e_EO = nanmean(ERP_wake.wake_e_ERP{1}(incl_sub,:),1);
sem_ERP.wake_e_EO = nanstd(ERP_wake.wake_e_ERP{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_ERP{1}(incl_sub,:),1));

m_ERP.wake_e_EC = nanmean(ERP_wake.wake_e_ERP{2}(incl_sub,:),1);
sem_ERP.wake_e_EC = nanstd(ERP_wake.wake_e_ERP{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_ERP{2}(incl_sub,:),1));

m_ERP.wake_m_EO = nanmean(ERP_wake.wake_m_ERP{1}(incl_sub,:),1);
sem_ERP.wake_m_EO = nanstd(ERP_wake.wake_m_ERP{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_ERP{1}(incl_sub,:),1));

m_ERP.wake_m_EC = nanmean(ERP_wake.wake_m_ERP{2}(incl_sub,:),1);
sem_ERP.wake_m_EC = nanstd(ERP_wake.wake_m_ERP{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_ERP{2}(incl_sub,:),1));


m_ERP.tonic_rand = nanmean(ERP_REM.REM_ERP_rand{1}(incl_sub,:),1);
sem_ERP.tonic_rand = nanstd(ERP_REM.REM_ERP_rand{1}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_ERP_rand{1}(incl_sub,:),1));

m_ERP.phasic_rand = nanmean(ERP_REM.REM_ERP_rand{2}(incl_sub,:),1);
sem_ERP.phasic_rand = nanstd(ERP_REM.REM_ERP_rand{2}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_ERP_rand{2}(incl_sub,:),1));

m_ERP.wake_e_EO_rand = nanmean(ERP_wake.wake_e_ERP_rand{1}(incl_sub,:),1);
sem_ERP.wake_e_EO_rand = nanstd(ERP_wake.wake_e_ERP_rand{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_ERP_rand{1}(incl_sub,:),1));

m_ERP.wake_e_EC_rand = nanmean(ERP_wake.wake_e_ERP_rand{2}(incl_sub,:),1);
sem_ERP.wake_e_EC_rand = nanstd(ERP_wake.wake_e_ERP_rand{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_ERP_rand{2}(incl_sub,:),1));

m_ERP.wake_m_EO_rand = nanmean(ERP_wake.wake_m_ERP_rand{1}(incl_sub,:),1);
sem_ERP.wake_m_EO_rand = nanstd(ERP_wake.wake_m_ERP_rand{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_ERP_rand{1}(incl_sub,:),1));

m_ERP.wake_m_EC_rand = nanmean(ERP_wake.wake_m_ERP_rand{2}(incl_sub,:),1);
sem_ERP.wake_m_EC_rand = nanstd(ERP_wake.wake_m_ERP_rand{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_ERP_rand{2}(incl_sub,:),1));


m_ERP.tonic_vol80 = nanmean(ERP_REM.REM_vol80{1}(incl_sub,:),1)
sem_ERP.tonic_vol80 = nanstd(ERP_REM.REM_vol80{1}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_vol80{1}(incl_sub,:),1));

m_ERP.tonic_vol85 = nanmean(ERP_REM.REM_vol85{1}(incl_sub,:),1)
sem_ERP.tonic_vol85 = nanstd(ERP_REM.REM_vol85{1}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_vol85{1}(incl_sub,:),1));

m_ERP.tonic_vol90 = nanmean(ERP_REM.REM_vol90{1}(incl_sub,:),1)
sem_ERP.tonic_vol90 = nanstd(ERP_REM.REM_vol90{1}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_vol90{1}(incl_sub,:),1));


m_ERP.phasic_vol80 = nanmean(ERP_REM.REM_vol80{2}(incl_sub,:),1)
sem_ERP.phasic_vol80 = nanstd(ERP_REM.REM_vol80{2}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_vol80{2}(incl_sub,:),1));

m_ERP.phasic_vol85 = nanmean(ERP_REM.REM_vol85{2}(incl_sub,:),1)
sem_ERP.phasic_vol85 = nanstd(ERP_REM.REM_vol85{2}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_vol85{2}(incl_sub,:),1));

m_ERP.phasic_vol90 = nanmean(ERP_REM.REM_vol90{2}(incl_sub,:),1)
sem_ERP.phasic_vol90 = nanstd(ERP_REM.REM_vol90{2}(incl_sub,:),1)./sqrt(size(ERP_REM.REM_vol90{2}(incl_sub,:),1));


m_ERP.wake_e_EC_vol80 = nanmean(ERP_wake.wake_e_vol80{2}(incl_sub,:),1);
sem_ERP.wake_e_vol80 = nanstd(ERP_wake.wake_e_vol80{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_vol80{2}(incl_sub,:),1));

m_ERP.wake_e_EC_vol85 = nanmean(ERP_wake.wake_e_vol85{2}(incl_sub,:),1);
sem_ERP.wake_e_vol85 = nanstd(ERP_wake.wake_e_vol85{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_vol85{2}(incl_sub,:),1));

m_ERP.wake_e_EC_vol90 = nanmean(ERP_wake.wake_e_vol90{2}(incl_sub,:),1);
sem_ERP.wake_e_vol90 = nanstd(ERP_wake.wake_e_vol90{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_vol90{2}(incl_sub,:),1));


m_ERP.wake_m_EC_vol80 = nanmean(ERP_wake.wake_m_vol80{2}(incl_sub,:),1);
sem_ERP.wake_m_vol80 = nanstd(ERP_wake.wake_m_vol80{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_vol80{2}(incl_sub,:),1));

m_ERP.wake_m_EC_vol85 = nanmean(ERP_wake.wake_m_vol85{2}(incl_sub,:),1);
sem_ERP.wake_m_vol85 = nanstd(ERP_wake.wake_m_vol85{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_vol85{2}(incl_sub,:),1));

m_ERP.wake_m_EC_vol90 = nanmean(ERP_wake.wake_m_vol90{2}(incl_sub,:),1);
sem_ERP.wake_m_vol90 = nanstd(ERP_wake.wake_m_vol90{2}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_vol90{2}(incl_sub,:),1));


m_ERP.wake_e_EO_vol80 = nanmean(ERP_wake.wake_e_vol80{1}(incl_sub,:),1);
sem_ERP.wake_e_vol80 = nanstd(ERP_wake.wake_e_vol80{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_vol80{1}(incl_sub,:),1));

m_ERP.wake_e_EO_vol85 = nanmean(ERP_wake.wake_e_vol85{1}(incl_sub,:),1);
sem_ERP.wake_e_vol85 = nanstd(ERP_wake.wake_e_vol85{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_vol85{1}(incl_sub,:),1));

m_ERP.wake_e_EO_vol90 = nanmean(ERP_wake.wake_e_vol90{1}(incl_sub,:),1);
sem_ERP.wake_e_vol90 = nanstd(ERP_wake.wake_e_vol90{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_e_vol90{1}(incl_sub,:),1));


m_ERP.wake_m_EO_vol80 = nanmean(ERP_wake.wake_m_vol80{1}(incl_sub,:),1);
sem_ERP.wake_m_vol80 = nanstd(ERP_wake.wake_m_vol80{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_vol80{1}(incl_sub,:),1));

m_ERP.wake_m_EO_vol85 = nanmean(ERP_wake.wake_m_vol85{1}(incl_sub,:),1);
sem_ERP.wake_m_vol85 = nanstd(ERP_wake.wake_m_vol85{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_vol85{1}(incl_sub,:),1));

m_ERP.wake_m_EO_vol90 = nanmean(ERP_wake.wake_m_vol90{1}(incl_sub,:),1);
sem_ERP.wake_m_vol90 = nanstd(ERP_wake.wake_m_vol90{1}(incl_sub,:),1)./sqrt(size(ERP_wake.wake_m_vol90{1}(incl_sub,:),1));


%%

for b = 1:size(psd_wake_EO_cluster,2)

  [h_wake(b) p_wake(b)] = ttest(psd_wake_EO_cluster(incl_sub,b),psd_wake_EC_cluster(incl_sub,b));

end

for b = 1:size(psd_wake_EO_cluster,2)

  [h_REM(b) p_REM(b)] = ttest(psd_OFF_phasic_mallcon(incl_sub,b),psd_OFF_tonic_mallcon(incl_sub,b));

end

%% PSD wake eyes open vs eyes closed

colors = linspecer(9);


fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])    

shadedErrorBar(f,mpsd_wake_EO_cluster,sem_psd_wake_EO_cluster,'lineProps',{'-','Color',[0.3 0.3 0.8],'LineWidth',3},'patchSaturation',.1);
hold on
shadedErrorBar(f,mpsd_wake_EC_cluster,sem_psd_wake_EC_cluster,'lineProps',{'-','Color',[0.5 0.7 0.9],'LineWidth',3},'patchSaturation',.1);
% plot(f,mpsd_wake_EO_cluster,'Color',[0.3 0.3 0.8],'LineWidth',3);
% hold on
% plot(f,mpsd_wake_EC_cluster,'Color',[0.5 0.7 0.9],'LineWidth',3);
xlabel('Frequency (Hz)');
ylabel('Power density (\muV^2/s)')
% ylabel('Log Power');
% title(['Mean allsub, ch ', EEG.chanlocs(ch).labels]);
xlim([0 18]);
xticks(0:2:18);
ylim([-0.5 6]);
% legend({'wake EO' 'wake EC'})
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square

sig_bins_wake = find(p_wake <= 0.05);
plot(f(sig_bins_wake),ones(length(sig_bins_wake),1)*-0.2,'*','Color','k');


saveas(fig,[Savefolder,'Figure2A_psd_EO_EC_',num2str(ch),'.svg']);


%% PSD phasic, tonic off

colors = linspecer(9);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])    

shadedErrorBar(f,mpsd_OFF_phasic_mallcon,sem_psd_OFF_phasic_mallcon,'lineProps',{'-','Color',[0.6365 0.3753 0.6753],'LineWidth',3},'patchSaturation',.1);
hold on
shadedErrorBar(f,mpsd_OFF_tonic_mallcon,sem_psd_OFF_tonic_mallcon,'lineProps',{'-','Color',[1 0.3 0.5],'LineWidth',3},'patchSaturation',.1);
% semilogy(f,mpsd_OFF_phasic_mallcon,'Color',[0.6365 0.3753 0.6753],'LineWidth',3);
% hold on
% semilogy(f,mpsd_OFF_tonic_mallcon,'Color',[1 0.3 0.5],'LineWidth',3);
xlabel('Frequency (Hz)');
ylabel('Power density (\muV^2/s)')
% ylabel('Log Power');
% title(['Mean allsub, ch ', EEG.chanlocs(ch).labels]);
xlim([0 18]);
xticks(0:2:18);
ylim([-0.5 6]);
% ylim([-2 8]);
% legend({'phasic' 'tonic'})
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square

sig_bins_REM = find(p_REM <= 0.05);
plot(f(sig_bins_REM),ones(length(sig_bins_REM),1)*-0.2,'*','Color','k');

saveas(fig,[Savefolder,'Figure2B_psd_phasic_tonic_',num2str(ch),'.svg']);


%% Compare EO vs EC

for s = 1:size(ERP_wake.wake_e_ERP{1},1)
    
    ERP_wake_EO(s,:) = nanmean(vertcat(ERP_wake.wake_e_ERP{1}(s,:),ERP_wake.wake_m_ERP{1}(s,:)),1);
    ERP_wake_EC(s,:) = nanmean(vertcat(ERP_wake.wake_e_ERP{2}(s,:),ERP_wake.wake_m_ERP{2}(s,:)),1);

end

mERP_wake_EO = nanmean(ERP_wake_EO,1);
sem_ERP_wake_EO = nanstd(ERP_wake_EO,1)./sqrt(length(incl_sub));

mERP_wake_EC = nanmean(ERP_wake_EC,1);
sem_ERP_wake_EC = nanstd(ERP_wake_EC,1)./sqrt(length(incl_sub));

for samp = 1:size(ERP_wake.wake_e_ERP_rand{1},2)

[h_EOEC(samp) p_EOEC(samp)] = ttest(ERP_wake_EO(incl_sub,samp), ERP_wake_EC(incl_sub,samp));
     
end


%% ERP wake eve EO vs EC

t = (-dt*fs:dt*fs-1)/fs*1000;

colors = linspecer(8);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar(t,mERP_wake_EO,sem_ERP_wake_EO,'lineProps',{'Color',[0.3 0.3 0.8],'LineWidth',3},'patchSaturation',.3);
% plot(t,m_ERP.wake_e_EO,'Color',[0.3 0.3 0.8],'LineWidth',3);
hold on
shadedErrorBar(t,mERP_wake_EC,sem_ERP_wake_EC,'lineProps',{'Color',[0.5 0.7 0.9],'LineWidth',3},'patchSaturation',.3);
% plot(t,m_ERP.wake_e_EC,'Color',[0.3 0.5 0.8],'LineWidth',3);
hold on
xline(0,'LineStyle','--','LineWidth',5);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'phasic' 'tonic' });
legend off
ylim([-3 4]);
xticks(-200:200:1000);

sig_bins_ttest = find(p_EOEC <= 0.05);
plot(t(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-2.5,'*','Color','k');

saveas(fig,[Savefolder,'Figure2C_ERP_EOEC_notmatched.svg']);


%% Compare phasic vs tonic (not matched)

for samp = 1:size(ERP_REM.REM_ERP_rand{1},2)

[h_phasictonic(samp) p_phasictonic(samp)] = ttest(ERP_REM.REM_ERP{1}(incl_sub,samp), ERP_REM.REM_ERP{2}(incl_sub,samp));
     
end

%% ERP for phasic and tonic (not matched)

t = (-dt*fs:dt*fs-1)/fs*1000;

colors = linspecer(8);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar(t,m_ERP.phasic,sem_ERP.phasic,'lineProps',{'Color',[0.6365 0.3753 0.6753],'LineWidth',1},'patchSaturation',.3);
% plot(t,m_ERP.phasic,'Color',[1 0.3 0.2],'LineWidth',3);
hold on
shadedErrorBar(t,m_ERP.tonic,sem_ERP.tonic,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',1},'patchSaturation',.3);
% plot(t,m_ERP.tonic,'Color',[1 0.3 0.5],'LineWidth',3);
hold on
xline(0,'LineStyle','--','LineWidth',5);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'phasic' 'tonic' });
legend off
ylim([-3 4]);
xticks(-200:200:1000);

% sig_bins_lme = find(p_state <= 0.05);
% plot(t(sig_bins_lme),ones(length(sig_bins_lme),1)*3.5,'*','Color','k');

sig_bins_ttest = find(p_phasictonic <= 0.05);
plot(t(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-2.5,'*','Color','k');


sig_times = t(sig_bins_ttest)

saveas(fig,[Savefolder,'Figure2D_ERP_phasictonic_notmatched.svg']);


%% Compare phasic vs tonic (matched)

for samp = 1:size(ERP_REM.REM_ERP_rand{1},2)

[h_phasictonic_matched(samp) p_phasictonic_matched(samp)] = ttest(ERP_REM.REM_ERP_rand{1}(incl_sub,samp), ERP_REM.REM_ERP_rand{2}(incl_sub,samp));
     
end

%% ERP for phasic and tonic (matched)

t = (-dt*fs:dt*fs-1)/fs*1000;

colors = linspecer(8);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar(t,m_ERP.phasic_rand,sem_ERP.phasic_rand,'lineProps',{'Color',[0.6365 0.3753 0.6753],'LineWidth',1},'patchSaturation',.3);
% plot(t,m_ERP.phasic,'Color',[1 0.3 0.2],'LineWidth',3);
hold on
shadedErrorBar(t,m_ERP.tonic_rand,sem_ERP.tonic_rand,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',1},'patchSaturation',.3);
% plot(t,m_ERP.tonic,'Color',[1 0.3 0.5],'LineWidth',3);
hold on
xline(0,'LineStyle','--','LineWidth',5);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'phasic' 'tonic' });
legend off
ylim([-3 4]);
xticks(-200:200:1000);

% sig_bins_lme = find(p_state <= 0.05);
% plot(t(sig_bins_lme),ones(length(sig_bins_lme),1)*3.5,'*','Color','k');

sig_bins_ttest = find(p_phasictonic_matched <= 0.05);
plot(t(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-2.5,'*','Color','k');

sig_times = t(sig_bins_ttest)
diff_ndx = find(diff(sig_bins_ttest)>1);

for d = 2:length(diff_ndx)
start_diff(d) = t(sig_bins_ttest(diff_ndx(d-1)+1));
end_diff(d) = t(sig_bins_ttest(diff_ndx(d)));

end

saveas(fig,[Savefolder,'Suppl_Figure2A_ERP_phasictonic_matched.svg']);


%% Compare eve vs mor

for s = 1:size(ERP_wake.wake_e_ERP{1},1)
    
    ERP_wake_eve(s,:) = nanmean(vertcat(ERP_wake.wake_e_ERP{1}(s,:),ERP_wake.wake_e_ERP{2}(s,:)),1);
    ERP_wake_mor(s,:) = nanmean(vertcat(ERP_wake.wake_m_ERP{1}(s,:),ERP_wake.wake_m_ERP{2}(s,:)),1);

end

mERP_wake_eve = nanmean(ERP_wake_eve,1);
sem_ERP_wake_eve = nanstd(ERP_wake_eve,1)./sqrt(length(incl_sub));

mERP_wake_mor = nanmean(ERP_wake_mor,1);
sem_ERP_wake_mor = nanstd(ERP_wake_mor,1)./sqrt(length(incl_sub));


for samp = 1:size(ERP_wake.wake_e_ERP_rand{1},2)

[h_evemor(samp) p_evemor(samp)] = ttest(ERP_wake_eve(incl_sub,samp), ERP_wake_mor(incl_sub,samp));
     
end


%% ERP eve vs mor

t = (-dt*fs:dt*fs-1)/fs*1000;

colors = linspecer(8);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar(t,mERP_wake_eve,sem_ERP_wake_eve,'lineProps',{'Color',[0.3 0.5 0.3],'LineWidth',3},'patchSaturation',.3);
% plot(t,m_ERP.wake_e_EO,'Color',[0.3 0.3 0.8],'LineWidth',3);
hold on
shadedErrorBar(t,mERP_wake_mor,sem_ERP_wake_mor,'lineProps',{'Color',[0.5 0.9 0.3],'LineWidth',3},'patchSaturation',.3);
% plot(t,m_ERP.wake_e_EC,'Color',[0.3 0.5 0.8],'LineWidth',3);
hold on
xline(0,'LineStyle','--','LineWidth',5);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'evening' 'morning' });
% legend off
ylim([-3 4]);
xticks(-200:200:1000);

sig_bins_ttest = find(p_EOEC <= 0.05);
plot(t(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-2.5,'*','Color','k');

saveas(fig,[Savefolder,'Suppl_Figure2B_ERP_evemor_notmatched.svg']);


%% Compare tonic volumes using a lme (not matched)

incl_sub = setdiff(1:19,12);

for samp = 1:size(ERP_REM.REM_vol80{1},2)
   
    table_vol_samp = [];
    
    ERP_tonic_vol80_samp = ERP_REM.REM_vol80{1}(incl_sub,samp);
    ERP_tonic_vol85_samp = ERP_REM.REM_vol85{1}(incl_sub,samp);
    ERP_tonic_vol90_samp = ERP_REM.REM_vol90{1}(incl_sub,samp);
    ERP_samp = vertcat(ERP_tonic_vol80_samp,ERP_tonic_vol85_samp,ERP_tonic_vol90_samp );
    
    vol = vertcat(repmat(80,length(incl_sub),1),repmat(85,length(incl_sub),1),repmat(90,length(incl_sub),1));
    sub_table = vertcat(incl_sub',incl_sub',incl_sub');
    table_vol_samp = table(sub_table,vol,ERP_samp,'VariableNames',{'sub','vol','ERP'});
    table_vol_samp.sub = categorical(table_vol_samp.sub);
    table_vol_samp.vol = categorical(table_vol_samp.vol);
    lme_vol = fitlme(table_vol_samp,'ERP ~ vol +(1|sub)','FitMethod','REML','DummyVarCoding','effects');
    stats = anova(lme_vol);
    p_vol_tonic(samp) = stats.pValue(2);
    F_vol_tonic(samp) = stats.FStat(2);
        
end

%% ERP for different volumes for tonic (not matched)

t = (-dt*fs:dt*fs-1)/fs*1000;

colors = colorGradient([0.9290 0.6940 0.1250],[0.6350 0.0780 0.1840],3);
% colors = linspecer(3);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
shadedErrorBar(t,m_ERP.tonic_vol80,sem_ERP.tonic_vol80,'lineProps',{'Color',colors(1,:),'LineWidth',3},'patchSaturation',.3);
hold on
shadedErrorBar(t,m_ERP.tonic_vol85,sem_ERP.tonic_vol85,'lineProps',{'Color',colors(2,:),'LineWidth',3},'patchSaturation',.3);
hold on
shadedErrorBar(t,m_ERP.tonic_vol90,sem_ERP.tonic_vol90,'lineProps',{'Color',colors(3,:),'LineWidth',3},'patchSaturation',.3);
hold on
xline(0,'LineStyle','--','LineWidth',5);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'vol50' 'vol55' 'vol60' });
% legend off
ylim([-4 4]);
xticks(-200:200:1000);

sig_bins_lme = find(p_vol_tonic <= 0.05);
plot(t(sig_bins_lme),ones(length(sig_bins_lme),1)*-3.75,'*','Color','k');

saveas(fig,[Savefolder,'Suppl_Figure2C_ERP_tonic_vol.svg']);

%% Compare phasic volumes using a lme (not matched)

incl_sub = setdiff(1:19,12);

for samp = 1:size(ERP_REM.REM_vol80{2},2)
   
    table_vol_samp = [];
    
    ERP_phasic_vol80_samp = ERP_REM.REM_vol80{2}(incl_sub,samp);
    ERP_phasic_vol85_samp = ERP_REM.REM_vol85{2}(incl_sub,samp);
    ERP_phasic_vol90_samp = ERP_REM.REM_vol90{2}(incl_sub,samp);
    ERP_samp = vertcat(ERP_phasic_vol80_samp,ERP_phasic_vol85_samp,ERP_phasic_vol90_samp );
    
    vol = vertcat(repmat(80,length(incl_sub),1),repmat(85,length(incl_sub),1),repmat(90,length(incl_sub),1));
    sub_table = vertcat(incl_sub',incl_sub',incl_sub');
    table_vol_samp = table(sub_table,vol,ERP_samp,'VariableNames',{'sub','vol','ERP'});
    table_vol_samp.sub = categorical(table_vol_samp.sub);
    table_vol_samp.vol = categorical(table_vol_samp.vol);
    lme_vol = fitlme(table_vol_samp,'ERP ~ vol +(1|sub)','FitMethod','REML','DummyVarCoding','effects');
    stats = anova(lme_vol);
    p_vol_phasic(samp) = stats.pValue(2);
    F_vol_phasic(samp) = stats.FStat(2);
        
end

%% ERP for different volumes for phasic (not matched)


t = (-dt*fs:dt*fs-1)/fs*1000;

colors = colorGradient([0.9290 0.6940 0.1250],[0.6350 0.0780 0.1840],3);
% colors = linspecer(3);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
shadedErrorBar(t,m_ERP.phasic_vol80,sem_ERP.phasic_vol80,'lineProps',{'Color',colors(1,:),'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(t,m_ERP.phasic_vol85,sem_ERP.phasic_vol85,'lineProps',{'Color',colors(2,:),'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(t,m_ERP.phasic_vol90,sem_ERP.phasic_vol90,'lineProps',{'Color',colors(3,:),'LineWidth',1},'patchSaturation',.3);
hold on
xline(0,'LineStyle','--','LineWidth',5);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'vol50' 'vol55' 'vol60' });
% legend off
ylim([-4 4]);
xticks(-200:200:1000);

sig_bins_lme = find(p_vol_phasic <= 0.05);
plot(t(sig_bins_lme),ones(length(sig_bins_lme),1)*-3.75,'*','Color','k');

saveas(fig,[Savefolder,'Suppl_Figure2C_ERP_phasic_vol.svg']);
