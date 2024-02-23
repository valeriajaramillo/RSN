% clear all;
% close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

addpath(genpath('/user/HS301/m17462/matlab/ScientificColourMaps7'));
addpath(genpath('/user/HS301/m17462/matlab/kispi'));

load('/vol/research/nemo/datasets/RSN/data/analysis/power_allsub/power_allsub_mICA_avref_09-Mar-2023.mat');

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/power_allsub/';

load('/user/HS301/m17462/matlab/Scripts/RSN/analyses/psd/topo/EEG_chanlocs.mat');

incl_sub = setdiff(1:19,12);

%% Average across on and off blocks and calculate change
ch = [2 34 65 94]; % Fz, AFz, AFF1h, AFF2h
% ch = 16:18; % O1, Oz, Oz
lower_freq = 7;
higher_freq = 12;

f = 0:0.1:40;

%% wake eve EC

colors = linspecer(10);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for s = 1:length(incl_sub)
    
    mfft_wake_e_cluster = squeeze(nanmean(mfft_wake_EC_e(incl_sub(s),ch,:),2));
%     npower = mfft_wake_e_cluster/sum(mfft_wake_e_cluster); %mfft_wake_e_cluster; %
    power_wake_e_EC(s,:) = mfft_wake_e_cluster;
    zpower = (mfft_wake_e_cluster-nanmean(mfft_wake_e_cluster))/nanstd(mfft_wake_e_cluster); %mfft_wake_e_cluster; %
    zpower_wake_e_EC(s,:) = zpower;

    plot(f,zpower,'Color',colors(1,:),'LineWidth',2);
    hold on

    [pks,locs] = findpeaks(power_wake_e_EC(s,:));
%     [pks,locs] = findpeaks(zpower);
    
    indiv_peaks(s,locs) = 1;
            
    IAPF_wake_e_EC(s) = NaN;
    
    if find(f(locs)>=7 & f(locs)<=higher_freq)
        
        [row,col] = find(f(locs)>=lower_freq & f(locs)<=higher_freq);
        
        if length(col)>1
           
            [r,c] = find(pks(col)==max(pks(col)));
            
            IAPF_wake_e_EC(s) = f(locs(col(r)));
            scatter(IAPF_wake_e_EC(s),pks(col(r)),50,colors(1,:),'filled')
        else
            
            IAPF_wake_e_EC(s) = f(locs(col));
            scatter(IAPF_wake_e_EC(s),pks(col),50,colors(1,:),'filled')
        end
        
        
    end


end

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Normalized Power (z-score)');
title(['Wake eve EC, ch ', EEG.chanlocs(ch).labels]);
xlim([3 15]);
% saveas(fig,[Savefolder,'indsub_wake_eve_EC_ch',EEG.chanlocs(ch).labels,'.svg']);

%% wake mor EC

colors = linspecer(12);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for s = 1:length(incl_sub)
    
    mfft_wake_m_cluster = squeeze(nanmean(mfft_wake_EC_m(incl_sub(s),ch,:),2));
%     power = mfft_wake_m_cluster; %/sum(mfft_wake_e_cluster); %mfft_wake_e_cluster; %
    power_wake_m_EC(s,:) = mfft_wake_m_cluster;
    zpower = (mfft_wake_m_cluster-nanmean(mfft_wake_m_cluster))/nanstd(mfft_wake_m_cluster); %mfft_wake_e_cluster; %
    zpower_wake_m_EC(s,:) = zpower;

    plot(f,zpower,'Color',colors(2,:),'LineWidth',2);
    hold on

    [pks,locs] = findpeaks(power_wake_m_EC(s,:));
%     [pks,locs] = findpeaks(zpower);
    
    indiv_peaks(s,locs) = 1;
            
    IAPF_wake_m_EC(s) = NaN;
    
    if find(f(locs)>=lower_freq & f(locs)<=higher_freq)
        
        [row,col] = find(f(locs)>=lower_freq & f(locs)<=higher_freq);
        
        if length(col)>1
           
            [r,c] = find(pks(col)==max(pks(col)));
            
            IAPF_wake_m_EC(s) = f(locs(col(r)));
            scatter(IAPF_wake_m_EC(s),pks(col(r)),50,colors(2,:),'filled')
        else
            
            IAPF_wake_m_EC(s) = f(locs(col));
            scatter(IAPF_wake_m_EC(s),pks(col),50,colors(2,:),'filled')
        end
        
        
    end


end

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Normalized Power (z-score)');
title(['Wake mor EC, ch ', EEG.chanlocs(ch).labels]);
xlim([3 15]);
% saveas(fig,[Savefolder,'indsub_wake_mor_EC_ch',EEG.chanlocs(ch).labels,'.svg']);


%% wake eve EO

colors = linspecer(10);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for s = 1:length(incl_sub)
    
    mfft_wake_e_cluster = squeeze(nanmean(mfft_wake_EO_e(incl_sub(s),ch,:),2));
%     npower = log(mfft_wake_e_cluster); %/sum(mfft_wake_e_cluster); %mfft_wake_e_cluster; %
    power_wake_e_EO(s,:) = mfft_wake_e_cluster;
    zpower = (mfft_wake_e_cluster-nanmean(mfft_wake_e_cluster))/nanstd(mfft_wake_e_cluster); %mfft_wake_e_cluster; %
    zpower_wake_e_EO(s,:) = zpower;

    plot(f,zpower,'Color',colors(1,:),'LineWidth',2);
    hold on

    [pks,locs] = findpeaks(power_wake_e_EO(s,:));
%     [pks,locs] = findpeaks(zpower);
    
    indiv_peaks(s,locs) = 1;
            
    IAPF_wake_e_EO(s) = NaN;
    
    if find(f(locs)>=7 & f(locs)<=higher_freq)
        
        [row,col] = find(f(locs)>=lower_freq & f(locs)<=higher_freq);
        
        if length(col)>1
           
            [r,c] = find(pks(col)==max(pks(col)));
            
            IAPF_wake_e_EO(s) = f(locs(col(r)));
            scatter(IAPF_wake_e_EO(s),pks(col(r)),50,colors(1,:),'filled')
        else
            
            IAPF_wake_e_EO(s) = f(locs(col));
            scatter(IAPF_wake_e_EO(s),pks(col),50,colors(1,:),'filled')
        end
        
        
    end


end

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Normalized Power (z-score)');
title(['Wake eve EO, ch ', EEG.chanlocs(ch).labels]);
xlim([3 15]);
% saveas(fig,[Savefolder,'indsub_wake_eve_EO_ch',EEG.chanlocs(ch).labels,'.svg']);

%% wake mor EO

colors = linspecer(12);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for s = 1:length(incl_sub)
    
    mfft_wake_m_cluster = squeeze(nanmean(mfft_wake_EO_m(incl_sub(s),ch,:),2));
%     npower = log(mfft_wake_m_cluster); %/sum(mfft_wake_m_cluster); %mfft_wake_e_cluster; %
    power_wake_m_EO(s,:) = mfft_wake_m_cluster;
    zpower = (mfft_wake_m_cluster-nanmean(mfft_wake_m_cluster))/nanstd(mfft_wake_m_cluster); %mfft_wake_e_cluster; %
    zpower_wake_m_EO(s,:) = zpower;

    plot(f,zpower,'Color',colors(2,:),'LineWidth',2);
    hold on

    [pks,locs] = findpeaks(power_wake_m_EO(s,:));
%     [pks,locs] = findpeaks(zpower);
    
    indiv_peaks(s,locs) = 1;
            
    IAPF_wake_m_EO(s) = NaN;
    
    if find(f(locs)>=lower_freq & f(locs)<=higher_freq)
        
        [row,col] = find(f(locs)>=lower_freq & f(locs)<=higher_freq);
        
        if length(col)>1
           
            [r,c] = find(pks(col)==max(pks(col)));
            
            IAPF_wake_m_EO(s) = f(locs(col(r)));
            scatter(IAPF_wake_m_EO(s),pks(col(r)),50,colors(2,:),'filled')
        else
            
            IAPF_wake_m_EO(s) = f(locs(col));
            scatter(IAPF_wake_m_EO(s),pks(col),50,colors(2,:),'filled')
        end
        
        
    end


end

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Normalized Power (z-score)');
title(['Wake mor EO, ch ', EEG.chanlocs(ch).labels]);
xlim([3 15]);
% saveas(fig,[Savefolder,'indsub_wake_mor_EO_ch',EEG.chanlocs(ch).labels,'.svg']);


%% phasic REM

colors = linspecer(12);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for s = 1:length(incl_sub)
    
    mfft_phasic_cluster = squeeze(nanmean(mfft_phasic(incl_sub(s),ch,:),2));
%     npower = log(mfft_phasic_cluster); %/sum(mfft_phasic_cluster); %mfft_wake_e_cluster; %
    power_phasic(s,:) = mfft_phasic_cluster;
    zpower = (mfft_phasic_cluster-nanmean(mfft_phasic_cluster))/nanstd(mfft_phasic_cluster); %mfft_wake_e_cluster; %
    zpower_phasic(s,:) = zpower;

    plot(f,zpower,'Color',colors(9,:),'LineWidth',2);
    hold on

    [pks,locs] = findpeaks(zpower_phasic(s,:));
%     [pks,locs] = findpeaks(zpower);
    
    indiv_peaks(s,locs) = 1;
            
    IAPF_phasic(s) = NaN;
    
    if find(f(locs)>=lower_freq & f(locs)<=higher_freq)
        
        [row,col] = find(f(locs)>=lower_freq & f(locs)<=higher_freq);
        
        if length(col)>1
           
            [r,c] = find(pks(col)==max(pks(col)));
            
            IAPF_phasic(s) = f(locs(col(r)));
            scatter(IAPF_phasic(s),pks(col(r)),50,colors(9,:),'filled')
        else
            
            IAPF_phasic(s) = f(locs(col));
            scatter(IAPF_phasic(s),pks(col),50,colors(9,:),'filled')
        end
        
        
    end


end

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Normalized Power (z-score)');
title(['Phasic REM, ch ', EEG.chanlocs(ch).labels]);
xlim([3 15]);
% saveas(fig,[Savefolder,'indsub_phasic_ch',EEG.chanlocs(ch).labels,'.svg']);

%% tonic REM

colors = linspecer(12);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for s = 1:length(incl_sub)
    
    mfft_tonic_cluster = squeeze(nanmean(mfft_tonic(incl_sub(s),ch,:),2));
%     npower = log(mfft_tonic_cluster); %/sum(mfft_tonic_cluster); %mfft_wake_e_cluster; %
    power_tonic(s,:) = mfft_tonic_cluster;
    zpower = (mfft_tonic_cluster-nanmean(mfft_tonic_cluster))/nanstd(mfft_tonic_cluster); %mfft_wake_e_cluster; %
    zpower_tonic(s,:) = zpower;

    plot(f,zpower,'Color',colors(10,:),'LineWidth',2);
    hold on

    [pks,locs] = findpeaks(zpower_tonic(s,:));
%     [pks,locs] = findpeaks(zpower); 

    indiv_peaks(s,locs) = 1;
            
    IAPF_tonic(s) = NaN;
    
    if find(f(locs)>=lower_freq & f(locs)<=higher_freq)
        
        [row,col] = find(f(locs)>=lower_freq & f(locs)<=higher_freq);
        
        if length(col)>1
           
            [r,c] = find(pks(col)==max(pks(col)));
            
            IAPF_tonic(s) = f(locs(col(r)));
            scatter(IAPF_tonic(s),pks(col(r)),50,colors(10,:),'filled')
        else
            
            IAPF_tonic(s) = f(locs(col));
            scatter(IAPF_tonic(s),pks(col),50,colors(10,:),'filled')
        end
        
        
    end


end

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Normalized Power (z-score)');
title(['Tonic REM, ch ', EEG.chanlocs(ch).labels]);
xlim([3 15]);
% saveas(fig,[Savefolder,'indsub_tonic_ch',EEG.chanlocs(ch).labels,'.svg']);

%% plot means - absolute

sem_power_wake_e_EC = nanstd(power_wake_e_EC,1)./sqrt(size(power_wake_e_EC,1));
sem_power_wake_m_EC = nanstd(power_wake_m_EC,1)./sqrt(size(power_wake_m_EC,1));
sem_power_wake_e_EO = nanstd(power_wake_e_EO,1)./sqrt(size(power_wake_e_EO,1));
sem_power_wake_m_EO = nanstd(power_wake_m_EO,1)./sqrt(size(power_wake_m_EO,1));
sem_power_phasic = nanstd(power_phasic,1)./sqrt(size(power_phasic,1));
sem_power_tonic = nanstd(power_tonic,1)./sqrt(size(power_tonic,1));

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
shadedErrorBar(f,nanmean(power_wake_e_EC,1),sem_power_wake_e_EC,'lineProps',{'Color',[0.3 0.5 0.3],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(power_wake_e_EO,1),sem_power_wake_e_EO,'lineProps',{'Color',[0.3 0.3 0.8],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(power_wake_m_EC,1),sem_power_wake_m_EC,'lineProps',{'Color',[0.5 0.9 0.3],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(power_wake_m_EO,1),sem_power_wake_m_EO,'lineProps',{'Color',[0.5 0.7 0.9],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(power_phasic,1),sem_power_phasic,'lineProps',{'Color',[1 0.3 0.2],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(power_tonic,1),sem_power_tonic,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',1},'patchSaturation',.3);

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Absolute Power');
% ylabel('Log Power');
title(['Mean allsub, ch ', EEG.chanlocs(ch).labels]);
xlim([0 15]);
% ylim([-2 8]);
% saveas(fig,[Savefolder,'power_msub_stages_ch',EEG.chanlocs(ch).labels,'.svg']);
legend({'wake eve EC' 'wake eve EO' 'wake mor EC' 'wake mor EO' 'phasic' 'tonic'})

%% plot means - log absolute

sem_logpower_wake_e_EC = nanstd(log(power_wake_e_EC),1)./sqrt(size(log(power_wake_e_EC),1));
sem_logpower_wake_m_EC = nanstd(log(power_wake_m_EC),1)./sqrt(size(log(power_wake_m_EC),1));
sem_logpower_wake_e_EO = nanstd(log(power_wake_e_EO),1)./sqrt(size(log(power_wake_e_EO),1));
sem_logpower_wake_m_EO = nanstd(log(power_wake_m_EO),1)./sqrt(size(log(power_wake_m_EO),1));
sem_logpower_phasic = nanstd(log(power_phasic),1)./sqrt(size(log(power_phasic),1));
sem_logpower_tonic = nanstd(log(power_tonic),1)./sqrt(size(log(power_tonic),1));

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
shadedErrorBar(f,nanmean(log(power_wake_e_EC),1),sem_logpower_wake_e_EC,'lineProps',{'Color',[0.3 0.5 0.3],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(log(power_wake_e_EO),1),sem_logpower_wake_e_EO,'lineProps',{'Color',[0.3 0.3 0.8],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(log(power_wake_m_EC),1),sem_logpower_wake_m_EC,'lineProps',{'Color',[0.5 0.9 0.3],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(log(power_wake_m_EO),1),sem_logpower_wake_m_EO,'lineProps',{'Color',[0.5 0.7 0.9],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(log(power_phasic),1),sem_logpower_phasic,'lineProps',{'Color',[1 0.3 0.2],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(log(power_tonic),1),sem_logpower_tonic,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',1},'patchSaturation',.3);

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Absolute Log Power');
% ylabel('Log Power');
title(['Mean allsub, ch ', EEG.chanlocs(ch).labels]);
xlim([0 15]);
% ylim([-2 8]);
saveas(fig,[Savefolder,'logpower_msub_stages_ch',EEG.chanlocs(ch).labels,'.svg']);
legend({'wake eve EC' 'wake eve EO' 'wake mor EC' 'wake mor EO' 'phasic' 'tonic'})

%% Plot means - z-scored

sem_zpower_wake_e_EC = nanstd(zpower_wake_e_EC,1)./sqrt(size(zpower_wake_e_EC,1));
sem_zpower_wake_m_EC = nanstd(zpower_wake_m_EC,1)./sqrt(size(zpower_wake_m_EC,1));
sem_zpower_wake_e_EO = nanstd(zpower_wake_e_EO,1)./sqrt(size(zpower_wake_e_EO,1));
sem_zpower_wake_m_EO = nanstd(zpower_wake_m_EO,1)./sqrt(size(zpower_wake_m_EO,1));
sem_zpower_phasic = nanstd(zpower_phasic,1)./sqrt(size(zpower_phasic,1));
sem_zpower_tonic = nanstd(zpower_tonic,1)./sqrt(size(zpower_tonic,1));

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
shadedErrorBar(f,nanmean(zpower_wake_e_EC,1),sem_zpower_wake_e_EC,'lineProps',{'Color',[0.3 0.5 0.3],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(zpower_wake_e_EO,1),sem_zpower_wake_e_EO,'lineProps',{'Color',[0.3 0.3 0.8],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(zpower_wake_m_EC,1),sem_zpower_wake_m_EC,'lineProps',{'Color',[0.5 0.9 0.3],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(zpower_wake_m_EO,1),sem_zpower_wake_m_EO,'lineProps',{'Color',[0.5 0.7 0.9],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(zpower_phasic,1),sem_zpower_phasic,'lineProps',{'Color',[1 0.3 0.2],'LineWidth',1},'patchSaturation',.3);
hold on
shadedErrorBar(f,nanmean(zpower_tonic,1),sem_zpower_tonic,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',1},'patchSaturation',.3);

set(gca,'Fontsize',35);
box on
axis square
xlabel('Frequency (Hz)');
ylabel('Z-score Power');
% ylabel('Log Power');
title(['Mean allsub, ch ', EEG.chanlocs(ch).labels]);
xlim([0 15]);
% ylim([-2 6]);
% saveas(fig,[Savefolder,'zpower_msub_stages_ch',EEG.chanlocs(ch).labels,'.svg']);
legend({'wake eve EC' 'wake eve EO' 'wake mor EC' 'wake mor EO' 'phasic' 'tonic'})

%%

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
violins = violinPlot([IAPF_wake_e_EC' IAPF_wake_e_EO' IAPF_wake_m_EC' IAPF_wake_m_EO']);
% violins = violinPlot([freq_change_alpha(:,1,ch) freq_change_alpha(:,2,ch) freq_change_alpha(:,3,ch) freq_change_alpha(:,4,ch)]);
ylabel('IAPF (Hz)');
set(gca,'Fontsize',45);
xticklabels({'eEC' 'eEO' 'mEC' 'mEO'});
box on
axis square
title(['Mean allsub, ch ', EEG.chanlocs(ch).labels]);

violins(1).ViolinColor = [0.3 0.5 0.3];
violins(2).ViolinColor = [0.3 0.3 0.8]; 
violins(3).ViolinColor = [0.5 0.9 0.3]; 
violins(4).ViolinColor = [0.5 0.7 0.9];

saveas(fig,[Savefolder,'violinplot_IAPF_evemor',EEG.chanlocs(ch).labels,'.svg']);

%%


IAPF_all = vertcat(IAPF_wake_e_EC',IAPF_wake_e_EO',IAPF_wake_m_EC',IAPF_wake_m_EO');
state = vertcat(repmat(1,length(incl_sub),1),repmat(2,length(incl_sub),1),repmat(3,length(incl_sub),1),repmat(4,length(incl_sub),1));

sub_table = vertcat(incl_sub',incl_sub',incl_sub',incl_sub');
table_IAPF_all = table(sub_table,state,IAPF_all,'VariableNames',{'sub','state','IAPF'});
table_IAPF_all.sub = categorical(table_IAPF_all.sub);
table_IAPF_all.state = categorical(table_IAPF_all.state);
    
lme_IAPF_all = fitlme(table_IAPF_all,'IAPF ~ state +(1|sub)','FitMethod','REML','DummyVarCoding','effects');
stats = anova(lme_IAPF_all);
p_IAPF_all = stats.pValue(2)
F_IAPF_all = stats.FStat(2)

[h p] = ttest(IAPF_wake_e_EC,IAPF_wake_m_EC)
[h p] = ttest(IAPF_wake_e_EO,IAPF_wake_m_EO)

[h p] = ttest(IAPF_wake_e_EC,IAPF_wake_e_EO)
[h p] = ttest(IAPF_wake_m_EC,IAPF_wake_m_EO)

%% Scatter wake eve - wake mor

nonan_sub = intersect(find(~isnan(IAPF_wake_e) == 1),find(~isnan(IAPF_wake_m) == 1));
[r p] = corr(IAPF_wake_e(nonan_sub)',IAPF_wake_m(nonan_sub)');

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(IAPF_wake_e,IAPF_wake_m,100,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('IAPF wake eve EC (Hz)');
ylabel('IAPF wake mor EC (Hz)');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
xlim([8 12]);
ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);


% nonan_sub = intersect(find(~isnan(IAPF_wake_e_EO) == 1),find(~isnan(IAPF_wake_m_EO) == 1));
% [r p] = corr(IAPF_wake_e_EO(nonan_sub)',IAPF_wake_m_EO(nonan_sub)');
% 
% fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% scatter(IAPF_wake_e_EO,IAPF_wake_m_EO,100,'k','LineWidth',3);
% l = lsline;
% l.LineWidth = 3;
% l.Color = 'k';
% xlabel('IAPF wake eve EO (Hz)');
% ylabel('IAPF wake mor EO (Hz)');
% title(['r = ', num2str(r), ', p = ',num2str(p)]);
% set(gca,'Fontsize',35);
% box on
% axis square
% xlim([8 12]);
% ylim([7 12]);
% % saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);

%%

nonan_sub = intersect(find(~isnan(IAPF_wake_e) == 1),find(~isnan(IAPF_phasic) == 1));
[r p] = corr(IAPF_wake_e(nonan_sub)',IAPF_phasic(nonan_sub)')

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(IAPF_wake_e,IAPF_phasic,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('IAPF wake eve EC (Hz)');
ylabel('IAPF phasic REM (Hz)');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
xlim([8 12]);
ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_phasic_ch',EEG.chanlocs(ch).labels,'.svg']);

%%

nonan_sub = intersect(find(~isnan(IAPF_wake_e) == 1),find(~isnan(IAPF_tonic) == 1));
[r p] = corr(IAPF_wake_e(nonan_sub)',IAPF_tonic(nonan_sub)')

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(IAPF_wake_e,IAPF_tonic,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('IAPF wake eve EC (Hz)');
ylabel('IAPF tonic REM (Hz)');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
xlim([8 12]);
ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_tonic_ch',EEG.chanlocs(ch).labels,'.svg']);

%% power correlations

alpha_ndx = find(f >= 7 & f <= 11);
zpower_wake_e_alpha = nanmean(zpower_wake_e(:,alpha_ndx),2);
zpower_wake_m_alpha = nanmean(zpower_wake_m(:,alpha_ndx),2);
zpower_tonic_alpha = nanmean(zpower_tonic(:,alpha_ndx),2);
zpower_phasic_alpha = nanmean(zpower_phasic(:,alpha_ndx),2);

[r p] = corr(zpower_wake_e_alpha,zpower_wake_m_alpha);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(zpower_wake_e_alpha, zpower_wake_m_alpha,100,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('Alpha power wake eve EC');
ylabel('Alpha power wake mor EC');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
% xlim([8 12]);
% ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);

[r p] = corr(zpower_tonic_alpha,zpower_wake_m_alpha);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(zpower_tonic_alpha, zpower_wake_m_alpha,100,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('Alpha power tonic');
ylabel('Alpha power wake mor EC');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
% xlim([8 12]);
% ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);


[r p] = corr(zpower_phasic_alpha,zpower_wake_m_alpha);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(zpower_phasic_alpha, zpower_wake_m_alpha,100,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('Alpha power phasic');
ylabel('Alpha power wake mor EC');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
% xlim([8 12]);
% ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);

zpower_alpha_change = zpower_wake_m_alpha-zpower_wake_e_alpha;

[r p] = corr(zpower_tonic_alpha,zpower_alpha_change);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(zpower_tonic_alpha, zpower_alpha_change,100,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('Alpha power tonic');
ylabel('Alpha power change');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
% xlim([8 12]);
% ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);

[r p] = corr(zpower_phasic_alpha,zpower_alpha_change);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(zpower_phasic_alpha, zpower_alpha_change,100,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('Alpha power phasic');
ylabel('Alpha power change');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
% xlim([8 12]);
% ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);


delta_ndx = find(f >= 1 & f <= 4);
zpower_nrem = nanmean(nanmean(mfft_nrem(incl_sub,ch,delta_ndx),2),3);

[r p] = corr(zpower_nrem,zpower_alpha_change);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
scatter(zpower_nrem, zpower_alpha_change,100,'k','LineWidth',3);
l = lsline;
l.LineWidth = 3;
l.Color = 'k';
xlabel('Delta power nrem');
ylabel('Alpha power change');
title(['r = ', num2str(r), ', p = ',num2str(p)]);
set(gca,'Fontsize',35);
box on
axis square
% xlim([8 12]);
% ylim([7 12]);
% saveas(fig,[Savefolder,'scatter_wake_e_wake_m_ch',EEG.chanlocs(ch).labels,'.svg']);

%%
% 
% IAPF_tonic_frontal = IAPF_tonic; 
% IAPF_phasic_frontal = IAPF_phasic; 
% 
% % IAPF_tonic_occipital = IAPF_tonic; 
% % IAPF_phasic_occipital = IAPF_phasic; 
% 
% for s = 1:length(IAPF_tonic_frontal)
% % figure
% plot([1 2],[IAPF_tonic_frontal(s) IAPF_tonic_occipital(s)]);
% hold on
% 
% end
% 
% ylim([7 11]);
% xlim([0 3]);


