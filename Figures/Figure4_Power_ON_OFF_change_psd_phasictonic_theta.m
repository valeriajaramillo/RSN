clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/';
sub_Folderpath = dir([Folderpath,'RSN*']);

% Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/';
Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/Figures/';

% load('/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/psd_allsub_mICA_avref_26-Aug-2022.mat');
load('/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/psd_allsub_mICA_avref_12-Mar-2023.mat');
% load('/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/psd_allsub_mICA_avref_rejtrials_26-Feb-2023');

% load('/user/HS301/m17462/matlab/Scripts/RSN/analyses/psd/topo/Cluster_el_9_10_Hz.mat');

load('/user/HS301/m17462/matlab/Scripts/RSN/analyses/psd/topo/EEG_chanlocs.mat');

%%

load('/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/statsresult/statsresult_alpha_7-8 Hz.mat');
alpha_cluster_el1 = statsresult.WhichCh_1_max_condition; 
clear statsresult

load('/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/statsresult/statsresult_alpha_10-11 Hz.mat');
alpha_cluster_el2 = statsresult.WhichCh_1_max_condition; 
clear statsresult

alpha_cluster_combined = unique([alpha_cluster_el1,alpha_cluster_el2]);

load('/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/statsresult/statsresult_theta_7-8 Hz.mat');
theta_cluster_el = statsresult.WhichCh_1_max_condition; 
clear statsresult


conditions = {'Peak'; 'Falling'; 'Trough'; 'Rising';
    'Peak'; 'Falling'; 'Trough'; 'Rising';}

%% Average across on and off blocks and calculate change - theta cluster

mpsd_con_ch_stage = mpsd_con_ch; % mpsd_con_ch: sub x ch x con x bins x ep
% mpsd_con_ch_stage = mpsd_con_ch_phasic;

% mpsd_con_ch_stage(3,:,7,:,:) = NaN; % should be NaN as this participant did not have this condition

on_block = 7:12;
off_block = 1:6;

ch = theta_cluster_el; %[2 34 65 94]; % Fz, AFz, AFF1h, AFF2h
% ch = theta_cluster_el;
% ch = 16:18; % O1, Oz, Oz
% ch = cluster_el

psd_ON_block = squeeze(nanmean(mpsd_con_ch_stage(:,:,:,:,on_block),5)); % average across on block 
psd_OFF_block = squeeze(nanmean(mpsd_con_ch_stage(:,:,:,:,off_block),5));

psd_ON = squeeze(nanmean(psd_ON_block(:,ch,:,:),2)); % average across channels
psd_OFF = squeeze(nanmean(psd_OFF_block(:,ch,:,:),2));

psd_ON_block_phasic = squeeze(nanmean(mpsd_con_ch_phasic(:,:,:,:,on_block),5)); % psd_ON: sub x ch x con x bins
psd_OFF_block_phasic = squeeze(nanmean(mpsd_con_ch_phasic(:,:,:,:,off_block),5));

psd_ON_block_tonic = squeeze(nanmean(mpsd_con_ch_tonic(:,:,:,:,on_block),5)); % psd_ON: sub x ch x con x bins
psd_OFF_block_tonic = squeeze(nanmean(mpsd_con_ch_tonic(:,:,:,:,off_block),5));

psd_ON_phasic = squeeze(nanmean(psd_ON_block_phasic(:,ch,:,:),2));
psd_OFF_phasic = squeeze(nanmean(psd_OFF_block_phasic(:,ch,:,:),2));

psd_ON_tonic = squeeze(nanmean(psd_ON_block_tonic(:,ch,:,:),2));
psd_OFF_tonic = squeeze(nanmean(psd_OFF_block_tonic(:,ch,:,:),2));

psd_ONOFF_change = (psd_ON - psd_OFF)./psd_OFF *100;
log_psd_ONOFF_change = log10(psd_ON./psd_OFF);
log_npsd_ONOFF_change = log_psd_ONOFF_change-nanmean(log_psd_ONOFF_change(:,5:8,:),2);

psd_ONOFF_change_phasic = (psd_ON_phasic - psd_OFF_phasic)./psd_OFF_phasic *100;
log_psd_ONOFF_change_phasic = log10(psd_ON_phasic./psd_OFF_phasic);
log_npsd_ONOFF_change_phasic = log_psd_ONOFF_change_phasic-nanmean(log_psd_ONOFF_change_phasic(:,5:8,:),2);

psd_ONOFF_change_tonic = (psd_ON_tonic - psd_OFF_tonic)./psd_OFF_tonic *100;
log_psd_ONOFF_change_tonic = log10(psd_ON_tonic./psd_OFF_tonic);
log_npsd_ONOFF_change_tonic = log_psd_ONOFF_change_tonic-nanmean(log_psd_ONOFF_change_tonic(:,5:8,:),2);

npsd_ONOFF_change = psd_ONOFF_change-nanmean(psd_ONOFF_change(:,5:8,:),2);

incl_sub = setdiff(1:19,12);
mpsd_ONOFF_change = squeeze(nanmean(psd_ONOFF_change(incl_sub,:,:),1));
sem_psd_ONOFF_change = squeeze(nanstd(psd_ONOFF_change(incl_sub,:,:))./sqrt(size(psd_ONOFF_change(incl_sub,:,:),1)));

m_npsd_ONOFF_change = squeeze(nanmean(npsd_ONOFF_change(incl_sub,:,:),1));
sem_npsd_ONOFF_change = squeeze(nanstd(npsd_ONOFF_change(incl_sub,:,:))./sqrt(size(npsd_ONOFF_change(incl_sub,:,:),1)));

%% Compare On to Off using a t-test

% two sample
% for con = 1:size(psd_ON,2)
%     
%     for b = 1:size(psd_ON,3)
%     
%     [h(con,b) p(con,b) ci stats] = ttest(psd_ON(incl_sub,con,b),psd_OFF(incl_sub,con,b));
%     tvalue(con,b) = stats.tstat;
%     cohens_d(con,b) = nanmean(psd_ON(incl_sub,con,b)-psd_OFF(incl_sub,con,b))/nanstd(psd_ON(incl_sub,con,b)-psd_OFF(incl_sub,con,b));
%     
%     end
% 
% end

% one sample
for con = 1:size(psd_ON,2)
    
    for b = 1:size(psd_ON,3)
    
%     [h(con,b) p_abs(con,b) ci stats] = ttest(psd_ONOFF_change(incl_sub,con,b));
    [h(con,b) p_abs(con,b) ci stats] = ttest(log_psd_ONOFF_change(incl_sub,con,b));
    tvalue_abs(con,b) = stats.tstat;
    cohens_d_abs(con,b) = nanmean(psd_ON(incl_sub,con,b)- psd_OFF(incl_sub,con,b))/nanstd(psd_ON(incl_sub,con,b)-psd_OFF(incl_sub,con,b));
%     cohens_d_abs2(con,b) = (nanmean(psd_ON(incl_sub,con,b))- nanmean(psd_OFF(incl_sub,con,b)))/nanstd(psd_ON(incl_sub,con,b)-psd_OFF(incl_sub,con,b)); % gives the same results as line above  
    
    end

end

%% Compare normalized On to Off using a t-test

% one sample
for con = 1:size(psd_ON,2)
    
    for b = 1:size(psd_ON,3)
    
    [h(con,b) p_norm(con,b) ci stats] = ttest(npsd_ONOFF_change(incl_sub,con,b));
    tvalue_norm(con,b) = stats.tstat;
    cohens_d_norm(con,b) = nanmean(psd_ON(incl_sub,con,b)-psd_OFF(incl_sub,con,b))/nanstd(psd_ON(incl_sub,con,b)-psd_OFF(incl_sub,con,b));
    
    end

end

%% Compare conditions using a lme - theta

% clear sub
% 
% for s = 1:size(sub_Folderpath,1)
%    
%     sub{s} = sub_Folderpath(s).name;
%     
% end
% 
% 
% for b = 1:size(psd_ONOFF_change,3)
%    
%     psd_ON_bin = psd_ON(incl_sub,:,b);
%     psd_OFF_bin = psd_OFF(incl_sub,:,b);
% 
%     table_allcon_theta = [];
%     
%     for con = 5:8
%         
%         sub_table = vertcat(sub(incl_sub)',sub(incl_sub)');
%         cond = repmat(con,length(incl_sub)*2,1);
%         wind = vertcat(repmat(1,length(incl_sub),1),repmat(0,length(incl_sub),1));
%         psd_bin_con = vertcat(psd_ON_bin(:,con),psd_OFF_bin(:,con));
%         table_bin_con_theta = table(sub_table,cond,wind,psd_bin_con,'VariableNames',{'sub','condition','window','power'});
%         table_allcon_theta = vertcat(table_allcon_theta,table_bin_con_theta);
%         
%     end
%     
%    lme = fitlme(table_allcon_theta,'power ~ condition * window + (1|sub)','FitMethod','REML','DummyVarCoding','effects');
% %    [pVal_alpha(b),F_alhpa(b),DF1_alpha(b),DF2_alpha(b)] = coefTest(lme);
% %    [beta,betanames,stats] = fixedEffects(lme);
% %    p_con_alpha(b) = stats.pValue(2);
% %    p_win_alpha(b) = stats.pValue(3);
% %    p_con_win_alpha(b) = stats.pValue(4);
% 
%    stats = anova(lme);
%    p_con_theta(b) = stats.pValue(2);
%    p_win_theta(b) = stats.pValue(3);
%    p_con_win_theta(b) = stats.pValue(4);
% 
%     
% end


clear sub

for s = 1:size(sub_Folderpath,1)
   
    sub{s} = sub_Folderpath(s).name;
    
end


for b = 1:size(psd_ONOFF_change,3)
   
    psd_change_phasic_bin = log_psd_ONOFF_change_phasic(incl_sub,:,b);
    psd_change_tonic_bin = log_psd_ONOFF_change_tonic(incl_sub,:,b);

    table_allcon_theta = [];
    
    for con = 5:8
        
        sub_table = vertcat(sub(incl_sub)',sub(incl_sub)');
        cond = repmat(con,length(incl_sub)*2,1);
        substage = vertcat(repmat(1,length(incl_sub),1),repmat(2,length(incl_sub),1));
        psd_change_bin_con = vertcat(psd_change_phasic_bin(:,con),psd_change_tonic_bin(:,con));
        table_bin_con_theta = table(sub_table,cond,substage,psd_change_bin_con,'VariableNames',{'sub','condition','substage','power_change'});
        table_allcon_theta = vertcat(table_allcon_theta,table_bin_con_theta);
        
    end
    
    table_allcon_theta.substage = categorical(table_allcon_theta.substage);
    table_allcon_theta.condition = categorical(table_allcon_theta.condition);
    table_allcon_theta.sub = categorical(table_allcon_theta.sub);
    
   lme = fitlme(table_allcon_theta,'power_change ~ substage * condition + (1|sub)','FitMethod','REML','DummyVarCoding','effects');
%    [pVal_theta(b),F_alhpa(b),DF1_theta(b),DF2_theta(b)] = coefTest(lme);
%    [beta,betanames,stats] = fixedEffects(lme);
%    p_con_theta(b) = stats.pValue(2);
%    p_win_theta(b) = stats.pValue(3);
%    p_con_win_theta(b) = stats.pValue(4);

   stats = anova(lme);
   p_con_theta(b) = stats.pValue(2);
   p_substage_theta(b) = stats.pValue(3);
   p_con_substage_theta(b) = stats.pValue(4);

    
end

%% Compare conditions using a lme - theta normalized


clear sub

for s = 1:size(sub_Folderpath,1)
   
    sub{s} = sub_Folderpath(s).name;
    
end


for b = 1:size(npsd_ONOFF_change,3)
   
    npsd_change_phasic_bin = log_npsd_ONOFF_change_phasic(incl_sub,:,b);
    npsd_change_tonic_bin = log_npsd_ONOFF_change_tonic(incl_sub,:,b);

    table_allcon_theta = [];
    
    for con = 5:8
        
        sub_table = vertcat(sub(incl_sub)',sub(incl_sub)');
        cond = repmat(con,length(incl_sub)*2,1);
        substage = vertcat(repmat(1,length(incl_sub),1),repmat(2,length(incl_sub),1));
        npsd_change_bin_con = vertcat(npsd_change_phasic_bin(:,con),npsd_change_tonic_bin(:,con));
        table_bin_con_theta = table(sub_table,cond,substage,npsd_change_bin_con,'VariableNames',{'sub','condition','substage','power_change'});
        table_allcon_theta = vertcat(table_allcon_theta,table_bin_con_theta);
        
    end
    
    table_allcon_theta.substage = categorical(table_allcon_theta.substage);
    table_allcon_theta.condition = categorical(table_allcon_theta.condition);
    table_allcon_theta.sub = categorical(table_allcon_theta.sub);
    
   lme = fitlme(table_allcon_theta,'power_change ~ substage * condition + (1|sub)','FitMethod','REML','DummyVarCoding','effects');
%    [pVal_theta(b),F_alhpa(b),DF1_theta(b),DF2_theta(b)] = coefTest(lme);
%    [beta,betanames,stats] = fixedEffects(lme);
%    p_con_theta(b) = stats.pValue(2);
%    p_win_theta(b) = stats.pValue(3);
%    p_con_win_theta(b) = stats.pValue(4);

   stats = anova(lme);
   p_con_theta_norm(b) = stats.pValue(2);
   p_substage_theta_norm(b) = stats.pValue(3);
   p_con_substage_theta_norm(b) = stats.pValue(4);

    
end

%% Plot power spectrum ON-OFF change for Theta stimulation

colors = linspecer(4);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for con = 5:8
    
%     shadedErrorBar(f,mpsd_ONOFF_change(con,:),sem_psd_ONOFF_change(con,:),'lineProps',{'Color',colors(con-4,:),'LineWidth',1},'patchSaturation',.3);
    
    plot(f,mpsd_ONOFF_change(con,:),'Color',colors(con-4,:),'LineWidth',7);

    hold on
    
%     shadedErrorBar(f,mnpsd_ONOFF_change(con,:),sem_psd_ONOFF_change(con,:),'lineProps',{'Color',colors(con-4,:),'LineWidth',1},'patchSaturation',.3);
%     hold on
    errorbar(15+(con-4)/2,0,nanmean(sem_psd_ONOFF_change(con,:),2),'Color',colors(con-4,:),'LineWidth',3);
 
end



for con = 5:8
 
%  sig_bins = find(p(con,:) <= 0.05);
%  plot(f(sig_bins),ones(length(sig_bins),1)*-25-1*con,'*','Color',colors(con-4,:));
%  hold on

end

con = 5;
sig_bins_ttest = find(p_abs(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-18,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

con = 6;
sig_bins_ttest = find(p_abs(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-20,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

con = 7;
sig_bins_ttest = find(p_abs(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-22,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

con = 8;
sig_bins_ttest = find(p_abs(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-24,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

% sig_bins_lme = find(p_con_win_theta <= 0.05);
% plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-28,'*','Color','k');

sig_bins_lme = find(p_con_theta <= 0.05);
plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-28,'square','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);

% sig_bins_lme = find(p_substage_theta <= 0.05);
% plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-30,'*','Color',[0.5 0.5 0.5]);
% 
% sig_bins_lme = find(p_con_substage_theta <= 0.05);
% plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-32,'*','Color','r');

% legend(nm.condition{5:8});
xlim([0 18]);
ylim([-30 30])
xlabel('Frequency (Hz)');
ylabel('Power change (%)');
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
xticks(0:2:18);
% set(gca, 'XTickMode', 'manual', 'XTick', [0:2:40], 'xlim', [0,40])
% set(gca,'Fontsize',15);
yline(0,'LineWidth',2);
yline(-26,'LineWidth',2);
  area([4.5 7.5],[30 30],-30,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.1,'EdgeColor','none','LineStyle','none')
%   xline(4.5,'LineWidth',1,'LineStyle','--')
  xline(4.5,'LineWidth',1,'LineStyle','--')
  xline(7.5,'LineWidth',1,'LineStyle','--')

box off
axis square

saveas(fig,[Savefolder,'Figure4_psd_change_theta_ttest_statsresult_cluster.svg']);


%% 7-8 Hz

% bins = find(f >= 7 & f <= 8);
% 
% fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% violins = violinPlot([nanmean(psd_ONOFF_change(incl_sub,5,bins),3) nanmean(psd_ONOFF_change(incl_sub,6,bins),3) nanmean(psd_ONOFF_change(incl_sub,7,bins),3) nanmean(psd_ONOFF_change(incl_sub,8,bins),3)]);
% % ylabel('Frequency change (Hz)');
% ylabel('Power change (%)');
% set(gca,'Fontsize',35);
% xticklabels(conditions(5:8));
% xtickangle(45);
% box on
% axis square
% title('7-8 Hz');
% % ylim([-0.2 0.35])
% % yticks(-0.2:0.1:0.3);
% xlim([0 5])
% 
% 
%   for cond = 1:4
%       violins(cond).ViolinColor = colors(cond,:); 
%       violins(cond).ShowMean = 1;
%       violins(cond).ShowData = 1;
%   end
% 
% saveas(fig,[Savefolder,'Figure4_violinplots_thetapower_7_8_ttest_statsresult_cluster.svg']);
% 
% [h p_theta_peak] = ttest(nanmean(psd_ONOFF_change(incl_sub,5,bins),3))
% [h p_theta_falling] = ttest(nanmean(psd_ONOFF_change(incl_sub,6,bins),3))
% [h p_theta_trough] = ttest(nanmean(psd_ONOFF_change(incl_sub,7,bins),3))
% [h p_theta_rising] = ttest(nanmean(psd_ONOFF_change(incl_sub,8,bins),3))


%% Plot normalized power spectrum ON-OFF change for Theta stimulation

colors = linspecer(4);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for con = 5:8
    
%     shadedErrorBar(f,mpsd_ONOFF_change(con,:),sem_psd_ONOFF_change(con,:),'lineProps',{'Color',colors(con-4,:),'LineWidth',1},'patchSaturation',.3);
    
    plot(f,m_npsd_ONOFF_change(con,:),'Color',colors(con-4,:),'LineWidth',7);

    hold on
    
%     shadedErrorBar(f,mnpsd_ONOFF_change(con,:),sem_psd_ONOFF_change(con,:),'lineProps',{'Color',colors(con-4,:),'LineWidth',1},'patchSaturation',.3);
%     hold on
%     errorbar(9+(con-4)/4,0,nanmean(sem_npsd_ONOFF_change(con,:),2),'Color',colors(con-4,:),'LineWidth',3);
 
end



for con = 5:8
 
%  sig_bins = find(p(con,:) <= 0.05);
%  plot(f(sig_bins),ones(length(sig_bins),1)*-25-1*con,'*','Color',colors(con-4,:));
%  hold on

end

con = 5;
sig_bins_ttest = find(p_norm(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-10,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

con = 6;
sig_bins_ttest = find(p_norm(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-11,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

con = 7;
sig_bins_ttest = find(p_norm(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-12,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

con = 8;
sig_bins_ttest = find(p_norm(con,:) <= 0.05);
plot(f(sig_bins_ttest),ones(length(sig_bins_ttest),1)*-13,'square','MarkerFaceColor',colors(con-4,:),'MarkerEdgeColor',colors(con-4,:),'MarkerSize',5);

% sig_bins_lme = find(p_con_win_theta <= 0.05);
% plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-28,'*','Color','k');

sig_bins_lme = find(p_con_theta_norm <= 0.05);
plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-15,'square','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);

% sig_bins_lme = find(p_substage_theta <= 0.05);
% plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-30,'*','Color',[0.5 0.5 0.5]);
% 
% sig_bins_lme = find(p_con_substage_theta <= 0.05);
% plot(f(sig_bins_lme),ones(length(sig_bins_lme),1)*-32,'*','Color','r');

xlim([3.5 12.5]);
ylim([-16 10])
xticks(4:1:12);
% xline(10,'LineWidth',2);
% plot([0 40],[0 0],'Color',[0.5 0.5 0.5]);
% hold on
% plot([0 40],[-10 -10],'Color',[0.5 0.5 0.5]);
% hold on
% plot([0 40],[10 10],'Color',[0.5 0.5 0.5]);
xlabel('Frequency (Hz)');
ylabel('Normalized power change (%)');
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
yline(0,'LineWidth',2);
yline(-14,'LineWidth',2);
% yline(-24);
  area([4.5 7.5],[10 10],-16,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',.1,'EdgeColor','none','LineStyle','none')
%   xline(4.5,'LineWidth',1,'LineStyle','--')
  xline(4.5,'LineWidth',1,'LineStyle','--')
  xline(7.5,'LineWidth',1,'LineStyle','--')
  
  
box off
axis square

saveas(fig,[Savefolder,'Figure4_npsd_change_theta_ttest_statsresult_cluster.svg']);


%% 6.5-8.7 Hz

f_sig_bins_lme = f(sig_bins_lme);
freq_ndx = find(f >= 6.5 & f <= 8.7);

psd_on_ch = squeeze(nanmean(mpsd_con_ch_stage(incl_sub,ch,:,:,1:12),2));
m_psd_off_ch = nanmean(squeeze(nanmean(mpsd_con_ch_stage(incl_sub,ch,:,:,off_block),2)),4); % mean across cluster and off block
perc_change = (psd_on_ch./m_psd_off_ch-1)*100;
log_perc_change = log10(psd_on_ch./m_psd_off_ch);

perc_change_freq = squeeze(nanmean(perc_change(:,:,freq_ndx,:),3));
log_perc_change_freq = squeeze(nanmean(log_perc_change(:,:,freq_ndx,:),3));

movmean_sec = 2;

cond = 5;
nperc_change_cond5 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond5 = movmean(nperc_change_cond1,[movmean_sec],2);
m_nperc_change_cond5 = nanmean(nperc_change_cond1,1);
sem_nperc_change_cond5 = nanstd(nperc_change_cond1,1)./sqrt(size(nperc_change_cond1,1));
log_perc_change_freq_cond5 = squeeze(log_perc_change_freq(:,cond,:)-nanmean(log_perc_change_freq(:,5:8,:),2));

cond = 6;
nperc_change_cond6 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond6 = movmean(nperc_change_cond6,[movmean_sec],2);
m_nperc_change_cond6 = nanmean(nperc_change_cond6,1);
sem_nperc_change_cond6 = nanstd(nperc_change_cond6,1)./sqrt(size(nperc_change_cond6,1));
log_perc_change_freq_cond6 = squeeze(log_perc_change_freq(:,cond,:)-nanmean(log_perc_change_freq(:,5:8,:),2));

cond = 7;
nperc_change_cond7 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond7= movmean(nperc_change_cond7,[movmean_sec],2);
m_nperc_change_cond7 = nanmean(nperc_change_cond7,1);
sem_nperc_change_cond7 = nanstd(nperc_change_cond7,1)./sqrt(size(nperc_change_cond7,1));
log_perc_change_freq_cond7 = squeeze(log_perc_change_freq(:,cond,:)-nanmean(log_perc_change_freq(:,5:8,:),2));

cond = 8;
nperc_change_cond8 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond8 = movmean(nperc_change_cond8,[movmean_sec],2);
m_nperc_change_cond8 = nanmean(nperc_change_cond8,1);
sem_nperc_change_cond8 = nanstd(nperc_change_cond8,1)./sqrt(size(nperc_change_cond8,1));
log_perc_change_freq_cond8 = squeeze(log_perc_change_freq(:,cond,:)-nanmean(log_perc_change_freq(:,5:8,:),2));


psd_tf_con5_minus_con7 = nperc_change_cond5 - nperc_change_cond7;
mpsd_tf_con5_minus_con7 = nanmean(psd_tf_con5_minus_con7,1);
sem_psd_tf_con5_minus_con7 = nanstd(psd_tf_con5_minus_con7)./sqrt(size(psd_tf_con5_minus_con7,1));
log_psd_tf_con5_minus_con7 = log_perc_change_freq_cond5 - log_perc_change_freq_cond7;

psd_tf_con6_minus_con8 = nperc_change_cond6 - nperc_change_cond8;
mpsd_tf_con6_minus_con8 = nanmean(psd_tf_con6_minus_con8,1);
sem_psd_tf_con6_minus_con8 = nanstd(psd_tf_con6_minus_con8)./sqrt(size(psd_tf_con6_minus_con8,1));
log_psd_tf_con6_minus_con8 = log_perc_change_freq_cond6 - log_perc_change_freq_cond8;


clear p_freq1_con5_con7 p_freq1_con6_con8
for sec = 1:12
    
  [h p_freq1_con5_con7(sec)] = ttest(nperc_change_cond5(:,sec),nperc_change_cond7(:,sec))
  [h p_freq1_con6_con8(sec)] = ttest(nperc_change_cond6(:,sec),nperc_change_cond8(:,sec))

end

sig_times_con5_con7 = find(p_freq1_con5_con7 <= 0.05);
trend_times_con5_con7 = find(p_freq1_con5_con7 <= 0.1 & p_freq1_con5_con7 > 0.05);

sig_times_con6_con8 = find(p_freq1_con6_con8 <= 0.05);
trend_times_con6_con8 = find(p_freq1_con6_con8 <= 0.1 & p_freq1_con6_con8 > 0.05);

%%
fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar([0.5:1:11.5],m_nperc_change_cond5,sem_nperc_change_cond5,'lineProps',{'Color',colors(1,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond1,'Color',colors(1,:),'LineWidth',7);
hold on
shadedErrorBar([0.5:1:11.5],m_nperc_change_cond7,sem_nperc_change_cond7,'lineProps',{'Color',colors(3,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond3,'Color',colors(3,:),'LineWidth',7);
hold on
% shadedErrorBar([0.5:1:11.5],mpsd_tf_con1_minus_con3,sem_psd_tf_con1_minus_con3,'lineProps',{'Color','k','LineWidth',1},'patchSaturation',.3);
plot(sig_times_con5_con7-0.5,ones(length(sig_times_con5_con7),1,1)*18,'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
hold on
plot(trend_times_con5_con7-0.5,ones(length(trend_times_con5_con7),1,1)*18,'+','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
xlim([0.5 11.5]);
xline(6,'LineStyle','--','LineWidth',2);
ylim([-20 20]);
xticks(0:2:12);
xlabel('Time (s)');
ylabel('Normalized power change (%)');
title('6.5-8.7 Hz');
set(gca,'Fontsize',35,'LineWidth',3);
box off
axis square

% saveas(fig,[Savefolder,'Figure4_npsd_change_theta_ttest_statsresult_cluster_time_freq1_peak_trough.svg']);

%%
fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar([0.5:1:11.5],m_nperc_change_cond6,sem_nperc_change_cond6,'lineProps',{'Color',colors(2,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond2,'Color',colors(2,:),'LineWidth',7);
hold on
shadedErrorBar([0.5:1:11.5],m_nperc_change_cond8,sem_nperc_change_cond8,'lineProps',{'Color',colors(4,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond4,'Color',colors(4,:),'LineWidth',7);
hold on
% shadedErrorBar([0.5:1:11.5],mpsd_tf_con2_minus_con4,sem_psd_tf_con2_minus_con4,'lineProps',{'Color','k','LineWidth',1},'patchSaturation',.3);
hold on
plot(sig_times_con6_con8-0.5,ones(length(sig_times_con6_con8),1,1)*18,'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
hold on
plot(trend_times_con6_con8-0.5,ones(length(trend_times_con6_con8),1,1)*18,'+','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
xlim([0.5 11.5]);
xline(6,'LineStyle','--');
xlim([0.5 11.5]);
xline(6,'LineStyle','--','LineWidth',2);
ylim([-20 20]);
xticks(0:2:12);
xlabel('Time (s)');
ylabel('Normalized power change (%)');
title('6.5-8.7 Hz');
set(gca,'Fontsize',35,'LineWidth',3);
box off
axis square

% saveas(fig,[Savefolder,'Figure4_npsd_change_theta_ttest_statsresult_cluster_time_freq1_falling_rising.svg']);


%% 10.1-11.5 Hz

f_sig_bins_lme = f(sig_bins_lme);
freq_ndx = find(f >= 10.1 & f <= 11.5);

psd_on_ch = squeeze(nanmean(mpsd_con_ch_stage(incl_sub,ch,:,:,1:12),2));
m_psd_off_ch = nanmean(squeeze(nanmean(mpsd_con_ch_stage(incl_sub,ch,:,:,off_block),2)),4); % mean across cluster and off block
perc_change = (psd_on_ch./m_psd_off_ch-1)*100;
log_perc_change = log10(psd_on_ch./m_psd_off_ch);

perc_change_freq = squeeze(nanmean(perc_change(:,:,freq_ndx,:),3));
log_perc_change_freq = squeeze(nanmean(log_perc_change(:,:,freq_ndx,:),3));

movmean_sec = 2;

cond = 5;
nperc_change_cond5 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond5 = movmean(nperc_change_cond1,[movmean_sec],2);
m_nperc_change_cond5 = nanmean(nperc_change_cond1,1);
sem_nperc_change_cond5 = nanstd(nperc_change_cond1,1)./sqrt(size(nperc_change_cond1,1));
log_perc_change_freq_cond5 = squeeze(log_perc_change_freq(:,cond,:));%-nanmean(log_perc_change_freq(:,5:8,:),2));

cond = 6;
nperc_change_cond6 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond6 = movmean(nperc_change_cond6,[movmean_sec],2);
m_nperc_change_cond6 = nanmean(nperc_change_cond6,1);
sem_nperc_change_cond6 = nanstd(nperc_change_cond6,1)./sqrt(size(nperc_change_cond6,1));
log_perc_change_freq_cond6 = squeeze(log_perc_change_freq(:,cond,:));%-nanmean(log_perc_change_freq(:,5:8,:),2));

cond = 7;
nperc_change_cond7 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond7= movmean(nperc_change_cond7,[movmean_sec],2);
m_nperc_change_cond7 = nanmean(nperc_change_cond7,1);
sem_nperc_change_cond7 = nanstd(nperc_change_cond7,1)./sqrt(size(nperc_change_cond7,1));
log_perc_change_freq_cond7 = squeeze(log_perc_change_freq(:,cond,:));%-nanmean(log_perc_change_freq(:,5:8,:),2));

cond = 8;
nperc_change_cond8 = squeeze(perc_change_freq(:,cond,:)-nanmean(perc_change_freq(:,5:8,:),2));
nperc_change_cond8 = movmean(nperc_change_cond8,[movmean_sec],2);
m_nperc_change_cond8 = nanmean(nperc_change_cond8,1);
sem_nperc_change_cond8 = nanstd(nperc_change_cond8,1)./sqrt(size(nperc_change_cond8,1));
log_perc_change_freq_cond8 = squeeze(log_perc_change_freq(:,cond,:));%-nanmean(log_perc_change_freq(:,5:8,:),2));


psd_tf_con5_minus_con7 = nperc_change_cond5 - nperc_change_cond7;
mpsd_tf_con5_minus_con7 = nanmean(psd_tf_con5_minus_con7,1);
sem_psd_tf_con5_minus_con7 = nanstd(psd_tf_con5_minus_con7)./sqrt(size(psd_tf_con5_minus_con7,1));
log_psd_tf_con5_minus_con7 = log_perc_change_freq_cond5 - log_perc_change_freq_cond7;

psd_tf_con6_minus_con8 = nperc_change_cond6 - nperc_change_cond8;
mpsd_tf_con6_minus_con8 = nanmean(psd_tf_con6_minus_con8,1);
sem_psd_tf_con6_minus_con8 = nanstd(psd_tf_con6_minus_con8)./sqrt(size(psd_tf_con6_minus_con8,1));
log_psd_tf_con6_minus_con8 = log_perc_change_freq_cond6 - log_perc_change_freq_cond8;


clear p_freq1_con5_con7 p_freq1_con6_con8
for sec = 1:12
    
  [h p_freq1_con5_con7(sec)] = ttest(nperc_change_cond5(:,sec),nperc_change_cond7(:,sec))
  [h p_freq1_con6_con8(sec)] = ttest(nperc_change_cond6(:,sec),nperc_change_cond8(:,sec))

end

sig_times_con5_con7 = find(p_freq1_con5_con7 <= 0.05);
trend_times_con5_con7 = find(p_freq1_con5_con7 <= 0.1 & p_freq1_con5_con7 > 0.05);

sig_times_con6_con8 = find(p_freq1_con6_con8 <= 0.05);
trend_times_con6_con8 = find(p_freq1_con6_con8 <= 0.1 & p_freq1_con6_con8 > 0.05);

%%
fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar([0.5:1:11.5],m_nperc_change_cond5,sem_nperc_change_cond5,'lineProps',{'Color',colors(1,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond1,'Color',colors(1,:),'LineWidth',7);
hold on
shadedErrorBar([0.5:1:11.5],m_nperc_change_cond7,sem_nperc_change_cond7,'lineProps',{'Color',colors(3,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond3,'Color',colors(3,:),'LineWidth',7);
hold on
% shadedErrorBar([0.5:1:11.5],mpsd_tf_con1_minus_con3,sem_psd_tf_con1_minus_con3,'lineProps',{'Color','k','LineWidth',1},'patchSaturation',.3);
plot(sig_times_con5_con7-0.5,ones(length(sig_times_con5_con7),1,1)*18,'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
hold on
plot(trend_times_con5_con7-0.5,ones(length(trend_times_con5_con7),1,1)*18,'+','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
xlim([0.5 11.5]);
xline(6,'LineStyle','--','LineWidth',2);
ylim([-20 20]);
xticks(0:2:12);
xlabel('Time (s)');
ylabel('Normalized power change (%)');
title('10.1-11.5 Hz');
set(gca,'Fontsize',35,'LineWidth',3);
box off
axis square

% saveas(fig,[Savefolder,'Figure4_npsd_change_theta_ttest_statsresult_cluster_time_freq2_peak_trough.svg']);

%%
fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

shadedErrorBar([0.5:1:11.5],m_nperc_change_cond6,sem_nperc_change_cond6,'lineProps',{'Color',colors(2,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond2,'Color',colors(2,:),'LineWidth',7);
hold on
shadedErrorBar([0.5:1:11.5],m_nperc_change_cond8,sem_nperc_change_cond8,'lineProps',{'Color',colors(4,:),'LineWidth',1},'patchSaturation',.3);
% plot([0.5:1:11.5],m_nperc_change_cond4,'Color',colors(4,:),'LineWidth',7);
hold on
% shadedErrorBar([0.5:1:11.5],mpsd_tf_con2_minus_con4,sem_psd_tf_con2_minus_con4,'lineProps',{'Color','k','LineWidth',1},'patchSaturation',.3);
hold on
plot(sig_times_con6_con8-0.5,ones(length(sig_times_con6_con8),1,1)*18,'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
hold on
plot(trend_times_con6_con8-0.5,ones(length(trend_times_con6_con8),1,1)*18,'+','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20,'LineWidth',2);
xlim([0.5 11.5]);
xline(6,'LineStyle','--');
xlim([0.5 11.5]);
xline(6,'LineStyle','--','LineWidth',2);
ylim([-20 20]);
xticks(0:2:12);
xlabel('Time (s)');
ylabel('Normalized power change (%)');
title('10.1-11.5 Hz');
set(gca,'Fontsize',35,'LineWidth',3);
box off
axis square

% saveas(fig,[Savefolder,'Figure4_npsd_change_theta_ttest_statsresult_cluster_time_freq2_falling_rising.svg']);
