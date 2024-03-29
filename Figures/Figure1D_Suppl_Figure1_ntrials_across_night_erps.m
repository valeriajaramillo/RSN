clear all;
close all;

addpath(genpath('S:\projects\RSN\matlab\matlab\Henry\useful_functions')); % contains linspecer function, circular statistics toolbox functions, echt function, shadedErrorBar function, see README on where to find this
addpath(genpath('S:\projects\RSN\matlab\matlab\DataViz'));  % Dataviz toolbox, see README on where to find this

load('D:\Valeria\RSN\data\for_sharing\data_to_make_figures\ERP_ntrials_allsub.mat')
load('D:\Valeria\RSN\data\for_sharing\data_to_make_figures\nm_ntrials_allsub.mat')

incl_sub = setdiff(1:19,12);

bins = {'All' '1st third' '2nd third' '3rd third'};

Savefolder = 'D:\Valeria\RSN\data\for_sharing\data_to_make_figures\Figures\';

%%

ERPs_ntrials_allcon_all = sum(ERPs_ntrials_all(:,1:4),2);
ERPs_ntrials_allcon_tonic = sum(ERPs_ntrials_tonic(:,1:4),2);
ERPs_ntrials_allcon_phasic = sum(ERPs_ntrials_phasic(:,1:4),2);

ERPs_ntrials_allcon_all_bins = squeeze(sum(ERPs_ntrials_bins(:,1:4,:),2));

m_ERPs_ntrials_allcon_all = nanmean(ERPs_ntrials_allcon_all(incl_sub),1);
sd_ERPs_ntrials_allcon_all = nanstd(ERPs_ntrials_allcon_all(incl_sub),1);

m_ERPs_ntrials_bins = nanmean(ERPs_ntrials_allcon_all_bins(incl_sub,:));
sd_ERPs_ntrials_bins = nanstd(ERPs_ntrials_allcon_all_bins(incl_sub,:));


ntrials_ERPs_table = table(bins',vertcat(m_ERPs_ntrials_allcon_all,m_ERPs_ntrials_bins'),vertcat(sd_ERPs_ntrials_allcon_all,sd_ERPs_ntrials_bins'),'VariableNames',{'Bins' 'Mean' 'SD'});
% writetable(ntrials_ERPs_table,[Savefolder,'ntrials_ERPs.xlsx']);

m_nm_ntrials_allcon_all = nanmean(nm_ntrials_all(incl_sub,:),1);
sd_nm_ntrials_allcon_all = nanstd(nm_ntrials_all(incl_sub,:),1);

m_nm_ntrials_bins = squeeze(nanmean(nm_ntrials_bins(incl_sub,:,:),1));
sd_nm_ntrials_bins = squeeze(nanstd(nm_ntrials_bins(incl_sub,:,:),1));

ntrials_nm_table = table(bins',vertcat(m_ERPs_ntrials_allcon_all,m_ERPs_ntrials_bins'),vertcat(sd_ERPs_ntrials_allcon_all,sd_ERPs_ntrials_bins'),'VariableNames',{'Bins' 'Mean' 'SD'});

%% AEP boxplot

colors = linspecer(5);

ERPs_data = ERPs_ntrials_allcon_all(incl_sub);
ERPs_phasic = ERPs_ntrials_allcon_phasic(incl_sub);

condition_names = {'AEP'};
group_names = {'Random'};

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

% h = daboxplot(ERPs_data,'groups',ERPs_inx);

% non-filled boxplots and cutomized medians
h = daboxplot(ERPs_phasic,'outliers',0,'outsymbol','kx',...
    'xtlabels', condition_names,'fill',0,'color',colors(5,:),...
    'whiskers',0,'boxalpha',0.1);
ylabel('Number of trials');
xl = xlim; xlim([xl(1), xl(2)+1]);     % make more space for the legend
set(h.md,'LineWidth',1.5); % customize median lines
hold on
j = daboxplot(ERPs_data,'outliers',1,'outsymbol','kx',...
    'xtlabels', condition_names,'fill',1,'color',colors(5,:),...
    'whiskers',1,'boxalpha',0.7);
ylabel('Number of trials');
xl = xlim; xlim([xl(1), xl(2)+1]);     % make more space for the legend
set(j.md,'LineWidth',1.5); % customize median lines
hold on
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
xtickangle(45);
box off
axis square

saveas(fig,[Savefolder,'Figure1D_Number_of_trials_AEPs.svg']);

%% Alpha and Theta boxplot

alpha_allcon = vertcat(nm_ntrials_all(incl_sub,1),nm_ntrials_all(incl_sub,2),nm_ntrials_all(incl_sub,3),nm_ntrials_all(incl_sub,4));
theta_allcon = vertcat(nm_ntrials_all(incl_sub,5),nm_ntrials_all(incl_sub,6),nm_ntrials_all(incl_sub,7),nm_ntrials_all(incl_sub,8));

alpha_allcon_phasic = vertcat(nm_ntrials_phasic(incl_sub,1),nm_ntrials_phasic(incl_sub,2),nm_ntrials_phasic(incl_sub,3),nm_ntrials_phasic(incl_sub,4));
theta_allcon_phasic = vertcat(nm_ntrials_phasic(incl_sub,5),nm_ntrials_phasic(incl_sub,6),nm_ntrials_phasic(incl_sub,7),nm_ntrials_phasic(incl_sub,8));

alpha_allcon_tonic = vertcat(nm_ntrials_tonic(incl_sub,1),nm_ntrials_tonic(incl_sub,2),nm_ntrials_tonic(incl_sub,3),nm_ntrials_tonic(incl_sub,4));
theta_allcon_tonic = vertcat(nm_ntrials_tonic(incl_sub,5),nm_ntrials_tonic(incl_sub,6),nm_ntrials_tonic(incl_sub,7),nm_ntrials_tonic(incl_sub,8));

all_data = horzcat(alpha_allcon,theta_allcon);
all_data_phasic = horzcat(alpha_allcon_phasic,theta_allcon_phasic);
all_data_tonic = horzcat(alpha_allcon_tonic,theta_allcon_tonic);

group_inx = vertcat(repmat(1,length(nm_ntrials_all(incl_sub,5)),1),repmat(2,length(nm_ntrials_all(incl_sub,5)),1),repmat(3,length(nm_ntrials_all(incl_sub,5)),1),repmat(4,length(nm_ntrials_all(incl_sub,5)),1));

condition_names = {'Alpha' 'Theta'};
group_names = {'Peak' 'Falling' 'Trough' 'Rising'};

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

% h = daboxplot(ERPs_data,'groups',ERPs_inx);

% non-filled boxplots and cutomized medians
h = daboxplot(all_data_phasic,'groups',group_inx,'outliers',0,'outsymbol','kx',...
    'xtlabels', condition_names,'fill',0,'color',colors,...
    'whiskers',0,'boxalpha',0.1);
ylabel('Number of trials');
xl = xlim; xlim([xl(1), xl(2)+1]);     % make more space for the legend
set(h.md,'LineWidth',1.5); % customize median lines
hold on
j = daboxplot(all_data,'groups',group_inx,'outliers',1,'outsymbol','kx',...
    'xtlabels', condition_names,'fill',1,'color',colors,...
    'whiskers',1,'boxalpha',0.7);
ylabel('Number of trials');
xl = xlim; xlim([xl(1), xl(2)+1]);     % make more space for the legend
set(j.md,'LineWidth',1.5); % customize median lines
set(j.ot,'SizeData',50); % customize median lines
hold on

% hold on
% dabarplot(all_data,'groups',group_inx,'fill',1,'color',colors,'baralpha',0.7,'xtlabels', condition_names);
% hold on
% dabarplot(all_data_phasic,'groups',group_inx,'fill',1,'color',colors,'xtlabels', condition_names);
% ylabel('Number of trials');

% filled boxplots, different color scheme, non-jittered scatter underneath
% h = daboxplot(all_data,'groups',group_inx,'outsymbol','k+',...
%     'xtlabels', condition_names,'legend',group_names,'color',colors,...
%     'whiskers',1,'scatter',0,'jitter',0,'scattersize',13);
% ylabel('Number of trials');
% xl = xlim; xlim([xl(1), xl(2)+1]);    % make more space for the legend
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
xtickangle(45);
box off
axis square

saveas(fig,[Savefolder,'Figure1D_Number_of_trials_CLAS.svg']);



%% Alpha thirds

alpha_1st = vertcat(nm_ntrials_bins(incl_sub,1,1),nm_ntrials_bins(incl_sub,2,1),nm_ntrials_bins(incl_sub,3,1),nm_ntrials_bins(incl_sub,4,1));
alpha_2nd = vertcat(nm_ntrials_bins(incl_sub,1,2),nm_ntrials_bins(incl_sub,2,2),nm_ntrials_bins(incl_sub,3,2),nm_ntrials_bins(incl_sub,4,2));
alpha_3rd = vertcat(nm_ntrials_bins(incl_sub,1,3),nm_ntrials_bins(incl_sub,2,3),nm_ntrials_bins(incl_sub,3,3),nm_ntrials_bins(incl_sub,4,3));

all_data = horzcat(alpha_1st,alpha_2nd,alpha_3rd);

group_inx = vertcat(repmat(1,length(nm_ntrials_all(incl_sub,5)),1),repmat(2,length(nm_ntrials_all(incl_sub,5)),1),repmat(3,length(nm_ntrials_all(incl_sub,5)),1),repmat(4,length(nm_ntrials_all(incl_sub,5)),1));

condition_names = {'1st' '2nd' '3rd'};
group_names = {'Peak' 'Falling' 'Trough' 'Rising'};

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

% h = daboxplot(ERPs_data,'groups',ERPs_inx);

% non-filled boxplots and cutomized medians
j = daboxplot(all_data,'groups',group_inx,'outliers',1,'outsymbol','kx',...
    'xtlabels', condition_names,'fill',1,'color',colors,...
    'whiskers',1,'boxalpha',0.7);
ylabel('Number of trials');
xl = xlim; xlim([xl(1), xl(2)+1]);     % make more space for the legend
set(j.md,'LineWidth',1.5); % customize median lines
set(j.ot,'SizeData',50); % customize median lines
hold on

% hold on
% dabarplot(all_data,'groups',group_inx,'fill',1,'color',colors,'baralpha',0.7,'xtlabels', condition_names);
% hold on
% dabarplot(all_data_phasic,'groups',group_inx,'fill',1,'color',colors,'xtlabels', condition_names);
% ylabel('Number of trials');

% filled boxplots, different color scheme, non-jittered scatter underneath
% h = daboxplot(all_data,'groups',group_inx,'outsymbol','k+',...
%     'xtlabels', condition_names,'legend',group_names,'color',colors,...
%     'whiskers',1,'scatter',0,'jitter',0,'scattersize',13);
% ylabel('Number of trials');
% xl = xlim; xlim([xl(1), xl(2)+1]);    % make more space for the legend
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
xtickangle(45);
box off
axis square
ylim([0 60]);

saveas(fig,[Savefolder,'Suppl_Figure1_Number_of_trials_CLAS_alpha_thirds.svg']);

%% Theta thirds

theta_1st = vertcat(nm_ntrials_bins(incl_sub,5,1),nm_ntrials_bins(incl_sub,6,1),nm_ntrials_bins(incl_sub,7,1),nm_ntrials_bins(incl_sub,8,1));
theta_2nd = vertcat(nm_ntrials_bins(incl_sub,5,2),nm_ntrials_bins(incl_sub,6,2),nm_ntrials_bins(incl_sub,7,2),nm_ntrials_bins(incl_sub,8,2));
theta_3rd = vertcat(nm_ntrials_bins(incl_sub,5,3),nm_ntrials_bins(incl_sub,6,3),nm_ntrials_bins(incl_sub,7,3),nm_ntrials_bins(incl_sub,8,3));

all_data = horzcat(theta_1st,theta_2nd,theta_3rd);

group_inx = vertcat(repmat(1,length(nm_ntrials_all(incl_sub,5)),1),repmat(2,length(nm_ntrials_all(incl_sub,5)),1),repmat(3,length(nm_ntrials_all(incl_sub,5)),1),repmat(4,length(nm_ntrials_all(incl_sub,5)),1));

condition_names = {'1st' '2nd' '3rd'};
group_names = {'Peak' 'Falling' 'Trough' 'Rising'};

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

% h = daboxplot(ERPs_data,'groups',ERPs_inx);

% non-filled boxplots and cutomized medians
j = daboxplot(all_data,'groups',group_inx,'outliers',1,'outsymbol','kx',...
    'xtlabels', condition_names,'fill',1,'color',colors,...
    'whiskers',1,'boxalpha',0.7);
ylabel('Number of trials');
xl = xlim; xlim([xl(1), xl(2)+1]);     % make more space for the legend
set(j.md,'LineWidth',1.5); % customize median lines
set(j.ot,'SizeData',50); % customize median lines
hold on

% hold on
% dabarplot(all_data,'groups',group_inx,'fill',1,'color',colors,'baralpha',0.7,'xtlabels', condition_names);
% hold on
% dabarplot(all_data_phasic,'groups',group_inx,'fill',1,'color',colors,'xtlabels', condition_names);
% ylabel('Number of trials');

% filled boxplots, different color scheme, non-jittered scatter underneath
% h = daboxplot(all_data,'groups',group_inx,'outsymbol','k+',...
%     'xtlabels', condition_names,'legend',group_names,'color',colors,...
%     'whiskers',1,'scatter',0,'jitter',0,'scattersize',13);
% ylabel('Number of trials');
% xl = xlim; xlim([xl(1), xl(2)+1]);    % make more space for the legend
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
xtickangle(45);
box off
axis square
ylim([0 60]);

saveas(fig,[Savefolder,'Suppl_Figure1_Number_of_trials_CLAS_theta_thirds.svg']);



