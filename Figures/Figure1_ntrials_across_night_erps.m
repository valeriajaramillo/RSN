clear all;
close all;

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/';
sub_Folderpath = dir([Folderpath,'RSN*']);

addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));
addpath(genpath('/user/HS301/m17462/matlab/DataViz'));

hours = 11;
nbins = 3;
dt = 1;
fs = 500;

incl_sub = setdiff(1:19,12);


Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/Figures/';


%% Load erp data

for s = 1:length(sub_Folderpath)
    
    display(sub_Folderpath(s).name);

    erp_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_erp_good.mat']);
    load([Folderpath,sub_Folderpath(s).name,filesep,erp_good_file(1).name],'ERP');
    
%     goodREM_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_goodREM.mat']);
%     load([Folderpath,sub_Folderpath(s).name,filesep,goodREM_file(1).name]);

    nm_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_nm_good_psd_allch.mat']);
    load([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file(1).name],'nm');

 %% ERPs
 
    for con = 1:size(ERP.trigs_con,2)
 
        ERP_trigs = ERP.trigs;
        ERP_trigs_good = ERP_trigs(ERP.trigs_state_good_all{con});
            
        ERPs_ntrials_all(s,con) = length(ERP.trigs_state_good_all{con});
        ERPs_ntrials_tonic(s,con) = length(ERP.trigs_state_good{1,con});
        ERPs_ntrials_phasic(s,con) = length(ERP.trigs_state_good{2,con});
        

      for b = 1:nbins
          bin_duration = 10/nbins;
          bin_samps = (b-1)*bin_duration*60*60*fs+1:b*bin_duration*60*60*fs;
          ERPs_ntrials_bins(s,con,b) = length(intersect(ERP_trigs_good,bin_samps));
          ERPs_perc_trials_bins(s,con,b) = ERPs_ntrials_bins(s,con,b)/length(ERP_trigs_good)*100;
      end

    clear ERP_trigs_good 
    
    end
    
    
      %% Neuromodulation
    
    for con = 1:size(nm.ON_start_good,2)
        
        ON_start = nm.ON_start_good{con};
        ON_start_good = ON_start(nm.good_trialndx{con});
%         ON_start_good_hours = ceil(ON_start_good/(fs*60*60));
        
        nm_ntrials_all(s,con) = nm.ntrials_good{con};
        nm_ntrials_phasic(s,con) = nm.ntrials_phasic{con};
        nm_ntrials_tonic(s,con) = nm.ntrials_tonic{con};
      
%         for h = 1:hours
%             ntrials_hours(s,con,h) = length(find(ON_start_good_hours == h));
%             perc_trials_hours(s,con,h) = ntrials_hours(s,con,h)/length(ON_start_good)*100;
%         end

        for b = 1:nbins
            bin_duration = 10/nbins;
            bin_samps = (b-1)*bin_duration*60*60*fs+1:b*bin_duration*60*60*fs;
            nm_ntrials_bins(s,con,b) = length(intersect(ON_start_good,bin_samps));
            nm_perc_trials_bins(s,con,b) = nm_ntrials_bins(s,con,b)/length(ON_start_good)*100;
            
        end
       
   
    end
    
    clear ON_start*
    
end  


%% 

ERPs_ntrials_allcon_all = sum(ERPs_ntrials_all(:,1:4),2);
ERPs_ntrials_allcon_tonic = sum(ERPs_ntrials_tonic(:,1:4),2);
ERPs_ntrials_allcon_phasic = sum(ERPs_ntrials_phasic(:,1:4),2);

ERPs_ntrials_allcon_all_bins = squeeze(sum(ERPs_ntrials_bins(:,1:4,:),2));

m_ERPs_ntrials_allcon_all = nanmean(ERPs_ntrials_allcon_all(incl_sub),1);
sd_ERPs_ntrials_allcon_all = nanstd(ERPs_ntrials_allcon_all(incl_sub),1);

m_ERPs_ntrials_bins = nanmean(ERPs_ntrials_allcon_all_bins(incl_sub,:));
sd_ERPs_ntrials_bins = nanstd(ERPs_ntrials_allcon_all_bins(incl_sub,:));

bins = {'All' '1st third' '2nd third' '3rd third'};

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

% saveas(fig,[Savefolder,'Figure1_Number_of_trials_AEPs.svg']);

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

% saveas(fig,[Savefolder,'Figure1_Number_of_trials_CLAS.svg']);



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



