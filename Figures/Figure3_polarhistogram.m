clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/Figures/';

%% Load file

% load('/vol/research/nemo/datasets/RSN/data/analysis/phase_allsub/phase_allsub_mICA_avref_24-Feb-2023.mat');
load('/vol/research/nemo/datasets/RSN/data/analysis/phase_allsub/phase_allsub_mICA_avref_alphathetafilt03-Aug-2023.mat');

incl_sub = setdiff(1:19,12);

%% Alpha stim - alphafilt

% RSN_015 a bit different than others

ch = 2;

colors = linspecer(4);

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for con = 1:4
    
    for sub = 1:length(incl_sub)

%     polarhistogram(m(:,ch,con),100,'FaceColor',colors(con,:),'LineWidth',.75)  % Polar Histogram   
%     polarplot([0 m(incl_sub(sub),ch,con)],[0 r(incl_sub(sub),ch,con)],'Color',[colors(con,:) .4],'LineWidth',3)
    polarplot([0 m_alphafilt(incl_sub(sub),ch,con)],[0 r_alphafilt(incl_sub(sub),ch,con)],'Color',[colors(con,:) .4],'LineWidth',3)

    hold on
    % title({['RSN ',mICA_file(1).name(5:7)] ['Alpha (Fz data)']});
    
    end
    
%     polarplot([0 circ_mean(m(incl_sub,ch,con))],[0 nanmean(r(incl_sub,ch,con))],'Color',colors(con,:),'LineWidth',5)
    polarplot([0 circ_mean(m_alphafilt(incl_sub,ch,con))],[0 nanmean(r_alphafilt(incl_sub,ch,con))],'Color',colors(con,:),'LineWidth',5)
    rlim([0 0.9])
   
end

  set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
  
% print([Savefolder,'poloarplot_alpha_allsub'],'-dpng');
% saveas(gcf,[Savefolder,'Figure3_polarplot_alpha_allsub.svg']);

%% Theta stim - alphafilt

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for con = 5:8
    
    for sub = 1:length(incl_sub)

%     polarhistogram(m(:,ch,con),100,'FaceColor',colors(con,:),'LineWidth',.75)  % Polar Histogram   
%     polarplot([0 m(incl_sub(sub),ch,con)],[0 r(incl_sub(sub),ch,con)],'Color',[colors(con,:) .4],'LineWidth',3)
    polarplot([0 m_alphafilt(incl_sub(sub),ch,con)],[0 r_alphafilt(incl_sub(sub),ch,con)],'Color',[colors(con-4,:) .4],'LineWidth',3)

    hold on
    % title({['RSN ',mICA_file(1).name(5:7)] ['Alpha (Fz data)']});
    
    end
    
%     polarplot([0 circ_mean(m(incl_sub,ch,con))],[0 nanmean(r(incl_sub,ch,con))],'Color',colors(con,:),'LineWidth',5)
    polarplot([0 circ_mean(m_alphafilt(incl_sub,ch,con))],[0 nanmean(r_alphafilt(incl_sub,ch,con))],'Color',colors(con-4,:),'LineWidth',5)
    rlim([0 0.9])

end

set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);


% saveas(gcf,[Savefolder,'Suppl_Figure3_polarplot_thetastim_alphafilt_allsub.svg']);


%% Theta stim - thetafilt
ch = 2;

colors = linspecer(4);

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for con = 5:8
    
    for sub = 1:length(incl_sub)

%     polarhistogram(m(:,ch,con),100,'FaceColor',colors(con,:),'LineWidth',.75)  % Polar Histogram   
%     polarplot([0 m(incl_sub(sub),ch,con)],[0 r(incl_sub(sub),ch,con)],'Color',[colors(con-4,:) .4],'LineWidth',3)
    polarplot([0 m_thetafilt(incl_sub(sub),ch,con)],[0 r_thetafilt(incl_sub(sub),ch,con)],'Color',[colors(con-4,:) .4],'LineWidth',3)

    hold on
    % title({['RSN ',mICA_file(1).name(5:7)] ['Alpha (Fz data)']});
    
    end
    
    polarplot([0 circ_mean(m_thetafilt(incl_sub,ch,con))],[0 nanmean(r_thetafilt(incl_sub,ch,con))],'Color',colors(con-4,:),'LineWidth',5)
    rlim([0 0.9])

end

set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);


% print([Savefolder,'poloarplot_theta_allsub'],'-dpng');
% saveas(gcf,[Savefolder,'Figure3_polarplot_theta_allsub.svg']);

%% Alpha stim - thetafilt

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
for con = 1:4
    
    for sub = 1:length(incl_sub)

%     polarhistogram(m(:,ch,con),100,'FaceColor',colors(con,:),'LineWidth',.75)  % Polar Histogram   
%     polarplot([0 m(incl_sub(sub),ch,con)],[0 r(incl_sub(sub),ch,con)],'Color',[colors(con-4,:) .4],'LineWidth',3)
    polarplot([0 m_thetafilt(incl_sub(sub),ch,con)],[0 r_thetafilt(incl_sub(sub),ch,con)],'Color',[colors(con,:) .4],'LineWidth',3)

    hold on
    % title({['RSN ',mICA_file(1).name(5:7)] ['Alpha (Fz data)']});
    
    end
    
    polarplot([0 circ_mean(m_thetafilt(incl_sub,ch,con))],[0 nanmean(r_thetafilt(incl_sub,ch,con))],'Color',colors(con,:),'LineWidth',5)
    rlim([0 0.9])
end

set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);

% print([Savefolder,'poloarplot_theta_allsub'],'-dpng');
% saveas(gcf,[Savefolder,'Suppl_Figure3_polarplot_alphastim_thetafilt_allsub.svg']);

%%

m_resultant_alpha = squeeze(nanmean(nanmean(r_alphafilt(incl_sub,ch,1:4),1),3))
std_resultant_alpha = squeeze(nanmean(nanstd(r_alphafilt(incl_sub,ch,1:4),1),3))

m_resultant_theta = squeeze(nanmean(nanmean(r_thetafilt(incl_sub,ch,5:8),1),3))
std_resultant_theta = squeeze(nanmean(nanstd(r_thetafilt(incl_sub,ch,5:8),1),3))

