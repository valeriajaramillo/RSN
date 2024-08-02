clear all;
close all;

addpath(genpath('/users/nemo/software/eeglab')); % eeglab toolbox, see README on where to find this
addpath(genpath('/users/nemo/software/Henry/useful_functions')); % contains linspecer function, circular statistics toolbox functions, echt function, shadedErrorBar function, see README on where to find this

Savefolder = '/parallel_scratch/nemo/RSN/analysis/analysis/Figures/';

%% Load file

load('/parallel_scratch/nemo/RSN/analysis/analysis/phase_allsub/phase_allsub_mICA_avref_alphathetafilt_notecht_19-Jun-2024.mat');

incl_sub = setdiff(1:19,12);

%% Alpha stim - alphafilt

% RSN_015 a bit different than others

ch = 2;

colors = linspecer(4);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
fig.WindowState = 'maximized';

for con = 1:4
    
    for sub = 1:length(incl_sub)

    polarplot([0 m_alphafilt(incl_sub(sub),ch,con)],[0 r_alphafilt(incl_sub(sub),ch,con)],'Color',[colors(con,:) .4],'LineWidth',3)
    hold on
    
    end
    
    polarplot([0 circ_mean(m_alphafilt(incl_sub,ch,con))],[0 nanmean(r_alphafilt(incl_sub,ch,con))],'Color',colors(con,:),'LineWidth',5)
    rlim([0 0.5])
   
end

  set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
  
saveas(gcf,[Savefolder,'Figure3D_polarplot_alpha_allsub_notecht.svg']);

%% Theta stim - alphafilt

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
fig.WindowState = 'maximized';

for con = 5:8
    
    for sub = 1:length(incl_sub)

    polarplot([0 m_alphafilt(incl_sub(sub),ch,con)],[0 r_alphafilt(incl_sub(sub),ch,con)],'Color',[colors(con-4,:) .4],'LineWidth',3)

    hold on
    
    end
    
    polarplot([0 circ_mean(m_alphafilt(incl_sub,ch,con))],[0 nanmean(r_alphafilt(incl_sub,ch,con))],'Color',colors(con-4,:),'LineWidth',5)
    rlim([0 0.5])

end

set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);


saveas(gcf,[Savefolder,'Suppl_Figure4C_polarplot_thetastim_alphafilt_allsub_notecht.svg']);


%% Theta stim - thetafilt
ch = 2;

colors = linspecer(4);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
fig.WindowState = 'maximized';

for con = 5:8
    
    for sub = 1:length(incl_sub)

    polarplot([0 m_thetafilt(incl_sub(sub),ch,con)],[0 r_thetafilt(incl_sub(sub),ch,con)],'Color',[colors(con-4,:) .4],'LineWidth',3)

    hold on
    
    end
    
    polarplot([0 circ_mean(m_thetafilt(incl_sub,ch,con))],[0 nanmean(r_thetafilt(incl_sub,ch,con))],'Color',colors(con-4,:),'LineWidth',5)
    rlim([0 0.5])

end

set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);


saveas(gcf,[Savefolder,'Figure3I_polarplot_theta_allsub_notecht.svg']);

%% Alpha stim - thetafilt

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
fig.WindowState = 'maximized';

for con = 1:4
    
    for sub = 1:length(incl_sub)

    polarplot([0 m_thetafilt(incl_sub(sub),ch,con)],[0 r_thetafilt(incl_sub(sub),ch,con)],'Color',[colors(con,:) .4],'LineWidth',3)
    hold on
    
    end
    
    polarplot([0 circ_mean(m_thetafilt(incl_sub,ch,con))],[0 nanmean(r_thetafilt(incl_sub,ch,con))],'Color',colors(con,:),'LineWidth',5)
    rlim([0 0.5])
end

set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);

saveas(gcf,[Savefolder,'Suppl_Figure4A_polarplot_alphastim_thetafilt_allsub_notecht.svg']);

%%

m_resultant_alpha = squeeze(nanmean(nanmean(r_alphafilt(incl_sub,ch,1:4),1),3))
std_resultant_alpha = squeeze(nanmean(nanstd(r_alphafilt(incl_sub,ch,1:4),1),3))

m_resultant_theta = squeeze(nanmean(nanmean(r_thetafilt(incl_sub,ch,5:8),1),3))
std_resultant_theta = squeeze(nanmean(nanstd(r_thetafilt(incl_sub,ch,5:8),1),3))

