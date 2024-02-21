clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

addpath(genpath('/user/HS301/m17462/matlab/ScientificColourMaps7'));

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/Figures/';

%% Load file

% load('/vol/research/nemo/datasets/RSN/data/analysis/phase_allsub/phase_allsub_mICA_avref_26-Aug-2022.mat');
% load('/vol/research/nemo/datasets/RSN/data/analysis/phase_allsub/phase_allsub_mICA_avref_24-Feb-2023.mat');
load('/vol/research/nemo/datasets/RSN/data/analysis/phase_allsub/phase_allsub_mICA_avref_alphathetafilt03-Aug-2023.mat');

condition = {'Alpha Phase 0' 'Alpha Phase 90' 'Alpha Phase 180' 'Alpha Phase 270' ...
    'Theta Phase 0' 'Theta Phase 90' 'Theta Phase 180' 'Theta Phase 270'};

target_phase = [0 90 180 270 0 90 180 270];

%% topoplot layout 

close all

% Computes electrode locations 

cfg = []; 
cfg.layout = 'EEG1005.lay';
   
layout = ft_prepare_layout(cfg);

figure
hold on
ylim([min(layout.pos(:,2)) max(layout.pos(:,2))])
xlim([min(layout.pos(:,1)) max(layout.pos(:,1))])

load('/user/HS301/m17462/matlab/Scripts/RSN/preprocessing/sleep/chans.mat');

for i = 1:128

%scatter(layout.pos(i,1),layout.pos(i,2))

clear indx
indx = strcmp(chans(i),layout.label);
text(layout.pos(indx,1),layout.pos(indx,2),chans(i))
layout2.pos(i,:) = layout.pos(indx,:);
layout2.pos(i,:) = layout.pos(indx,:);
layout2.width(i,:) = layout.width(indx,:);
layout2.height(i,:) = layout.height(indx,:);
layout2.label{i} = chans{i};

end

layout2.mask = layout.mask;
layout2.outline = layout.outline;

% Scalp Current Density - create electrode positions

[~,ftpath] = ft_version; 
elec = ft_read_sens(strcat(ftpath, '/template/electrode/standard_1005.elc' )); 

clear elec2

for i = 1:128

%scatter(layout.pos(i,1),layout.pos(i,2))

clear indx

indx = strcmp(chans(i),elec.label);

%text(layout.pos(indx,1),layout.pos(indx,2),chans(i))

elec2.chanpos(i,:) = elec.chanpos(indx,:);
elec2.elecpos(i,:) = elec.elecpos(indx,:);
elec2.label(i,:) = elec.label(indx,:);


end

elec2.chantype = elec.chantype;
elec2.unit = elec.unit;
elec2.chanunit = elec.chanunit;
elec2.type = elec.type;

%

clc

cfg = [];
cfg.layout = layout2;
cfg.method = 'triangulation';
cfg.compress = 'yes';
cfg.feedback = 'yes';%'yes';
neighbours = ft_prepare_neighbours(cfg);

%% Alpha stim alphafilt - Plot r - averaged across conditions and significant electrodes for all conditions

incl_sub = setdiff(1:19,[12]);     


figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r(incl_sub,:,1:4),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 0.8],'gridscale',300);
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r_alphafilt(incl_sub,:,1:4),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','clim',[0 0.8],'gridscale',300);
hold on

for con = 1:4

% nonan_sub = intersect(find(~isnan(m(:,1,con)) == 1),incl_sub);
nonan_sub = intersect(find(~isnan(m_alphafilt(:,1,con)) == 1),incl_sub);

for ch = 1:128
[pval(con,ch) vval(con,ch)] = circ_vtest(m_alphafilt(nonan_sub,ch,con), deg2rad(target_phase(con)), r_alphafilt(nonan_sub,ch,con));
m_dev_target(con,ch) = circ_mean(m_alphafilt(nonan_sub,ch,con)-deg2rad(target_phase(con)));
std_dev_target(con,ch) = circ_std(m_alphafilt(nonan_sub,ch,con)-deg2rad(target_phase(con)));

end

%sgtitle('Bandpower Change From Off')

end

sig_ch = find(sum(pval(1:4,:) < 0.05/128,1) == 4);

% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% cb = colorbar;

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',20);
axis off
axis square

% hold on
% scatter(layout2.pos(mask_bands_alpha(band,:),1),layout2.pos(mask_bands_alpha(band,:),2),5,'w','filled');

hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(sig_ch,1),layout2.pos(sig_ch,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
% end

% print([Savefolder,'topo_resultant_allcon'],'-dpng');
% saveas(gcf,[Savefolder,'Figure3_alphastim_topo_phase_accuracy_',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);
clear sig_ch

ch = 2
m_dev_target_alpha = rad2deg(circ_mean(m_dev_target(1:4,ch)))
std_dev_target_alpha = rad2deg(circ_mean(std_dev_target(1:4,ch)))

p_val_alpha = pval(1:4,ch)
v_val_alpha = vval(1:4,ch)

%% Theta stim alphafilt

incl_sub = setdiff(1:19,[12]);     


figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r(incl_sub,:,1:4),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 0.8],'gridscale',300);
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r_alphafilt(incl_sub,:,5:8),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','clim',[0 0.8],'gridscale',300);
hold on

for con = 5:8

% nonan_sub = intersect(find(~isnan(m(:,1,con)) == 1),incl_sub);
nonan_sub = intersect(find(~isnan(m_alphafilt(:,1,con)) == 1),incl_sub);

for ch = 1:128
[pval(con,ch) vval(con,ch)] = circ_vtest(m_alphafilt(nonan_sub,ch,con), deg2rad(target_phase(con)), r_alphafilt(nonan_sub,ch,con));

end

%sgtitle('Bandpower Change From Off')

end

sig_ch = find(sum(pval(5:8,:) < 0.05/128,1) == 4);

% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% cb = colorbar;

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',20);
axis off
axis square

% hold on
% scatter(layout2.pos(mask_bands_alpha(band,:),1),layout2.pos(mask_bands_alpha(band,:),2),5,'w','filled');

hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(sig_ch,1),layout2.pos(sig_ch,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
axis off
axis square

% end

% print([Savefolder,'topo_resultant_allcon'],'-dpng');
% saveas(gcf,[Savefolder,'Suppl_Figure3_thetastim_alphafilt_topo_phase_accuracy_',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);


%% Theta stim thetafilt - Plot r - averaged across conditions and significant electrodes for all conditions

incl_sub = setdiff(1:19,[12]);     


figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r(incl_sub,:,5:8),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 0.8],'gridscale',300);
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r_thetafilt(incl_sub,:,5:8),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','clim',[0 0.8],'gridscale',300);
hold on

clear pval vval
for con = 5:8

% nonan_sub = intersect(find(~isnan(m(:,1,con)) == 1),incl_sub);
nonan_sub = intersect(find(~isnan(m_thetafilt(:,1,con)) == 1),incl_sub);

for ch = 1:128
[pval(con,ch) vval(con,ch)] = circ_vtest(m_thetafilt(nonan_sub,ch,con), deg2rad(target_phase(con)), r_thetafilt(nonan_sub,ch,con));

m_dev_target(con,ch) = circ_mean(m_thetafilt(nonan_sub,ch,con)-deg2rad(target_phase(con)));
std_dev_target(con,ch) = circ_std(m_thetafilt(nonan_sub,ch,con)-deg2rad(target_phase(con)));

end

%sgtitle('Bandpower Change From Off')

end

sig_ch = find(sum(pval(5:8,:) < 0.05/128,1) == 4);

% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% cb = colorbar;

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',20);
axis off
axis square

% hold on
% scatter(layout2.pos(mask_bands_alpha(band,:),1),layout2.pos(mask_bands_alpha(band,:),2),5,'w','filled');

hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(sig_ch,1),layout2.pos(sig_ch,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
axis off
axis square

% end

ch = 2;
m_dev_target_theta = rad2deg(circ_mean(m_dev_target(5:8,ch)))
std_dev_target_theta = rad2deg(circ_mean(std_dev_target(5:8,ch)))

p_val_theta = pval(5:8,ch)
v_val_theta = vval(5:8,ch)

% print([Savefolder,'topo_resultant_allcon'],'-dpng');
% saveas(gcf,[Savefolder,'Figure3_thetastim_topo_phase_accuracy_',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%% Alpha stim thetafilt 

figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r(incl_sub,:,5:8),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 0.8],'gridscale',300);
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),nanmean(nanmean(r_thetafilt(incl_sub,:,1:4),1),3),'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','clim',[0 0.8],'gridscale',300);
hold on

clear pval vval
for con = 1:4

% nonan_sub = intersect(find(~isnan(m(:,1,con)) == 1),incl_sub);
nonan_sub = intersect(find(~isnan(m_thetafilt(:,1,con)) == 1),incl_sub);

for ch = 1:128
[pval(con,ch) vval(con,ch)] = circ_vtest(m_thetafilt(nonan_sub,ch,con), deg2rad(target_phase(con)), r_thetafilt(nonan_sub,ch,con));

end

%sgtitle('Bandpower Change From Off')

end

sig_ch = find(sum(pval(1:4,:) < 0.05/128,1) == 4);

% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% cb = colorbar;

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',20);
axis off
axis square

% hold on
% scatter(layout2.pos(mask_bands_alpha(band,:),1),layout2.pos(mask_bands_alpha(band,:),2),5,'w','filled');

hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(sig_ch,1),layout2.pos(sig_ch,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
axis off
axis square

% end

% print([Savefolder,'topo_resultant_allcon'],'-dpng');
% saveas(gcf,[Savefolder,'alphastim_thetafilt_topo_phase_accuracy_',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);
