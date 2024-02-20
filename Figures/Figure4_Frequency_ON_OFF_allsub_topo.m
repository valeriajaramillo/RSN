clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));
addpath(genpath('/user/HS301/m17462/matlab/ScientificColourMaps7'));
addpath(genpath('/user/HS301/m17462/matlab/DataViz'));

% load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/freqalphatheta_allsub_27-Mar-2023');

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/freqalphatheta_allsub_23-Jun-2023.mat');


Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/Figures/';

% Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/';
% Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/topo_allsub_lapref/';


%% Average across on and off blocks and calculate change

on_block = 7:12;
off_block = 1:6;
          
conditions = {'Alpha Phase 0'; 'Alpha Phase 90'; 'Alpha Phase 180'; 'Alpha Phase 270';...
            'Theta Phase 0'; 'Theta Phase 90'; 'Theta Phase 180'; 'Theta Phase 270';};
        
band_name = {'7-12 Hz';'4-7 Hz'};

incl_sub = setdiff(1:19,12);

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

%% Alpha stim - Topo lme condition


close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 1;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_alphastim_alphafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_condition > statsresult.p95_clus_condition
      if length(statsresult.WhichCh_1_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Figure4_alphastim_topo_lme_condition_alphaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%%
% close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 2;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_alphastim_thetafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_condition > statsresult.p95_clus_condition
      if length(statsresult.WhichCh_1_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 8]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
% 
% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Figure4_alphastim_topo_lme_condition_thetaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%% Alpha stim - Topo lme substage

close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 1;%:15

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_alphastim_alphafreq.mat')
  
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_substage > statsresult.p95_clus_substage
      if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
      PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 9]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])    

saveas(gcf,[Savefolder,'Suppl_Figure4_alphastim_topo_lme_substage_alphaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%%
% close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 2;%:15
    
load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_alphastim_thetafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_substage > statsresult.p95_clus_substage
      if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
      PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);


ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 9]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])    


saveas(gcf,[Savefolder,'Suppl_Figure4_alphastim_topo_lme_substage_thetaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);


%% Alpha stim - Topo lme condition*substage

close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 1; %:15
load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_alphastim_alphafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
      if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
      PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 7]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])    

saveas(gcf,[Savefolder,'Suppl_Figure4_alphastim_topo_lme_substagecondition_alphaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%%
% close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 2; %:15

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_alphastim_thetafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
      if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
      PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 7]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])    

saveas(gcf,[Savefolder,'Suppl_Figure4_alphastim_topo_lme_substagecondition_thetaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);


%% Theta stim - Topo lme condition

close all
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 1;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_thetastim_alphafreq.mat')
 
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_condition > statsresult.p95_clus_condition
      if length(statsresult.WhichCh_1_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 8]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Figure4_thetastim_topo_lme_condition_alphaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%%

% close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 2;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_thetastim_thetafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_condition(band,:) < 0.05);
  if statsresult.clus_max_condition > statsresult.p95_clus_condition
      if length(statsresult.WhichCh_1_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_condition) > statsresult.p95_clus_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 7]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Figure4_thetastim_topo_lme_condition_thetaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%% Theta stim - Topo lme substage

close all
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 1;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_thetastim_alphafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);
  if statsresult.clus_max_substage > statsresult.p95_clus_substage
      if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
  
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 7]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% % PlotChans2 = find(pvalue_substage(band,:) < 0.05);
% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Suppl_Figure4_thetastim_topo_lme_substage_alphaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%%
% close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 2;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_thetastim_thetafreq.mat')

figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);
  if statsresult.clus_max_substage > statsresult.p95_clus_substage
      if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
  
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 7]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% % PlotChans2 = find(pvalue_substage(band,:) < 0.05);
% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Suppl_Figure4_thetastim_topo_lme_substage_thetaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);


%% Theta stim - Topo lme substage *condition

close all
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 1;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_thetastim_alphafreq.mat')
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);
  if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
      if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
  
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 7]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% % PlotChans2 = find(pvalue_substage(band,:) < 0.05);
% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Suppl_Figure4_thetastim_topo_lme_substage_condition_alphaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

%%
% close all
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

band = 2;

load('/vol/research/nemo/datasets/RSN/data/analysis/frequency_allsub/statsresult/statsresult_thetastim_thetafreq.mat')
   
figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(3,5,band)
% hold on

% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat_bands_alpha(band,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 35],'gridscale',300);  

load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
PlotChans2 = [];
% PlotChans2 = find(pvalue_substage(band,:) < 0.05);
  if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
      if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
  else
        PlotChans2 = [];
  end
  
ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300,'clim',[0 7]);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% % PlotChans2 = find(pvalue_substage(band,:) < 0.05);
% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% hold on
% topoplottest3(statsresult.Fvalue_substage_condition,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
%     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band}])   

saveas(gcf,[Savefolder,'Suppl_Figure4_thetastim_topo_lme_substage_condition_thetaband_N',num2str(length(incl_sub)),'_lapaz_colorbar.svg']);

