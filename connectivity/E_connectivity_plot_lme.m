clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));

addpath(genpath('/vol/research/nemo/software/automaticanalysis/'));
aa_ver5

addpath(genpath('/user/HS301/m17462/matlab/ScientificColourMaps7'));

band_name = {'1-4 Hz' '4-7 Hz' '7-12 Hz' '13-30 Hz'};

conditions = {'Peak'; 'Falling'; 'Trough'; 'Rising';
    'Peak'; 'Falling'; 'Trough'; 'Rising';}

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/connectivity/';

statsresult_folder = '/vol/research/nemo/datasets/RSN/data/analysis/connectivity/';     

%% topoplot layout 

close all

% Computes electrode location

cfg = []; 
cfg.layout = 'EEG1005.lay';
   
layout = ft_prepare_layout(cfg);

figure
hold on
ylim([min(layout.pos(:,2)) max(layout.pos(:,2))])
xlim([min(layout.pos(:,1)) max(layout.pos(:,1))])

load('/user/HS301/m17462/matlab/Scripts/RSN/preprocessing/sleep/chans.mat');

for i = 1:127

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

for i = 1:127

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

% load('EEG_chanlocs.mat');

%% PLV alpha stim - condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'.mat']);
        
        
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
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end
   
%% PLV alpha stim - substage

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage > statsresult.p95_clus_substage
        if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end


%% PLV alpha stim - substage * condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
        if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end

%% PLI alpha stim - condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_alphastim_PLI_',band_name{band1},'.mat']);
        
        
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
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end
   
%% PLI alpha stim - substage

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_alphastim_PLI_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage > statsresult.p95_clus_substage
        if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end


%% PLI alpha stim - substage * condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_alphastim_PLI_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
        if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end


%% PLV theta stim - condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_thetastim_PLV_',band_name{band1},'.mat']);
        
        
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
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end
   
%% PLV theta stim - substage

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_thetastim_PLV_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage > statsresult.p95_clus_substage
        if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end


%% PLV theta stim - substage * condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_thetastim_PLV_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
        if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end

%% PLI theta stim - condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_thetastim_PLI_',band_name{band1},'.mat']);
        
        
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
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end
   
%% PLI theta stim - substage

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_thetastim_PLI_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage > statsresult.p95_clus_substage
        if length(statsresult.WhichCh_1_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage) > statsresult.p95_clus_substage
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage]; %find(pvalue_substage(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end


%% PLI theta stim - substage * condition

incl_sub = setdiff(1:19,12);

for band1 = 1:4
    
%     for band2 = band1:4
            

      load([statsresult_folder,'statsresult_thetastim_PLI_',band_name{band1},'.mat']);
        
        
      if statsresult.clus_max_substage_condition > statsresult.p95_clus_substage_condition
        if length(statsresult.WhichCh_1_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_1_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_2_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_2_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_3_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_3_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_4_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_4_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      if length(statsresult.WhichCh_5_max_substage_condition) > statsresult.p95_clus_substage_condition
        PlotChans2 = [PlotChans2 statsresult.WhichCh_5_max_substage_condition]; %find(pvalue_substage_condition(band,:) < 0.05);
      end
      Plotchans2 = unique(PlotChans2);  
      else
        PlotChans2 = [];
      end
      
      if ~isempty(PlotChans2)
          
%         close all
        figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
        hold on
          
        ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
        'interplim','mask','gridscale',300);
        hold on
        % scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
        % 
        % hold on
        % scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

        scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

        hold on
        scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);
        
        load('lapaz.mat');
        colormap(lapaz)
        colorbar
        set(gca,'FontSize',18);
        axis off
        axis square
        title([band_name{band1}])  

      end

       PlotChans2_all{band1} = PlotChans2;
%      PlotChans2_all{band1,band2} = PlotChans2;
      
      
%     end
    
end