clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

addpath(genpath('/vol/research/nemo/software/automaticanalysis/'));
aa_ver5

addpath(genpath('/user/HS301/m17462/matlab/ScientificColourMaps7'));

band_name = {'1-4 Hz' '4-7 Hz' '7-12 Hz' '13-30 Hz'};

conditions = {'Peak'; 'Falling'; 'Trough'; 'Rising';
    'Peak'; 'Falling'; 'Trough'; 'Rising';}

load('/vol/research/nemo/datasets/RSN/data/analysis/connectivity/connectivity_allsub_25-Aug-2023.mat');

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/Figures/';

statsresult_folder = '/vol/research/nemo/datasets/RSN/data/analysis/connectivity/statsresult_sensor/';

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

%% alphastim lme PLV - condition

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
band = 1;

for band1 = 1:4
% band1 = 3;
% band1 = 2;    
    for band2 = band1:4
% band2 = 3;

% figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
subplot(3,5,band)
band = band+1;

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);

clear statsresult
PlotChans2 = [];

 load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300);
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Fvalue_condition(band1,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','gridscale',300);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'Figure5_lme_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end

% saveas(gcf,[Savefolder,'Suppl_Figure5_lme_alphastim_condition_PLV_allbands.svg'])



%% PLV Alpha stim on vs off (t-test)

incl_sub = setdiff(1:19,12);

colors = linspecer(4);

Y = -.2.*zeros(127,1);

% comps = NaN(127,127);
% comps(:,1) = [1:127];
% 
% close all
% 
% for i = 2:127
%     
%     comps(i:127,i) = i:127;
%     
% end

% for band1 = 2 %1:4
% band1 = 3;
band1 = 2;
%     for band2 = 3 %band1:4
band2 = 3;

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

% tileplot = tiledlayout(1,4);

clear statsresult
PlotChans2 = [];

 load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.mat']);
 
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

for cond = 1:4
 
%     figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

%     nexttile(cond)
    subplot(1,4,cond)
    
%     subplot(4,4,(band-1)*4+cond)
    hold on
    %
    ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Y,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','clim',[-1 1]);
  
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu')))
%cb = colorbar;
    %}
    
    for seed = 1:127 %[2 34 65 94]
       
%     tmp_comps = comps(~isnan(comps(:,seed)),seed);
    
    P =0;
    T = 0;
    T2 = 0;
    
    for ch = 1:length(Plotchans2) %1:127 %length(tmp_comps)
        
        chan = Plotchans2(ch);
        
         x = squeeze(PLV_on.band1{band1}.band2{band2}.cond{cond}(chan,seed,:));
         y = squeeze(PLV_off.band1{band1}.band2{band2}.cond{cond}(chan,seed,:)); 
        
%         x = squeeze(PLV_on.band{band}.cond{cond}(chan,seed,:));
%         y = squeeze(PLV_off.band{band}.cond{cond}(chan,seed,:)); 

%         x = squeeze(PLV_on.band{band}.cond{cond}(tmp_comps(chan),seed,:));  
%         y = squeeze(PLV_off.band{band}.cond{cond}(tmp_comps(chan),seed,:));  

        change = log10(x./y);
                
%         [h,p,ci,stats] = ttest(x(incl_sub),y(incl_sub));
        [h,p,ci,stats] = ttest(change(incl_sub));

        if ~isnan(stats.tstat)
            T2 = T2 + stats.tstat;
        end

        if p <=.01
          if ~isnan(stats.tstat)
            T = T+stats.tstat; 
          end
        end
          
            %P = P+1;
            %
            if p <.01
            if stats.tstat>0 
            
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','r','Linewidth',1);
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));%,.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
            
            else
                
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','b','Linewidth',1);
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));%.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
                
            end
            end
            %}
        
        %}
        
    end  
    
    if T>0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),1*abs(T),colors(2,:),'filled')
%         scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'r','filled')
    end
        
    if T<0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),1*abs(T),colors(1,:),'filled')
%      scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'b','filled')

    end
    
%     Ts.band{band}(seed,cond) = T;
%     Ts2.band{band}(seed,cond) = T2;

    Ts.band1{band1}.band2{band2}(seed,cond) = T;
    Ts2.band1{band1}.band2{band2}(seed,cond) = T2;

    end
    
set(gca,'FontSize',18);
axis off
axis square
%scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
%scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')

    %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')

% title([conditions{cond}])

% saveas(gcf,[Savefolder,'Figure5_ttest_',band_name{band1},'-',band_name{band2},'_',conditions{cond},'_cluster_alpha_dots_PLV.svg']);
% close
end

% sgtitle([band_name{band},' alpha PLV'])
% cd(figPath)

saveas(gcf,[Savefolder,'Figure5_ttest_',band_name{band1},'-',band_name{band2},'cluster_alpha_dots_PLV.svg']);

%     end
% 
% end

clear pl_


%% alphastim lme PLV - substage

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
band = 1;

for band1 = 1:4
% band1 = 3;
% band1 = 2;    
    for band2 = band1:4
% band2 = 3;

% figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
subplot(3,5,band)
band = band+1;

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);

clear statsresult
PlotChans2 = [];

 load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300);
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Fvalue_substage(band1,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','gridscale',300);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'Figure5_lme_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end

saveas(gcf,[Savefolder,'Suppl_Figure5_lme_alphastim_substage_PLV_allbands.svg'])

%% alphastim lme PLV - substage * condition

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
band = 1;

for band1 = 1:4
% band1 = 3;
% band1 = 2;    
    for band2 = band1:4
% band2 = 3;

% figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
subplot(3,5,band)
band = band+1;

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);

clear statsresult
PlotChans2 = [];

 load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300);
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Fvalue_substage_condition(band1,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','gridscale',300);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'Figure5_lme_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end

saveas(gcf,[Savefolder,'Suppl_Figure5_lme_alphastim_substage_condition_PLV_allbands.svg'])





%% alphastim lme PLI - condition

% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
band = 1;

% for band1 = 1:4
% band1 = 1;
band1 = 2;    
%     for band2 = band1:4
band2 = 3;

figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
% subplot(3,5,band)
band = band+1;

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);

clear statsresult
PlotChans2 = [];

 load([statsresult_folder,'statsresult_alphastim_PLI_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300);
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Fvalue_condition(band1,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','gridscale',300);
hold on
% scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

saveas(gcf,[Savefolder,'Figure5_lme_alphastim_PLI_',band_name{band1},'-',band_name{band2},'.svg'])

%     end

% end

% saveas(gcf,[Savefolder,'Suppl_Figure5_lme_alphastim_condition_PLI_allbands.svg'])
   
%% alphastim lme PLI - substage

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
band = 1;

for band1 = 1:4
% band1 = 3;
% band1 = 2;    
    for band2 = band1:4
% band2 = 3;

% figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
subplot(3,5,band)
band = band+1;

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);

clear statsresult
PlotChans2 = [];

 load([statsresult_folder,'statsresult_alphastim_PLI_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300);
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Fvalue_substage(band1,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','gridscale',300);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'Figure5_lme_alphastim_PLI_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end

saveas(gcf,[Savefolder,'Suppl_Figure5_lme_alphastim_substage_PLI_allbands.svg'])

%% alphastim lme PLI - substage * condition

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
band = 1;

for band1 = 1:4
% band1 = 3;
% band1 = 2;    
    for band2 = band1:4
% band2 = 3;

% figure('Renderer','painters','units','normalized','outerposition',[0 0 0.5 0.5])

% subplot(1,2,1)
subplot(3,5,band)
band = band+1;

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);

clear statsresult
PlotChans2 = [];

 load([statsresult_folder,'statsresult_alphastim_PLI_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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

ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),statsresult.Fvalue_substage_condition,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','gridscale',300);
% ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Fvalue_substage_condition(band1,:),'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','gridscale',300);
hold on
scatter(layout2.pos(2,1),layout2.pos(2,2),50,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),20,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'Figure5_lme_alphastim_PLI_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end

saveas(gcf,[Savefolder,'Suppl_Figure5_lme_alphastim_substage_condition_PLI_allbands.svg'])

