clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

addpath(genpath('/user/HS301/m17462/matlab/ScientificColourMaps7'));

% load('/vol/research/nemo/datasets/RSN/data/analysis/connectivity/connectivity_allsub_15-Jun-2023.mat');

load('/vol/research/nemo/datasets/RSN/data/analysis/connectivity/connectivity_allsub_25-Aug-2023.mat');

band_name = {'1-4 Hz' '4-7 Hz' '7-12 Hz' '13-30 Hz' '8-12 Hz'};
          
% conditions = {'Alpha Phase 0'; 'Alpha Phase 90'; 'Alpha Phase 180'; 'Alpha Phase 270';...
%             'Theta Phase 0'; 'Theta Phase 90'; 'Theta Phase 180'; 'Theta Phase 270';}

conditions = {'Peak'; 'Falling'; 'Trough'; 'Rising';
    'Peak'; 'Falling'; 'Trough'; 'Rising';}

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

%% PLV off

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

for band = 1:4

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

for cond = 1:4
    
    subplot(1,4,cond)
    hold on
    %
    ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Y,'mask',layout2.mask,'outline',layout2.outline, ...
    'interplim','mask','clim',[-1 1]);

ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu')))
%cb = colorbar;
    %}
    
    for seed = [2 34 65 94]
       
%     tmp_comps = comps(~isnan(comps(:,seed)),seed);
    
    P =0;
    T = 0;
    T2 = 0;
    
    for chan = 1:127 %length(tmp_comps)
        
%         x = squeeze(PLV_on.band{band}.cond{cond}(chan,seed,:));
        y = squeeze(PLV_off.band{band}.cond{cond}(chan,seed,:)); 

%         x = squeeze(PLV_on.band{band}.cond{cond}(tmp_comps(chan),seed,:));  
%         y = squeeze(PLV_off.band{band}.cond{cond}(tmp_comps(chan),seed,:));    
                
        plv_thresh = 0.5;
        [h,p,ci,stats] = ttest(y(incl_sub),plv_thresh,'Tail','right');
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
%             if stats.tstat>0 
            
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','r','Linewidth',1);
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',[0.5 0.5 0.5],'Linewidth',.5);%,.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
            
%             else
                
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','b','Linewidth',1);
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(1,:),'Linewidth',.5);%.5*abs(stats.tstat));
% %             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));
% 
%             pl_.Color(4) = .2;
                
%             end
            end
            %}
        
        %}
        
    end  
    
    if T>0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.05*abs(T),[0.5 0.5 0.5],'filled')
%         scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'r','filled')
    end
        
    if T<0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.05*abs(T),[0.5 0.5 0.5],'filled')
%      scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'b','filled')

    end
Ts.band{band}(seed,cond) = T;
                Ts2.band{band}(seed,cond) = T2;

    end
    
axis off
axis square
%scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
%scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')

    %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')

title([conditions{cond}])
end

sgtitle([band_name{band},' alpha OFF PLV > ',num2str(plv_thresh)])
% cd(figPath)
saveas(gcf,[Savefolder,band_name{band},'allch_alpha_OFF_dots_PLV_',num2str(plv_thresh),'.svg']);


end

clear pl_



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

for band1 = 2 %1:4

    for band2 = 3 %band1:4
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

for cond = 1:4
    
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
    
    for chan = 1:127 %length(tmp_comps)
        
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
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(2,:),'Linewidth',.5);%,.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
            
            else
                
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','b','Linewidth',1);
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(1,:),'Linewidth',.5);%.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
                
            end
            end
            %}
        
        %}
        
    end  
    
    if T>0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(2,:),'filled')
%         scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'r','filled')
    end
        
    if T<0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(1,:),'filled')
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

title([conditions{cond}])
end

% sgtitle([band_name{band},' alpha PLV'])
% cd(figPath)

saveas(gcf,[Savefolder,'ttest_',band_name{band1},'-',band_name{band2},'allch_alpha_dots_PLV.svg']);

    end

end

clear pl_

% saveas(gcf,[Savefolder,band_name{band},'allch_alpha_dots_PLV_allbands.svg']);


% %% Alpha stim PLV from each channel to average of all other channels
% 
% incl_sub = setdiff(1:19,12);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% seeds = 1:127; %[2 34 65 94];
% % seeds = [2];
% 
% % for seed = 1
% 
% for band = 3 %1:5
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
% subplot(1,4,1)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 1:4
%          dat = [dat log10(squeeze(nanmean(PLV_on.band{band}.cond{cond}(:,seeds,incl_sub),2))./squeeze(nanmean(PLV_off.band{band}.cond{cond}(:,seeds,incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
% %     ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% %     colormap(flipud(brewermap(64,'RdBu')))
% %     cb = colorbar;
%     
%     PlotChans = [2];
%     hold on
%     scatter(layout2.pos(PlotChans,1),layout2.pos(PlotChans,2),200,'x','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',2)
%     hold on
%     scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),80,[0.9098 0.4588 0.4275],'filled')
%     load('lapaz.mat');
%     colormap(lapaz)
% %     colorbar
%     
% %     PlotChans = [2];
% %     topoplottest3(stat.stat,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
% %     'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
% %     hold on
% %     PlotChans2 = find(stat.mask == 1);
% %    topoplottest3(stat.stat,EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
% %     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
% %     load('lapaz.mat');
% %     colormap(lapaz)
% %     colorbar
% 
% 
% alpha_PLV_sig_el{band} = find(stat.mask == 1);
%         
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' alpha PLV ANOVA average allch'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_alpha_allch_PLV_ANOVA.svg'])
% end
% 
% 
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% % end
% 
% 
% clear pl_

%% Alpha stim PLV from each channel to average of all other channels - lme

% incl_sub = setdiff(1:19,12);
% 
% seeds = 1:127; %[2 34 65 94];
% % seeds = [2];
% 
% % for seed = 1
% 
% for band1 = 1:4
%     
%     band2 = band1;
%     
% %     table_allch = [];
%     
%     for ch = 1:127
%         
%         perc_change_allcon = [];   
%     
%             for cond = 1:4
% 
%     
%                 PLV_ON_table = squeeze(nanmean(PLV_on.band1{band1}.band2{band2}.cond{cond}(ch,seeds,incl_sub),2));
%                 PLV_OFF_table = squeeze(nanmean(PLV_off.band1{band1}.band2{band2}.cond{cond}(ch,seeds,incl_sub),2));
%         
%                 PLV_ON_phasic_table = squeeze(nanmean(PLV_on_phasic.band1{band1}.band2{band2}.cond{cond}(ch,seeds,incl_sub),2));
%                 PLV_OFF_phasic_table = squeeze(nanmean(PLV_off_phasic.band1{band1}.band2{band2}.cond{cond}(ch,seeds,incl_sub),2));
%                 
%                 PLV_ON_tonic_table = squeeze(nanmean(PLV_on_tonic.band1{band1}.band2{band2}.cond{cond}(ch,seeds,incl_sub),2));
%                 PLV_OFF_tonic_table = squeeze(nanmean(PLV_off_tonic.band1{band1}.band2{band2}.cond{cond}(ch,seeds,incl_sub),2));
%         
%                 perc_change_PLV = log10(PLV_ON_table./PLV_OFF_table); %ifq_alphaband_ON - ifq_alphaband_OFF;
%                 perc_change_PLV_phasic = log10(PLV_ON_phasic_table./PLV_OFF_phasic_table ); %ifq_alphaband_ON_phasic - ifq_alphaband_OFF_phasic;
%                 perc_change_PLV_tonic = log10(PLV_ON_tonic_table./PLV_OFF_tonic_table); %ifq_alphaband_ON_tonic - ifq_alphaband_OFF_tonic;
%        
%                 perc_change_phasictonic = vertcat(perc_change_PLV_phasic,perc_change_PLV_tonic);
% 
%                 perc_change_allcon = vertcat(perc_change_allcon,perc_change_phasictonic);
%         
%                 clear perc_change_PLV_phasic perc_change_PLV_tonic perc_change_phasictonic
%             
%             end
%             
%            
%            	substage_all = [1.*ones(1,length(incl_sub)) 2.*ones(1,length(incl_sub)) 1.*ones(1,length(incl_sub)) 2.*ones(1,length(incl_sub)) 1.*ones(1,length(incl_sub)) 2.*ones(1,length(incl_sub)) 1.*ones(1,length(incl_sub)) 2.*ones(1,length(incl_sub))];
%             cond_all =  [ones(1,length(incl_sub)*2)  2.*ones(1,length(incl_sub)*2) 3.*ones(1,length(incl_sub)*2) 4.*ones(1,length(incl_sub)*2)];
%             sub_all = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
%             electrode = repmat(ch,length(incl_sub)*4*2,1);
%             
%             table_all = table(perc_change_allcon,substage_all',cond_all',electrode,sub_all','VariableNames',{'PLV_change','substage','condition','electrode','sub'});
%             
%             table_all.substage = categorical(table_all.substage);
%             table_all.condition = categorical(table_all.condition);
%             table_all.sub = categorical(table_all.sub);
%             table_all.electrode = categorical(table_all.electrode);
%             
% %             table_allch = vertcat(table_allch,table_all);
% 
%             lme = fitlme(table_all,'PLV_change ~ substage * condition +(1|sub)','FitMethod','REML','DummyVarCoding','effects');
%             
%             stats = anova(lme);
%           
%             Fvalue_substage(band1,ch) = stats.FStat(2); 
%             Fvalue_condition(band1,ch) = stats.FStat(3); 
%             Fvalue_substage_condition(band1,ch) = stats.FStat(4);
%           
%             pvalue_substage(band1,ch) = stats.pValue(2); 
%             pvalue_condition(band1,ch) = stats.pValue(3); 
%             pvalue_substage_condition(band1,ch) = stats.pValue(4);
% 
%     end
%             
% end

%% alphastim lme PLV

for band1 = 1:4
    
    for band2 = band1:4

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)

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
% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

scatter(layout2.pos(2,1),layout2.pos(2,2),400,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),200,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'lme_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end


%% alphastim lme PLV - cluster to electrode

band1 = 3;
band2 = 3;

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);

clear statsresult
PlotChans2 = [];

%  load([statsresult_folder,'statsresult_alphastim_PLV_',band_name{band1},'-',band_name{band2},'.mat']);
   load('/vol/research/nemo/datasets/RSN/data/analysis/connectivity/statsresult_sensor/statsresult_alphastim_PLV_cluster_to_electrode7-12 Hz-7-12 Hz.mat');    
        
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
% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

scatter(layout2.pos(2,1),layout2.pos(2,2),400,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),200,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title(['cluster to electrode ',band_name{band1},'-',band_name{band2}])  

saveas(gcf,[Savefolder,'lme_alphastim_PLV_cluster_to_electrode_',band_name{band1},'-',band_name{band2},'.svg'])

%% PLI Alpha stim on vs off (t-test)

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

for band1 = 1 %1:4

    for band2 = 3 %band1:4
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

for cond = 1:4
    
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
    
    for chan = 1:127 %length(tmp_comps)
        
         x = squeeze(PLI_on.band1{band1}.band2{band2}.cond{cond}(chan,seed,:));
         y = squeeze(PLI_off.band1{band1}.band2{band2}.cond{cond}(chan,seed,:)); 
        
%         x = squeeze(PLI_on.band{band}.cond{cond}(chan,seed,:));
%         y = squeeze(PLI_off.band{band}.cond{cond}(chan,seed,:)); 

%         x = squeeze(PLI_on.band{band}.cond{cond}(tmp_comps(chan),seed,:));  
%         y = squeeze(PLI_off.band{band}.cond{cond}(tmp_comps(chan),seed,:));  

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
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(2,:),'Linewidth',.5);%,.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
            
            else
                
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','b','Linewidth',1);
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(1,:),'Linewidth',.5);%.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
                
            end
            end
            %}
        
        %}
        
    end  
    
    if T>0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(2,:),'filled')
%         scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'r','filled')
    end
        
    if T<0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(1,:),'filled')
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

title([conditions{cond}])
end

% sgtitle([band_name{band},' alpha PLI'])
% cd(figPath)

saveas(gcf,[Savefolder,'ttest_',band_name{band1},'-',band_name{band2},'allch_alpha_dots_PLI.svg']);

    end

end

clear pl_

%% alphastim lme PLI

for band1 = 1:4
    
    for band2 = band1:4

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);


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
% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

scatter(layout2.pos(2,1),layout2.pos(2,2),400,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),200,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

saveas(gcf,[Savefolder,'lme_alphastim_PLI_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end

%% Alpha stim PLV from each channel to average of significant cluster

% incl_sub = setdiff(1:19,12);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% % seeds = [2];
% 
% % for seed = 1
% 
% for band = 3 %1:5
%     
% seeds = alpha_PLV_sig_el{band}; %[2 34 65 94];
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
%     subplot(1,4,1)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 1:4
%          dat = [dat log10(squeeze(nanmean(PLV_on.band{band}.cond{cond}(:,seeds,incl_sub),2))./squeeze(nanmean(PLV_off.band{band}.cond{cond}(:,seeds,incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
%     PlotChans = [2];
%     hold on
%     scatter(layout2.pos(PlotChans,1),layout2.pos(PlotChans,2),200,'x','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',2)
%     hold on
%     scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),80,[0.9098 0.4588 0.4275],'filled')
%     load('lapaz.mat');
%     colormap(lapaz)
%     
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
% % 
% %             scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),50,'w','filled')
%             
%             seed_plot = zeros(127,1);
%             seed_plot(seeds) = 1;
%             seed_plot = logical(seed_plot);
%             scatter(layout2.pos(seed_plot,1),layout2.pos(seed_plot,2),80,[1 0.6 0.6],'filled')
% 
%         
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' alpha PLV ANOVA'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_alpha_cluster','_PLV_ANOVA.svg'])
% end
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% % end
% 
% 
% clear pl_


% %% PLI Alpha stim on vs off (t-test)
% 
% incl_sub = setdiff(1:19,12);
% 
% colors = linspecer(4);
% 
% Y = -.2.*zeros(127,1);
% 
% % comps = NaN(127,127);
% % comps(:,1) = [1:127];
% % 
% % close all
% % 
% % for i = 2:127
% %     
% %     comps(i:127,i) = i:127;
% %     
% % end
% 
% for band = 5 %1:4
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% for cond = 1:4
%     
%     subplot(1,4,cond)
%     hold on
%     %
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),Y,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[-1 1]);
% 
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% %cb = colorbar;
%     %}
%     
%     for seed = 1:127 %[2 34 65 94]
%        
% %     tmp_comps = comps(~isnan(comps(:,seed)),seed);
%     
%     P =0;
%     T = 0;
%     T2 = 0;
%     
%     for chan = 1:127 %length(tmp_comps)
%         
%         x = squeeze(PLI_on.band{band}.cond{cond}(chan,seed,:));
%         y = squeeze(PLI_off.band{band}.cond{cond}(chan,seed,:)); 
% 
% %         x = squeeze(PLV_on.band{band}.cond{cond}(tmp_comps(chan),seed,:));  
% %         y = squeeze(PLV_off.band{band}.cond{cond}(tmp_comps(chan),seed,:));    
%                 
%         [h,p,ci,stats] = ttest(x(incl_sub),y(incl_sub));
%         if ~isnan(stats.tstat)
%             T2 = T2 + stats.tstat;
%         end
% 
%         if p <=.01
%           if ~isnan(stats.tstat)
%             T = T+stats.tstat; 
%           end
%         end
%           
%             %P = P+1;
%             %
%             if p <.01
%             if stats.tstat>0 
%             
%             %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','r','Linewidth',1);
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(2,:),'Linewidth',.5);%,.5*abs(stats.tstat));
% %             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));
% 
%             pl_.Color(4) = .2;
%             
%             else
%                 
%             %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','b','Linewidth',1);
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(1,:),'Linewidth',.5);%.5*abs(stats.tstat));
% %             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));
% 
%             pl_.Color(4) = .2;
%                 
%             end
%             end
%             %}
%         
%         %}
%         
%     end  
%     
%     if T>0 
%         scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(2,:),'filled')
% %         scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'r','filled')
%     end
%         
%     if T<0 
%         scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(1,:),'filled')
% %      scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'b','filled')
% 
%     end
% Ts.band{band}(seed,cond) = T;
%                 Ts2.band{band}(seed,cond) = T2;
% 
%     end
%     
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% title([conditions{cond}])
% end
% 
% sgtitle([band_name{band},' alpha PLI'])
% % cd(figPath)
% % saveas(gcf,[Savefolder,band_name{band},'allch_alpha_dots_PLI.svg']);
% 
% 
% end
% 
% clear pl_


%% Alpha stim PLI from each channel to average of all other channels 

% incl_sub = setdiff(1:19,12);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% seeds = 1:127; %[2 34 65 94];
% % seeds = [2];
% 
% for seed = 1
% 
% for band = 1:5
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
%     subplot(1,4,1)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 1:4
%          dat = [dat log10(squeeze(nanmean(PLI_on.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2))./squeeze(nanmean(PLI_off.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
% 
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
% % 
% %             scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),50,'w','filled')
% 
%     PlotChans = [2];
%     hold on
%     scatter(layout2.pos(PlotChans,1),layout2.pos(PlotChans,2),200,'x','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',2)
%     hold on
%     scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),80,[0.9098 0.4588 0.4275],'filled')
%     load('lapaz.mat');
%     colormap(lapaz)
% 
%         
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' alpha PLI ANOVA average allch'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_alpha_allch_PLI_ANOVA.svg'])
% end
% 
% 
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% end
% 
% 
% clear pl_

% %% Alpha stim PLI from each channel to average of Fz channels
% 
% incl_sub = setdiff(1:19,12);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% seeds = [2 34 65 94];
% % seeds = [2];
% 
% for seed = 1
% 
% for band = 1:5
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
% %     subplot(1,2,cond)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 1:4
%          dat = [dat log10(squeeze(nanmean(PLI_on.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2))./squeeze(nanmean(PLI_off.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
% 
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% cb = colorbar;
% 
%             scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),50,'w','filled')
% 
%         
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' alpha PLI ANOVA'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_alpha_seed_',num2str(seed),'_PLI_ANOVA.svg'])
% end
% 
% 
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% end
% 
% 
% clear pl_

%% PLV Theta stim on vs off (t-test)

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

for band1 = 2 %1:4

    for band2 = 3 %band1:4
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

for cond = 5:8
    
    subplot(1,4,cond-4)
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
    
    for chan = 1:127 %length(tmp_comps)
        
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
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(2,:),'Linewidth',.5);%,.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
            
            else
                
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','b','Linewidth',1);
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(1,:),'Linewidth',.5);%.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
                
            end
            end
            %}
        
        %}
        
    end  
    
    if T>0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(2,:),'filled')
%         scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'r','filled')
    end
        
    if T<0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(1,:),'filled')
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

title([conditions{cond}])
end

% sgtitle([band_name{band},' alpha PLV'])
% cd(figPath)

saveas(gcf,[Savefolder,'ttest_',band_name{band1},'-',band_name{band2},'allch_theta_dots_PLV.svg']);

    end

end

clear pl_

% saveas(gcf,[Savefolder,band_name{band},'allch_alpha_dots_PLV_allbands.svg']);


%% Theta stim PLV from each channel to average of all other channels 

% incl_sub = setdiff(1:19,12);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% seeds = 1:127; %[2 34 65 94];
% % seeds = [2];
% 
% for seed = 1
% 
% for band = 2 %1:5
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
%     subplot(1,4,1)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 5:8
%          dat = [dat log10(squeeze(nanmean(PLV_on.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2))./squeeze(nanmean(PLV_off.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
% 
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
% % 
% %             scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),50,'w','filled')
% 
%     PlotChans = [2];
%     hold on
%     scatter(layout2.pos(PlotChans,1),layout2.pos(PlotChans,2),200,'x','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',2)
%     hold on
%     scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),80,[0.9098 0.4588 0.4275],'filled')
%     load('lapaz.mat');
%     colormap(lapaz)
% 
%         
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' theta PLV ANOVA average allch'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_theta_allch_PLV_ANOVA.svg'])
% end
% 
% 
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% end
% 
% 
% clear pl_

%% thetastim lme PLV

for band1 = 1:4
    
    for band2 = band1:4

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);


 load([statsresult_folder,'statsresult_thetastim_PLV_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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
% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

scatter(layout2.pos(2,1),layout2.pos(2,2),400,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),200,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'lme_thetastim_PLV_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end


%% thetastim lme PLI

for band1 = 1:4
    
    for band2 = band1:4

figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)

% load('EEG_chanlocs.mat');
addpath(genpath('/user/HS301/m17462/matlab/kispi'));
PlotChans = [2];
% PlotChans2 = find(pvalue_condition(band1,:) < 0.05);


 load([statsresult_folder,'statsresult_thetastim_PLI_',band_name{band1},'-',band_name{band2},'.mat']);
        
        
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
% scatter(layout2.pos(2,1),layout2.pos(2,2),200,[0.6350 0.0780 0.1840],'x','LineWidth',5);
% 
% hold on
% scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),50,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

scatter(layout2.pos(2,1),layout2.pos(2,2),400,[0.6350 0.0780 0.1840],'x','LineWidth',2);

hold on
scatter(layout2.pos(PlotChans2,1),layout2.pos(PlotChans2,2),200,[0.9098 0.4588 0.4275],'o','filled','LineWidth',2);

load('lapaz.mat');
colormap(lapaz)
colorbar
set(gca,'FontSize',18);
axis off
axis square
title([band_name{band1},'-',band_name{band2}])  

% saveas(gcf,[Savefolder,'lme_thetastim_PLI_',band_name{band1},'-',band_name{band2},'.svg'])

    end

end

% %% Theta stim PLV from each channel to average of Fz channels
% 
% incl_sub = setdiff(1:19,[3 12]);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% seeds = [2 34 65 94];
% % seeds = [2];
% 
% for seed = 1
% 
% for band = 1:5
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
% %     subplot(1,2,cond)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 5:8
%          dat = [dat log10(squeeze(nanmean(PLV_on.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2))./squeeze(nanmean(PLV_off.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
% 
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu')))
% cb = colorbar;
% 
%             scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),15,'w','filled')
% 
%         
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' theta PLV ANOVA'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_theta_seed_',num2str(seed),'_PLV_ANOVA.svg'])
% end
% 
% 
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% end
% 
% 
% clear pl_


%% PLI Theta stim on vs off (t-test)

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

for band1 = 2 %1:4

    for band2 = 4 %band1:4
    
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
hold on

for cond = 5:8
    
    subplot(1,4,cond-4)
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
    
    for chan = 1:127 %length(tmp_comps)
        
         x = squeeze(PLI_on.band1{band1}.band2{band2}.cond{cond}(chan,seed,:));
         y = squeeze(PLI_off.band1{band1}.band2{band2}.cond{cond}(chan,seed,:)); 
        
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
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(2,:),'Linewidth',.5);%,.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(2,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
            
            else
                
            %pl_ = plot3([elec2.elecpos(seed,1) elec2.elecpos(chan,1)],[elec2.elecpos(seed,2) elec2.elecpos(chan,2)],[elec2.elecpos(seed,3) elec2.elecpos(chan,3)],'Color','b','Linewidth',1);
            pl_ = plot([layout2.pos(seed,1) layout2.pos(chan,1)],[layout2.pos(seed,2) layout2.pos(chan,2)],'Color',colors(1,:),'Linewidth',.5);%.5*abs(stats.tstat));
%             pl_ = plot([layout2.pos(seed,1) layout2.pos(tmp_comps(chan),1)],[layout2.pos(seed,2) layout2.pos(tmp_comps(chan),2)],'Color',colors(1,:),'Linewidth',0.5*abs(stats.tstat));

            pl_.Color(4) = .2;
                
            end
            end
            %}
        
        %}
        
    end  
    
    if T>0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(2,:),'filled')
%         scatter3(elec2.elecpos(seed,1),elec2.elecpos(seed,2),elec2.elecpos(seed,3),.5*abs(T),'r','filled')
    end
        
    if T<0 
        scatter(layout2.pos(seed,1),layout2.pos(seed,2),.5*abs(T),colors(1,:),'filled')
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

title([conditions{cond}])
end

% sgtitle([band_name{band},' alpha PLV'])
% cd(figPath)

saveas(gcf,[Savefolder,'ttest_',band_name{band1},'-',band_name{band2},'allch_theta_dots_PLI.svg']);

    end

end

clear pl_

% saveas(gcf,[Savefolder,band_name{band},'allch_alpha_dots_PLV_allbands.svg']);

%% Theta stim PLI from each channel to average of all other channels 
% 
% incl_sub = setdiff(1:19,12);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% seeds = 1:127; %[2 34 65 94];
% % seeds = [2];
% 
% for seed = 1
% 
% for band = 2 %1:5
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
%     subplot(1,4,1)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 5:8
%          dat = [dat log10(squeeze(nanmean(PLI_on.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2))./squeeze(nanmean(PLI_off.band{band}.cond{cond}(:,seeds(seed,:),incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
% 
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
% 
% %     scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),50,'w','filled')
% 
% 
%     PlotChans = [2];
%     hold on
%     scatter(layout2.pos(PlotChans,1),layout2.pos(PlotChans,2),200,'x','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',2)
%     hold on
%     scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),80,[0.9098 0.4588 0.4275],'filled')
%     load('lapaz.mat');
%     colormap(lapaz)
%         
% axis off
% axis square
% 
%  theta_PLI_sig_el{band} = find(stat.mask == 1);
% 
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' theta PLI ANOVA average allch'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_theta_allch_PLI_ANOVA.svg'])
% end
% 
% 
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% end
% 
% 
% clear pl_


%% Theta stim PLI from each channel to average of significant cluster

% incl_sub = setdiff(1:19,[3 12]);
% 
% cfg = [];
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesFmultivariate';%dependent samples T-test
% cfg.alpha       = 0.05; %significance level
% cfg.numrandomization = 1000;%'all';%number of permutations
% cfg.design = [ones(1,length(incl_sub))  2.*ones(1,length(incl_sub)) 3.*ones(1,length(incl_sub)) 4.*ones(1,length(incl_sub))];%off or on
% cfg.design(2,:) = [1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub) 1:length(incl_sub)];
% cfg.ivar   = 1;%the row on which the independent variable is (in this case, on or off)
% cfg.uvar = 2;
% cfg.parameter = 'parameter';
% cfg.channel = layout2.label;    
% cfg.neighbours = neighbours;    
% cfg.dim = [127 1];
% cfg.correctm = 'cluster';
% %cfg.clusterstatistic = how to combine the single samples that belong to a cluster, 'maxsum', 'maxsize', 'wcm' (default = 'maxsum')
%     
%    
% cfg.clusterthreshold = 'nonparametric_common';%method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
%    %cfg.clusteralpha     = for either parametric or nonparametric thresholding per tail (default = 0.05)
% cfg.clustercritval   = .05;%for parametric thresholding (default is determined by the statfun)
%    %cfg.clustertail      = -1, 1 or 0 (default = 0)
% 
% close all
% 
% % seeds = [2 34 65 94];
% % seeds = [2];
% 
% % for seed = 1
% 
% for band = 2 %1:5
%     
% seeds = theta_PLI_sig_el{band}; %[2 34 65 94];
% 
% figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% hold on
% 
% % for cond = 1:4
%     
%     subplot(1,4,1)
% %     hold on
%   
%     dat = [];
%      
%      for cond = 5:8
%          dat = [dat log10(squeeze(nanmean(PLI_on.band{band}.cond{cond}(:,seeds,incl_sub),2))./squeeze(nanmean(PLI_off.band{band}.cond{cond}(:,seeds,incl_sub),2)))];
%      end
% 
%     stat = ft_statistics_montecarlo(cfg,dat,cfg.design);
% 
%     ft_plot_topo(layout2.pos(:,1),layout2.pos(:,2),stat.stat,'mask',layout2.mask,'outline',layout2.outline, ...
%     'interplim','mask','clim',[0 20],'gridscale',300);  
% 
% 
% % ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% % colormap(flipud(brewermap(64,'RdBu')))
% % cb = colorbar;
% % 
% %             scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),15,'w','filled')
% 
%     PlotChans = [2];
%     hold on
%     scatter(layout2.pos(PlotChans,1),layout2.pos(PlotChans,2),200,'x','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'LineWidth',2)
%     hold on
%     scatter(layout2.pos(stat.mask,1),layout2.pos(stat.mask,2),80,[0.9098 0.4588 0.4275],'filled')
%     load('lapaz.mat');
%     colormap(lapaz)
%             
%             
%             seed_plot = zeros(127,1);
%             seed_plot(seeds) = 1;
%             seed_plot = logical(seed_plot);
%             scatter(layout2.pos(seed_plot,1),layout2.pos(seed_plot,2),80,[1 0.6 0.6],'filled')
% 
%         
% axis off
% axis square
% %scatter(layout2.pos(2,1),layout2.pos(2,2),55,'k','filled')
% %scatter(layout2.pos(47,1),layout2.pos(47,2),55,'k','filled')
% 
%     %scatter3(elec2.elecpos(2,1),elec2.elecpos(2,2),elec2.elecpos(2,3),55,'k','filled')
% 
% % title([conditions{cond}])
% 
% %     end
% sgtitle([band_name{band},' theta PLI ANOVA'])
% % cd(figPath)
% saveas(gcf,[Savefolder,band_name{band},'_theta_cluster_','_PLI_ANOVA.svg'])
% end
% 
% 
% 
% 
% %saveas(gcf,[band_name{band},'dots_PLV.svg'])
% 
% % end
% 
% 
% clear pl_