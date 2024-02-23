clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));
addpath(genpath('/user/HS301/m17462/matlab/ScientificColourMaps7'));
addpath(genpath('/user/HS301/m17462/matlab/kispi'));

load('/vol/research/nemo/datasets/RSN/data/analysis/psd_allsub/psd_allsub_mICA_avref_12-Mar-2023.mat');

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/power_allsub/';

%%

load('/user/HS301/m17462/matlab/Scripts/RSN/analyses/psd/topo/EEG_chanlocs.mat');

f = [0:.1:30];%was [2:.25:30]

band_freq = [1 4;4 7;8 12;13 30;...
             8 9; 9 10; 10 11; 11 12;...
             4 5; 5 6; 6 7; 7 8; 9.5 10.5; 5.5 6.5];

band_name = {'delta (1-4 Hz)';'theta (4-7 Hz)';'alpha (8-12 Hz)';'beta (13-30 Hz)'; ...
              '8-9 Hz'; '9-10 Hz'; '10-11 Hz'; '11-12 Hz';...
              '4-5 Hz'; '5-6 Hz'; '6-7 Hz'; '7-8 Hz'; '9.5-10.5 Hz'; '5.5-6.5 Hz'};
          
incl_sub = setdiff(1:19,12);
          
%% rem topos

for band = 1:4
    
    band_ndx = find(f >= band_freq(band,1) & f <= band_freq(band,2));
    power_band_rem = squeeze(nansum(mfft_rem(incl_sub,:,band_ndx),3)); 
    
%     find(power_band_rem(2,:)> 100)

    figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
    suptitle(band_name{band});

    for s = 1:size(power_band_rem,1)
    
    subplot(3,6,s)
    PlotChans = [2];
    % PlotChans2 = find(mask_bands_alpha(band,:) == 1);
    topoplottest3(power_band_rem(s,:),EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
    'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
    % hold on
    % topoplottest3(stat_bands_alpha(band,:),EEG.chanlocs,'maplimits',[0 30],'conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
    %     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
    load('lapaz.mat');
    colormap(lapaz)
%     colorbar
    set(gca,'FontSize',18);
    
    end
    
    saveas(gcf,[Savefolder,'rem_',band_name{band},'lapaz_colorbar.svg']);

end

%% nrem topos

for band = 1:4
    
    band_ndx = find(f >= band_freq(band,1) & f <= band_freq(band,2));
    power_band_rem = squeeze(nansum(mfft_nrem(incl_sub,:,band_ndx),3)); 
    
%     find(power_band_rem(2,:)> 100)

    figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
    suptitle(band_name{band});

    for s = 1:size(power_band_rem,1)
    
    subplot(3,6,s)
    PlotChans = [2];
    % PlotChans2 = find(mask_bands_alpha(band,:) == 1);
    topoplottest3(power_band_rem(s,:),EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
    'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
    % hold on
    % topoplottest3(stat_bands_alpha(band,:),EEG.chanlocs,'maplimits',[0 30],'conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
    %     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
    load('lapaz.mat');
    colormap(lapaz)
%     colorbar
    set(gca,'FontSize',18);
    
    end
    
    saveas(gcf,[Savefolder,'nrem_',band_name{band},'lapaz_colorbar.svg']);

end


%% wake eve topos

for band = 1:4
    
    band_ndx = find(f >= band_freq(band,1) & f <= band_freq(band,2));
    power_band_rem = squeeze(nansum(mfft_wake_e(incl_sub,:,band_ndx),3)); 
    
%     find(power_band_rem(2,:)> 100)

    figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
    suptitle(band_name{band});

    for s = 1:size(power_band_rem,1)
    
    subplot(3,6,s)
    PlotChans = [2];
    % PlotChans2 = find(mask_bands_alpha(band,:) == 1);
    topoplottest3(power_band_rem(s,:),EEG.chanlocs,'maplimits','minmax','conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
    'emarker2',{PlotChans,'x',[0.6350 0.0780 0.1840],10,3});
    % hold on
    % topoplottest3(stat_bands_alpha(band,:),EEG.chanlocs,'maplimits',[0 30],'conv','on','intrad',0.5,'gridscale',200,'electrodes','off','whitebk','on',...
    %     'emarker2',{PlotChans2,'o',[0.9098 0.4588 0.4275],5,2});
    load('lapaz.mat');
    colormap(lapaz)
%     colorbar
    set(gca,'FontSize',18);
    
    end
    
    saveas(gcf,[Savefolder,'wake_e_',band_name{band},'lapaz_colorbar.svg']);

end