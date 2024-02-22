clear all;
close all;

addpath(genpath('\\surrey.ac.uk\personal\hs301\m17462\matlab\eeglab')); % eeglab toolbox, see README on where to find this
addpath(genpath('\\surrey.ac.uk\personal\hs301\m17462\matlab\Henry\useful_functions')); % contains linspecer function, circular statistics toolbox functions, echt function, shadedErrorBar function, see README on where to find this

%% Load data

Folderpath = 'S:\datasets\RSN\data\hdEEG\';
sub_Folderpath = dir([Folderpath,'RSN*']);

for s = 1:length(sub_Folderpath)
    
mICA_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*sleep*_mICA_avref.set']);
nm_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*sleep*_nm_good.mat']);
goodREM_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_czref_goodREM.mat']);

Savefolder = [Folderpath,sub_Folderpath(s).name];

EEG = pop_loadset('filename',[mICA_file(1).name],'filepath',[Folderpath,sub_Folderpath(s).name]);
load([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file(1).name]);
load([Folderpath,sub_Folderpath(s).name,filesep,goodREM_file(1).name]);

data = double(EEG.data);
clear EEG.data

%% Alpha - Phase analysis

filt_lf = 7.5;
filt_hf = 12.5;
fs = 500;

colors = linspecer(4);

for con = 1:8
    
    for ch = 1:size(data,1)
 
        % based on Fz of hdEEG
        for trig = 1:length(nm.ON_trigs_StimTrak_good{con})
       
            hilb_fz = echt(data(ch,nm.ON_trigs_StimTrak_good{con}(trig)-fs:nm.ON_trigs_StimTrak_good{con}(trig)), filt_lf, filt_hf, fs);    
            stim_phase_good(ch,trig) = angle(hilb_fz(end)); 
       
            clear hilb_fz
        end  
 
        if length(nm.ON_trigs_StimTrak_good{con})==0
           stim_phase_good = [];
        end
 
    end
    
    nm.stimphase_good_alphafilt{con} = stim_phase_good;
  
    clear stim_phase_echt stim_phase_good

end


%% Alpha - Make polar plot 

% Make polar plot based on echt data
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

ch = 2;

for con = 1:4

    stim_phase_good = nm.stimphase_good{con};
    
    if ~isempty(stim_phase_good)
    
    polarhistogram(stim_phase_good(ch,:),100,'FaceColor',colors(con,:),'LineWidth',.75)  % Polar Histogram   
    hold on
%     title({['RSN ',mICA_file(1).name(5:7)] ['Alpha (Fz data)']});
    
    end

    clear stim_phase_good
    
end

if ~exist([Savefolder,filesep,'QC_Plots'])
    mkdir([Savefolder,filesep,'QC_Plots']);   
end
% print([Savefolder,filesep,'QC_Plots/Alpha_phase_plot_fz_good'],'-dpng');

set(gca,'Fontsize',35);

saveas(gcf,[Savefolder,filesep,'QC_Plots/Alpha_phase_plot_fz_good.svg']);

%% Theta - Phase analysis

filt_lf = 4.5;
filt_hf = 7.5;

colors = linspecer(4);

for con = 1:8
    
    for ch = 1:size(data,1)

        % based on Fz of hdEEG
        for trig = 1:length(nm.ON_trigs_StimTrak_good{con})
       
            hilb_fz = echt(data(ch,nm.ON_trigs_StimTrak_good{con}(trig)-fs:nm.ON_trigs_StimTrak_good{con}(trig)), filt_lf, filt_hf, fs);    
            stim_phase_good(ch,trig) = angle(hilb_fz(end)); 
       
            clear hilb_fz
        end  
    
        if length(nm.ON_trigs_StimTrak_good{con})==0
        stim_phase_good = [];
        end

    end
    
    nm.stimphase_good_thetafilt{con} = stim_phase_good;
  
    clear stim_phase_echt stim_phase_good
end


%% Theta - Make polar plot 

% Make polar plot based on hdEEG Fz data
figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

ch = 2;

for con = 5:8

    stim_phase_good = nm.stimphase_good{con};
    
    if ~isempty(stim_phase_good)
    
    polarhistogram(stim_phase_good(ch,:),100,'FaceColor',colors(con-4,:),'LineWidth',.75)  % Polar Histogram   
    hold on
%     title({['RSN ',mICA_file(1).name(5:7)] ['Theta (Fz data)']});
    
    end

    clear stim_phase_good
    
end

if ~exist([Savefolder,filesep,'QC_Plots'])
    mkdir([Savefolder,filesep,'QC_Plots']);   
end
% print([Savefolder,filesep,'QC_Plots/Theta_phase_plot_fz_good'],'-dpng');

set(gca,'Fontsize',35);

saveas(gcf,[Savefolder,filesep,'QC_Plots/Theta_phase_plot_fz_good.svg']);


%%

save([Savefolder,filesep,mICA_file(1).name(1:end-4),'_nm_good_phase_allch.mat'],'nm','-v7.3');

clear nm

end





