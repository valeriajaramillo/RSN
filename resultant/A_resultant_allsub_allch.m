clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/';
sub_Folderpath = dir([Folderpath,'RSN*']);

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/phase_allsub/';

%% Load nm data

% r = NaN(length(sub_Folderpath),128,8);
% m = NaN(length(sub_Folderpath),128,8);
% std = NaN(length(sub_Folderpath),128,8);

r_alphafilt = NaN(length(sub_Folderpath),128,8);
m_alphafilt = NaN(length(sub_Folderpath),128,8);
std_alphafilt = NaN(length(sub_Folderpath),128,8);

r_thetafilt = NaN(length(sub_Folderpath),128,8);
m_thetafilt  = NaN(length(sub_Folderpath),128,8);
std_thetafilt = NaN(length(sub_Folderpath),128,8);

for s = 1:length(sub_Folderpath)
    
    display(sub_Folderpath(s).name);

%     nm_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_avref_nm_good_psd_allch.mat']);
    nm_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_avref_nm_good_phase_allch.mat']);

    load([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file(1).name]);

    for con = 1:length(nm.psd_ON)
        
%         stimphase_good_con = nm.stimphase_good{con}; % psd: bins x ep x trials x ch
        stimphase_good_alphafilt_con = nm.stimphase_good_alphafilt{con}; % psd: bins x ep x trials x ch
        stimphase_good_thetafilt_con = nm.stimphase_good_thetafilt{con}; % psd: bins x ep x trials x ch

        
%         for ch = 1:size(stimphase_good_con,1)
%             
%             if ~ isempty(stimphase_good_con)
%         
%                 stimphase_good_ch = stimphase_good_con(ch,:);
%                 r(s,ch,con) = circ_r(stimphase_good_ch');
%                 m(s,ch,con) = circ_mean(stimphase_good_ch');
%                 std(s,ch,con) = circ_std(stimphase_good_ch');
%        
%                 clear stimphase_good_ch
%             
%             end
%         
%         end

%         clear stimphase_good_con


        for ch = 1:size(stimphase_good_alphafilt_con,1)
            
            if ~ isempty(stimphase_good_alphafilt_con)
        
                stimphase_good_alphafilt_ch = stimphase_good_alphafilt_con(ch,:);
                r_alphafilt(s,ch,con) = circ_r(stimphase_good_alphafilt_ch');
                m_alphafilt(s,ch,con) = circ_mean(stimphase_good_alphafilt_ch');
                std_alphafilt(s,ch,con) = circ_std(stimphase_good_alphafilt_ch');
       
                clear stimphase_good_alphafilt_ch
            
            end
        
        end
        
        clear stimphase_good_alphafilt_con
        
        
        for ch = 1:size(stimphase_good_thetafilt_con,1)
            
            if ~ isempty(stimphase_good_thetafilt_con)
        
                stimphase_good_thetafilt_ch = stimphase_good_thetafilt_con(ch,:);
                r_thetafilt(s,ch,con) = circ_r(stimphase_good_thetafilt_ch');
                m_thetafilt(s,ch,con) = circ_mean(stimphase_good_thetafilt_ch');
                std_thetafilt(s,ch,con) = circ_std(stimphase_good_thetafilt_ch');
       
                clear stimphase_good_thetafilt_ch
            
            end
        
        end
        
        clear stimphase_good_thetafilt_con
        
                
    end
    
    
end

%%

if ~exist(Savefolder)
   mkdir(Savefolder) 
end

% save([Savefolder,'phase_allsub_mICA_avref_',date,'.mat'],'r','m','std');
save([Savefolder,'phase_allsub_mICA_avref_alphathetafilt',date,'.mat'],'r_alphafilt','m_alphafilt','std_alphafilt','r_thetafilt','m_thetafilt','std_thetafilt');

%%
ch = 2;

r_alpha = squeeze(nanmean(r(:,ch,1:4),3));
r_theta = squeeze(nanmean(r(:,ch,5:8),3));

m_r_alpha = nanmean(r_alpha);
m_r_theta = nanmean(r_theta);

sd_r_alpha = nanstd(r_alpha);
sd_r_theta = nanstd(r_theta);


