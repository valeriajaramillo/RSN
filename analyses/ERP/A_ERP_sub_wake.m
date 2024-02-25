clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));

addpath /user/HS301/m17462/matlab/Henry/useful_functions
addpath /user/HS301/m17462/matlab/colorGradient

erpfolder = '/vol/research/nemo/datasets/RSN/data/analysis/erp_allsub/';

%% Load data

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/';
sub_Folderpath = dir([Folderpath,'RSN*']);

dt = 1;
fs = 500;

%%
for s = 1:length(sub_Folderpath)
    
display(s); 

mICA_file_wake_e = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_wake_e_fil_czref_mICA_avref.set']);
nm_good_file_wake_e = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_wake_e_nm_good.mat']);
goodREM_file_wake_e = dir([Folderpath,sub_Folderpath(s).name,filesep,'*wake_e_fil_czref_goodwake.mat']);

mICA_file_wake_m = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_wake_m_fil_czref_mICA_avref.set']);
nm_good_file_wake_m = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_wake_m_nm_good.mat']);
goodREM_file_wake_m = dir([Folderpath,sub_Folderpath(s).name,filesep,'*wake_m_fil_czref_goodwake.mat']);

Savefolder = [Folderpath,sub_Folderpath(s).name];

ch = [2 34 65 94]; % Fz, AFz, AFF1h, AFF2h
ERP.channels = ch;
     
     %% Wake evening
     
    clear EEG data data_stim trigs vol* nm rem_goodsamp2 wake_goodsamp2 wake_good_trial trial_data*
     
    EEG = pop_loadset('filename',[mICA_file_wake_e(1).name],'filepath',[Folderpath,sub_Folderpath(s).name],'loadmode', ch);
    load([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file_wake_e(1).name]);
    load([Folderpath,sub_Folderpath(s).name,filesep,goodREM_file_wake_e(1).name]);

    data = double(EEG.data);
    clear EEG.data

    data_stim = nanmean(data,1);

  %%
    for state = 1:length(nm.trigs_StimTrak_ERP)
        
        trigs = nm.trigs_StimTrak_ERP{state};
        ERP.trigs{state} = trigs;
        
   %% phase of stimuli

    % alpha phase analysis
    filt_lf_alpha = 7.5;
    filt_hf_alpha = 12.5;

        for trig = 1:length(trigs)
      
            if ismember(trigs(trig),wake_goodsamp2)
       
                hilb_fz = echt(data_stim(trigs(trig)-fs:trigs(trig)), filt_lf_alpha, filt_hf_alpha, fs);    
                stim_phase_good_alpha(trig) = angle(hilb_fz(end)); 
       
                clear hilb_fz
      
            else
          
                stim_phase_good_alpha(trig) = NaN; 
          
            end
      
        end 
 
                % theta phase analysis 
     filt_lf_theta = 4.5;
     filt_hf_theta = 7.5;

        for trig = 1:length(trigs)
      
            if ismember(trigs(trig),wake_goodsamp2)
       
                hilb_fz = echt(data_stim(trigs(trig)-fs:trigs(trig)), filt_lf_theta, filt_hf_theta, fs);    
                stim_phase_good_theta(trig) = angle(hilb_fz(end)); 
       
                clear hilb_fz
      
            else
          
                stim_phase_good_theta(trig) = NaN; 
          
            end
      
        end 
  
  
  
%         polarhistogram(stim_phase_good_alpha,50,'LineWidth',.75)  % Polar Histogram  

        stim_phase_good_deg_alpha = rad2deg(stim_phase_good_alpha);
        stim_phase_good_deg_alpha(find(stim_phase_good_deg_alpha < 0)) = stim_phase_good_deg_alpha(find(stim_phase_good_deg_alpha < 0))+360;

        stim_phase_good_deg_theta = rad2deg(stim_phase_good_theta);
        stim_phase_good_deg_theta(find(stim_phase_good_deg_theta < 0)) = stim_phase_good_deg_theta(find(stim_phase_good_deg_theta < 0))+360;

        % histogram(stim_phase_good_deg)  

        % group them into 4 phase bins
        ERP.trigs_con{state,1} = find(stim_phase_good_deg_alpha >= 300 | stim_phase_good_deg_alpha <= 30); 
        ERP.trigs_con{state,2} = find(stim_phase_good_deg_alpha >= 30 & stim_phase_good_deg_alpha <= 120);
        ERP.trigs_con{state,3} = find(stim_phase_good_deg_alpha >= 120 & stim_phase_good_deg_alpha <= 210);
        ERP.trigs_con{state,4} = find(stim_phase_good_deg_alpha >= 210 & stim_phase_good_deg_alpha <= 300);

        ERP.trigs_con{state,5} = find(stim_phase_good_deg_theta >= 300 | stim_phase_good_deg_theta <= 30); 
        ERP.trigs_con{state,6} = find(stim_phase_good_deg_theta >= 30 & stim_phase_good_deg_theta <= 120);
        ERP.trigs_con{state,7} = find(stim_phase_good_deg_theta >= 120 & stim_phase_good_deg_theta <= 210);
        ERP.trigs_con{state,8} = find(stim_phase_good_deg_theta >= 210 & stim_phase_good_deg_theta <= 300);
        
%%
        
        if ~isempty(trigs)
        
        trial_data= [];
          
        for t = 1:length(trigs)
            trial_data(t,:) = data_stim(trigs(t)-dt*fs:trigs(t)+dt*fs-1); % 501. sample is trig sample

            for time = 1:2*fs %-dt*fs:dt*fs-1
                x_alpha = echt(data_stim(trigs(t)-dt*fs+time-1-fs:trigs(t)-dt*fs+time-1), filt_lf_alpha, filt_hf_alpha, fs); 
                trial_data_alphafilt_phase(t,time) = angle(x_alpha(end)); % 501. sample is trig sample
                trial_data_alphafilt_real(t,time) = real(x_alpha(end)); % get real component of hilbert-transformed signal
                x_theta = echt(data_stim(trigs(t)-dt*fs+time-1-fs:trigs(t)-dt*fs+time-1), filt_lf_theta, filt_hf_theta, fs);
                trial_data_thetafilt_phase(t,time) = angle(x_theta(end));
                trial_data_thetafilt_real(t,time) = real(x_theta(end)); % get real component of hilbert-transformed signal
                
                clear x_alpha x_theta               
            end
            
            clear trig_ndx
        end
        
        ERP.trial_data{state} = trial_data;
        ERP.trial_data_alphafilt_phase{state} = trial_data_alphafilt_phase;
        ERP.trial_data_alphafilt_real{state} = trial_data_alphafilt_real;
        ERP.trial_data_thetafilt_phase{state} = trial_data_thetafilt_phase;
        ERP.trial_data_thetafilt_real{state} = trial_data_thetafilt_real;
     
        else
            
        ERP.trial_data{state} = NaN(1,2*fs);
        ERP.trial_data_alphafilt_phase{state} = NaN(1,2*fs);
        ERP.trial_data_alphafilt_real{state} = NaN(1,2*fs);
        ERP.trial_data_thetafilt_phase{state} = NaN(1,2*fs);
        ERP.trial_data_thetafilt_real{state} = NaN(1,2*fs);
        
        end
        
    end
    
    save([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file_wake_e(1).name(1:16),'_erp.mat'],'ERP');



     %% Wake morning
     
    clear EEG data data_stim trigs vol* nm rem_goodsamp2 wake_goodsamp2 wake_good_trial trial_data*
     
    EEG = pop_loadset('filename',[mICA_file_wake_m(1).name],'filepath',[Folderpath,sub_Folderpath(s).name],'loadmode', ch);
    load([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file_wake_m(1).name]);
    load([Folderpath,sub_Folderpath(s).name,filesep,goodREM_file_wake_m(1).name]);

    data = double(EEG.data);
    clear EEG.data

    data_stim = nanmean(data,1);

  %%
    for state = 1:length(nm.trigs_StimTrak_ERP)
        
        trigs = nm.trigs_StimTrak_ERP{state};
        ERP.trigs{state} = trigs;
        
   %% phase of stimuli

    % alpha phase analysis
    filt_lf_alpha = 7.5;
    filt_hf_alpha = 12.5;

        for trig = 1:length(trigs)
      
            if ismember(trigs(trig),wake_goodsamp2)
       
                hilb_fz = echt(data_stim(trigs(trig)-fs:trigs(trig)), filt_lf_alpha, filt_hf_alpha, fs);    
                stim_phase_good_alpha(trig) = angle(hilb_fz(end)); 
       
                clear hilb_fz
      
            else
          
                stim_phase_good_alpha(trig) = NaN; 
          
            end
      
        end 
 
                % theta phase analysis 
     filt_lf_theta = 4.5;
     filt_hf_theta = 7.5;

        for trig = 1:length(trigs)
      
            if ismember(trigs(trig),wake_goodsamp2)
       
                hilb_fz = echt(data_stim(trigs(trig)-fs:trigs(trig)), filt_lf_theta, filt_hf_theta, fs);    
                stim_phase_good_theta(trig) = angle(hilb_fz(end)); 
       
                clear hilb_fz
      
            else
          
                stim_phase_good_theta(trig) = NaN; 
          
            end
      
        end 
  
  
  
%         polarhistogram(stim_phase_good_alpha,50,'LineWidth',.75)  % Polar Histogram  

        stim_phase_good_deg_alpha = rad2deg(stim_phase_good_alpha);
        stim_phase_good_deg_alpha(find(stim_phase_good_deg_alpha < 0)) = stim_phase_good_deg_alpha(find(stim_phase_good_deg_alpha < 0))+360;

        stim_phase_good_deg_theta = rad2deg(stim_phase_good_theta);
        stim_phase_good_deg_theta(find(stim_phase_good_deg_theta < 0)) = stim_phase_good_deg_theta(find(stim_phase_good_deg_theta < 0))+360;

        % histogram(stim_phase_good_deg)  

        % group them into 4 phase bins
        ERP.trigs_con{state,1} = find(stim_phase_good_deg_alpha >= 300 | stim_phase_good_deg_alpha <= 30); 
        ERP.trigs_con{state,2} = find(stim_phase_good_deg_alpha >= 30 & stim_phase_good_deg_alpha <= 120);
        ERP.trigs_con{state,3} = find(stim_phase_good_deg_alpha >= 120 & stim_phase_good_deg_alpha <= 210);
        ERP.trigs_con{state,4} = find(stim_phase_good_deg_alpha >= 210 & stim_phase_good_deg_alpha <= 300);

        ERP.trigs_con{state,5} = find(stim_phase_good_deg_theta >= 300 | stim_phase_good_deg_theta <= 30); 
        ERP.trigs_con{state,6} = find(stim_phase_good_deg_theta >= 30 & stim_phase_good_deg_theta <= 120);
        ERP.trigs_con{state,7} = find(stim_phase_good_deg_theta >= 120 & stim_phase_good_deg_theta <= 210);
        ERP.trigs_con{state,8} = find(stim_phase_good_deg_theta >= 210 & stim_phase_good_deg_theta <= 300);
        
%%
        
        if ~isempty(trigs)
        
        trial_data= [];
          
        for t = 1:length(trigs)
            trial_data(t,:) = data_stim(trigs(t)-dt*fs:trigs(t)+dt*fs-1); % 501. sample is trig sample

            for time = 1:2*fs %-dt*fs:dt*fs-1
                x_alpha = echt(data_stim(trigs(t)-dt*fs+time-1-fs:trigs(t)-dt*fs+time-1), filt_lf_alpha, filt_hf_alpha, fs); 
                trial_data_alphafilt_phase(t,time) = angle(x_alpha(end)); % 501. sample is trig sample
                trial_data_alphafilt_real(t,time) = real(x_alpha(end)); % get real component of hilbert-transformed signal
                x_theta = echt(data_stim(trigs(t)-dt*fs+time-1-fs:trigs(t)-dt*fs+time-1), filt_lf_theta, filt_hf_theta, fs);
                trial_data_thetafilt_phase(t,time) = angle(x_theta(end));
                trial_data_thetafilt_real(t,time) = real(x_theta(end)); % get real component of hilbert-transformed signal
                
                clear x_alpha x_theta               
            end
            
            clear trig_ndx
        end
        
        ERP.trial_data{state} = trial_data;
        ERP.trial_data_alphafilt_phase{state} = trial_data_alphafilt_phase;
        ERP.trial_data_alphafilt_real{state} = trial_data_alphafilt_real;
        ERP.trial_data_thetafilt_phase{state} = trial_data_thetafilt_phase;
        ERP.trial_data_thetafilt_real{state} = trial_data_thetafilt_real;
     
        else
            
        ERP.trial_data{state} = NaN(1,2*fs);
        ERP.trial_data_alphafilt_phase{state} = NaN(1,2*fs);
        ERP.trial_data_alphafilt_real{state} = NaN(1,2*fs);
        ERP.trial_data_thetafilt_phase{state} = NaN(1,2*fs);
        ERP.trial_data_thetafilt_real{state} = NaN(1,2*fs);
        
        end
        
    end
    
    save([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file_wake_m(1).name(1:16),'_erp.mat'],'ERP');

    
    end

   