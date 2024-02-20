clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
% addpath(genpath('/user/HS301/m17462/matlab/sprint'));
% addpath(genpath('/user/HS301/m17462/matlab/brainstorm3-master'));

addpath(genpath('/user/HS301/m17462/matlab/eBOSC'));

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/';
sub_Folderpath = dir([Folderpath,'RSN*']);

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/oscillation_detection/';
% Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/oscillation_detection/30_45_Hz/';

%% Average across on and off blocks and calculate change
lower_freq_alpha = 7.5;
higher_freq_alpha = 12.5;
lower_freq_theta = 4.5;
higher_freq_theta = 7.5;
%%

for s = [4 17 18 19]
    
mICA_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_sleep*_fil_czref_mICA_avref.set']);
EEG = pop_loadset('filename',[mICA_file(1).name],'filepath',[Folderpath,sub_Folderpath(s).name]);
    
goodREM_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_sleep*_fil_czref_goodREM.mat']);
load([Folderpath,sub_Folderpath(s).name,filesep,goodREM_file(1).name]);
    
nm_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_nm_good_psd_allch.mat']);
load([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file(1).name]); 

%% eBOSC settings
cfg = [];
cfg.eBOSC.channel = []; %[];%all channels
cfg.eBOSC.trial = 1;
cfg.eBOSC.trial_background = 1;

%cfg.eBOSC.F           = 1:.5:44;    % frequency sampling
cfg.eBOSC.F = 1:.25:30; %1:.25:30;%2.^[0:.0125:5];%1:.25:18;%
% f = cfg.eBOSC.F;
cfg.eBOSC.wavenumber  = 6;% wavelet family parameter (time-frequency tradeoff)
cfg.eBOSC.fsample     = fs;     % current sampling frequency of EEG data

cfg.eBOSC.pad.tfr_s = 1;            % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.eBOSC.pad.detection_s = .5;     % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.eBOSC.pad.background_s = 1;     % padding of segments for BG (only avoiding edge artifacts)
cfg.eBOSC.threshold.excludePeak = [];                                   % lower and upper bound of frequencies to be excluded during background fit (Hz) (previously: LowFreqExcludeBG HighFreqExcludeBG)
cfg.eBOSC.threshold.duration	= repmat(3, 1, numel(cfg.eBOSC.F));         % vector of duration thresholds at each frequency (previously: ncyc)
cfg.eBOSC.threshold.percentile  = .95;                                      % percentile of background fit for power threshold

cfg.eBOSC.postproc.use      = 'no';         % Post-processing of rhythmic eBOSC.episodes, i.e., wavelet 'deconvolution' (default = 'no')
cfg.eBOSC.postproc.method   = 'MaxBias';	% Deconvolution method (default = 'MaxBias', FWHM: 'FWHM')
cfg.eBOSC.postproc.edgeOnly = 'yes';        % Deconvolution only at on- and offsets of eBOSC.episodes? (default = 'yes')
cfg.eBOSC.postproc.effSignal= 'PT'; 

%%


% rem_epochs = find(hypno2 == 'R');
data = double(EEG.data);
waves_allepch = [];

for ep = 1:nepochs
    
    for ch = 2 %1:128

%     samp_ep = rem_epochs(ep)*epochl*fs+1:(rem_epochs(ep)+1)*epochl*fs;
    samp_ep = (ep-1)*epochl*fs+1:ep*epochl*fs;
    
    signal.label{1} =  'Fz';
    signal.trial{1} = data(ch,samp_ep);
    signal.time{1} = (1:length(signal.trial{1}))./fs;


    [eBOSC, ~] = eBOSC_wrapper(cfg, signal);
    
    tmp_abu = squeeze(eBOSC.abundance_ep);
%     if size(tmp_abu,1)~= 1
%        tmp_abu = tmp_abu';
%     end
    f = cfg.eBOSC.F;
    abundance_alpha(ep) = sum(tmp_abu(f>= lower_freq_alpha & f<= higher_freq_alpha)); %abundance
    abundance_theta(ep) = sum(tmp_abu(f>= lower_freq_theta & f<= higher_freq_theta)); %abundance
    
    criteria_alpha = eBOSC.episodes.FrequencyMean >= lower_freq_alpha & eBOSC.episodes.FrequencyMean <= higher_freq_alpha; %& eBOSC.episodes.DurationC >= 3; %eBOSC.episodes.FrequencyMean >= bandf(band,1) & eBOSC.episodes.FrequencyMean <= bandf(band,2)
    density_alpha(ep) = sum(criteria_alpha)/epochl*60; %density (number per minute)
    
    criteria_theta = eBOSC.episodes.FrequencyMean >= lower_freq_theta & eBOSC.episodes.FrequencyMean <= higher_freq_theta; %& eBOSC.episodes.DurationC >= 3; %eBOSC.episodes.FrequencyMean >= bandf(band,1) & eBOSC.episodes.FrequencyMean <= bandf(band,2)
    density_theta(ep) = sum(criteria_theta)/epochl*60; %density (number per minute)

   channel  = repmat(ch,length(eBOSC.episodes.Channel),1);
%    startsamp = (rem_epochs(ep)-1)*fs*epochl + eBOSC.episodes.Onset*fs;
%    endsamp = (rem_epochs(ep)-1)*fs*epochl + eBOSC.episodes.Offset*fs;
   startsamp = (ep-1)*fs*epochl + eBOSC.episodes.Onset*fs;  
   endsamp = (ep-1)*fs*epochl + eBOSC.episodes.Offset*fs;
   duration = (endsamp-startsamp)/fs;
   cycles = eBOSC.episodes.DurationC;
   power = eBOSC.episodes.PowerMean;
   SNR = eBOSC.episodes.SNRMean;
   frequency = eBOSC.episodes.FrequencyMean;
   nep = repmat(ep,length(startsamp),1);
   stage = repelem(hypno2(ep),length(startsamp),1);
   offset = repmat(eBOSC.static.pv(:,2),length(startsamp),1);
   exponent = repmat(eBOSC.static.pv(:,1),length(startsamp),1);
   
%    alpha_ndx = find(frequency > 0 & frequency < 30);
%    cycles_ndx = find(cycles > 3);
%    alpha_cycles_ndx = intersect(alpha_ndx,cycles_ndx);
% %    duration_ndx = find(duration > 0.5);
% %    alpha_cycles_duration_ndx = intersect(alpha_cycles_ndx,duration_ndx);
%    SNR_ndx = find(SNR > 1);
%    alpha_cycles_SNR_ndx = intersect(alpha_cycles_ndx,SNR_ndx);
   
%    subplot(2,1,1)
%    plot(signal.trial{1});
%    hold on
%    for a = 1:length(alpha_cycles_SNR_ndx)
%    plot_samps = round(eBOSC.episodes.Onset(alpha_cycles_SNR_ndx(a))*fs:eBOSC.episodes.Offset(alpha_cycles_SNR_ndx(a))*fs);
%    plot(plot_samps,signal.trial{1}(plot_samps),'r')
%    hold on
%    end
%    subplot(2,1,2)
%    plot(samp_ep,data(ch,samp_ep))
%    hold on
%    for a = 1:length(alpha_cycles_SNR_ndx)   
%    plot_samps_orig = round(startsamp(alpha_cycles_SNR_ndx(a)):endsamp(alpha_cycles_SNR_ndx(a)));
%    plot(plot_samps_orig,data(ch,plot_samps_orig),'r')
%    hold on
%    end
%    
%    histogram(frequency(alpha_cycles_SNR_ndx),10)
                
   waves_ep = table(channel,startsamp,endsamp,duration,cycles,power,SNR,frequency,nep,stage,offset,exponent);
   waves_allepch= vertcat(waves_allepch,waves_ep);
   clear waves_ep channel startsamp endsamp duration cycles power SNR frequency nep stage offset exponent

    end

end

save([Savefolder,sub_Folderpath(s).name,'_eBOSC_waves.mat'],'waves_allepch','abundance_alpha','abundance_theta','density_alpha','density_theta','-v7.3')

clear abundance_alpha abundance_theta density_alpha density_theta

end

%%

% waves_ep_3cyc = waves_allepch(waves_allepch.cycles > 3,:);
% waves_ep_3cyc_rem = waves_ep_3cyc(waves_ep_3cyc.stage == 'R',:);
% waves_ep_3cyc_wake = waves_ep_3cyc(waves_ep_3cyc.stage == 'W',:);
% waves_ep_3cyc_n1 = waves_ep_3cyc(waves_ep_3cyc.stage == '1',:);
% waves_ep_3cyc_n2_n3 = waves_ep_3cyc(waves_ep_3cyc.stage == '2'|waves_ep_3cyc.stage == '3',:);

%%
% % subplot(1,4,1)
% histogram(waves_ep_3cyc_wake.frequency)
% hold on
% % title('Wake');
% 
% % subplot(1,4,4)
% histogram(waves_ep_3cyc_rem.frequency)
% % title('REM')
% hold on
% 
% % subplot(1,4,2)
% histogram(waves_ep_3cyc_n1.frequency)
% % title('N1')
% hold on
% 
% % subplot(1,4,3)
% histogram(waves_ep_3cyc_n2_n3.frequency)
% % title('N2 and N3')
% % hold on
% 
% 
% 
% 
% 
% 
% 
% 
