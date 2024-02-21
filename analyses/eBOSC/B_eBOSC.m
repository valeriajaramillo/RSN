clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/eBOSC'));

addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/';
sub_Folderpath = dir([Folderpath,'RSN*']);

waves_folder = '/vol/research/nemo/datasets/RSN/data/analysis/oscillation_detection/';

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/oscillation_detection/';

load('/vol/research/nemo/datasets/RSN/data/analysis/ISI_allsub/ISI_echt_psd_allsub_14-Mar-2023.mat')


%% Average across on and off blocks and calculate change
lower_freq_alpha = 7.5;
higher_freq_alpha = 12.5;
lower_freq_theta = 4.5;
higher_freq_theta = 7.5;

%%

for s = 1:length(sub_Folderpath)   

mICA_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_sleep*_fil_czref_mICA_avref.set']);
EEG = pop_loadset('filename',[mICA_file(1).name],'filepath',[Folderpath,sub_Folderpath(s).name]);

auxch_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_sleep*_auxch_all.set']);
EEG_aux = pop_loadset('filename',[auxch_file(1).name],'filepath',[Folderpath,sub_Folderpath(s).name]);
    
goodREM_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_sleep*_fil_czref_goodREM.mat']);
load([Folderpath,sub_Folderpath(s).name,filesep,goodREM_file(1).name]);
%     
% nm_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_nm_good_psd_allch.mat']);
% load([Folderpath,sub_Folderpath(s).name,filesep,nm_good_file(1).name]); 

waves_file = dir([waves_folder,sub_Folderpath(s).name,'*_eBOSC_waves.mat']);
load([waves_folder,waves_file(1).name])

%%

 good_ndx = [];
   for w = 1:length(waves_allepch.startsamp)
       samps = waves_allepch.startsamp(w):waves_allepch.endsamp(w);
       samps_good = intersect(samps,rem_goodsamp2);
       if isequal(length(samps),length(samps_good))
           good_ndx = [good_ndx w]; 
       end 
   end
   
  waves_allepch_good = waves_allepch(good_ndx,:);
  
   data = double(EEG.data);
   lEOG = EEG_aux.data(3,:);
   rEOG = EEG_aux.data(4,:);
   EMG = EEG_aux.data(2,:);
   ECG = EEG_aux.data(5,:);

   order=4;
   Lp=(20/(fs/2)); Hp=(40/(fs/2));
   [b,a] = butter(order,[Lp Hp],'bandpass');
   
   EMG_filt = filtfilt(b,a,EMG);
  
   rem_epochs = find(hypno2 == 'R');
   nrem_epochs = find(hypno2 == '2'| hypno2 == '3');
   n1_epochs = find(hypno2 == '1');
   wake_epochs = find(hypno2 == 'W');
  
   %% Make example figure

   
%    % example plot with RSN_001 and epoch 5
%    epoch = 5; % 42,46 muscle twitches, 52 has very clear alpha
%    
%    ep = rem_epochs(epoch);
%    samp_ep = (ep-1)*epochl*fs+1:ep*epochl*fs;
% 
%    waves_ep =waves_allepch_good(waves_allepch_good.nep == ep,:);
%    startsamp = waves_ep.startsamp;
%    endsamp = waves_ep.endsamp;
%    duration = waves_ep.duration;
%    frequency = waves_ep.frequency;
%    cycles = waves_ep.cycles;
%    SNR = waves_ep.SNR;
%    
%    bad_samps = setdiff(samp_ep,rem_goodsamp2);
%      
%    theta_ndx = find(frequency > lower_freq_theta & frequency < higher_freq_theta);
%    alpha_ndx = find(frequency > lower_freq_alpha & frequency < higher_freq_alpha);
% %    sigma_ndx = find(frequency > 12.5 & frequency < 16.5);
%    cycles_ndx = find(cycles > 3);
%    duration_ndx = find(duration > 0.3);
% 
%    theta_cycles_ndx = intersect(theta_ndx,cycles_ndx);
%    alpha_cycles_ndx = intersect(alpha_ndx,cycles_ndx);
% %    sigma_cycles_ndx = intersect(sigma_ndx,cycles_ndx);
%    
%    theta_cycles_duration_ndx = intersect(theta_cycles_ndx,duration_ndx);
%    alpha_cycles_duration_ndx = intersect(alpha_cycles_ndx,duration_ndx);
% %    sigma_cycles_duration_ndx = intersect(sigma_cycles_ndx,duration_ndx);
% 
% %    SNR_ndx = find(SNR > 1);
% %    alpha_cycles_SNR_ndx = intersect(alpha_cycles_ndx,SNR_ndx);
%    
%    
%    ch = 2;
%    
%    fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% 
%    subplot(4,1,1)
%    plot(samp_ep,data(ch,samp_ep),'k');
%    hold on
%    plot(bad_samps,data(ch,bad_samps),'Color',[0.7 0.7 0.7]) %orange
%    hold on
%    for c = 1:length(theta_cycles_duration_ndx)   
%    plot_samps_orig = round(startsamp(theta_cycles_duration_ndx(c)):endsamp(theta_cycles_duration_ndx(c)));
%    plot(plot_samps_orig,data(ch,plot_samps_orig),'Color',[0.4660 0.6740 0.1880]) %green
%    hold on
%    end
%    hold on
%    for c = 1:length(alpha_cycles_duration_ndx)   
%    plot_samps_orig = round(startsamp(alpha_cycles_duration_ndx(c)):endsamp(alpha_cycles_duration_ndx(c)));
%    plot(plot_samps_orig,data(ch,plot_samps_orig),'Color',[0.3010 0.7450 0.9330]) %clear blue
%    hold on
%    end
%    hold on
%    xticks((ep-1)*epochl*fs+1:2*fs:ep*epochl*fs);
%    xticklabels(0:2:30);
% %    for c = 1:length(sigma_cycles_duration_ndx)   
% %    plot_samps_orig = round(startsamp(sigma_cycles_duration_ndx(c)):endsamp(sigma_cycles_duration_ndx(c)));
% %    plot(plot_samps_orig,data(ch,plot_samps_orig),'Color',[0.8500 0.3250 0.0980]) %orange
% %    hold on
% %    end
%    title('Fz')
%    ylim([-50 50])
%    set(gca,'Fontsize',15);
%    ylabel('Amplitude (\muV)');
%   
%     subplot(4,1,2)
%     plot(samp_ep,lEOG(samp_ep),'Color',[1 0 0])
%     hold on
%     plot(samp_ep,rEOG(samp_ep),'Color',[0.6350 0.0780 0.1840])
%     ylim([-50 50])
%     title('EOG')
% %     legend({'lEOG' 'rEOG'})
%     xticks((ep-1)*epochl*fs+1:2*fs:ep*epochl*fs);
%     xticklabels(0:2:30);
%     set(gca,'Fontsize',15);
%     ylabel('Amplitude (\muV)');
%   
%     subplot(4,1,3)
%     plot(samp_ep,EMG_filt(samp_ep),'Color','b')
%     ylim([-50 50])
%     title('EMG')
%     xticks((ep-1)*epochl*fs+1:2*fs:ep*epochl*fs);
%     xticklabels(0:2:30);
%     set(gca,'Fontsize',15);
%     xlabel('Time (s)');
%     ylabel('Amplitude (\muV)');
%     
%     subplot(4,1,4)
%     plot(samp_ep,-ECG(samp_ep),'Color','b')
%     ylim([-200 200])
%     title('EMG')
%     xticks((ep-1)*epochl*fs+1:2*fs:ep*epochl*fs);
%     xticklabels(0:2:30);
%     set(gca,'Fontsize',15);
%     xlabel('Time (s)');
%     ylabel('Amplitude (\muV)');
% 
% 
% %    suptitle(sub_Folderpath(s).name)
% %    
% %    saveas(fig,[Savefolder,sub_Folderpath(s).name,'_waves_examples.svg']);
   
 %% Plot abundance and density across night
 
%  plot(1:length(hypno2),abundance_alpha,'k')
%  hold on
%  plot(wake_epochs,abundance_alpha(wake_epochs),'Color','b')
%  hold on
%  plot(n1_epochs,abundance_alpha(n1_epochs),'Color','cyan')
%  hold on
%  plot(nrem_epochs,abundance_alpha(nrem_epochs),'Color','m')
%  hold on
%  plot(rem_epochs,abundance_alpha(rem_epochs),'Color','green')
%  
%  plot(1:length(hypno2),abundance_theta,'k')
%  hold on
%  plot(wake_epochs,abundance_theta(wake_epochs),'Color','b')
%  hold on
%  plot(n1_epochs,abundance_theta(n1_epochs),'Color','cyan')
%  hold on
%  plot(nrem_epochs,abundance_theta(nrem_epochs),'Color','m')
%  hold on
%  plot(rem_epochs,abundance_theta(rem_epochs),'Color','green')
% 
%  plot(1:length(hypno2),density_alpha,'k')
%  hold on
%  plot(wake_epochs,density_alpha(wake_epochs),'Color','b')
%  hold on
%  plot(n1_epochs,density_alpha(n1_epochs),'Color','cyan')
%  hold on
%  plot(nrem_epochs,density_alpha(nrem_epochs),'Color','m')
%  hold on
%  plot(rem_epochs,density_alpha(rem_epochs),'Color','green')
%  
%  plot(1:length(hypno2),density_theta,'k')
%  hold on
%  plot(wake_epochs,density_theta(wake_epochs),'Color','b')
%  hold on
%  plot(n1_epochs,density_theta(n1_epochs),'Color','cyan')
%  hold on
%  plot(nrem_epochs,density_theta(nrem_epochs),'Color','m')
%  hold on
%  plot(rem_epochs,density_theta(rem_epochs),'Color','green')
%  
 %%
abundance_rem_alpha(s) = nanmean(abundance_alpha(rem_epochs)); 
abundance_rem_theta(s) = nanmean(abundance_theta(rem_epochs)); 

abundance_nrem_alpha(s) = nanmean(abundance_alpha(nrem_epochs)); 
abundance_nrem_theta(s) = nanmean(abundance_theta(nrem_epochs)); 

density_rem_alpha(s) = nanmean(density_alpha(rem_epochs)); 
density_rem_theta(s) = nanmean(density_theta(rem_epochs)); 

density_nrem_alpha(s) = nanmean(density_alpha(nrem_epochs)); 
density_nrem_theta(s) = nanmean(density_theta(nrem_epochs)); 

waves_ep_3cyc = waves_allepch_good(waves_allepch_good.cycles > 3,:);
waves_ep_3cyc_dur = waves_ep_3cyc(waves_ep_3cyc.duration > 0.3,:);
waves_ep_3cyc_dur_alphatheta = waves_ep_3cyc_dur(waves_ep_3cyc_dur.frequency > 0 & waves_ep_3cyc_dur.frequency < 12.5,:);
waves_ep_3cyc_dur_alpha = waves_ep_3cyc_dur(waves_ep_3cyc_dur.frequency > lower_freq_alpha & waves_ep_3cyc_dur.frequency < higher_freq_alpha,:);
waves_ep_3cyc_dur_theta = waves_ep_3cyc_dur(waves_ep_3cyc_dur.frequency > lower_freq_theta & waves_ep_3cyc_dur.frequency < higher_freq_theta,:);

waves_rem = waves_ep_3cyc_dur_alphatheta(waves_ep_3cyc_dur_alphatheta.stage == 'R',:);
waves_rem_alpha = waves_ep_3cyc_dur_alpha(waves_ep_3cyc_dur_alpha.stage == 'R',:);
waves_rem_theta = waves_ep_3cyc_dur_theta(waves_ep_3cyc_dur_theta.stage == 'R',:);

% mdn_wake_alpha(s) = nanmedian(waves_wake_alpha.frequency);
mdn_rem_alpha(s) = nanmedian(waves_rem_alpha.frequency);
m_rem_alpha(s) = nanmean(waves_rem_alpha.frequency);
n_rem_alpha(s) = length(waves_rem_alpha.frequency);


% mdn_wake_theta(s) = nanmedian(waves_wake_theta.frequency);
mdn_rem_theta(s) = nanmedian(waves_rem_theta.frequency);
m_rem_theta(s) = nanmean(waves_rem_theta.frequency);
n_rem_theta(s) = length(waves_rem_theta.frequency);

for ep = 1:nepochs
    
    waves_ep_ndx = find(waves_allepch.nep == ep);
    offset_allep(ep) = waves_allepch.offset(waves_ep_ndx(1));
    exponent_allep(ep) = waves_allepch.exponent(waves_ep_ndx(1));
    
    clear waves_ep_ndx
    
end

m_offset_wake(s) = nanmean(offset_allep(find(hypno2 == 'W')));
m_offset_n1(s) = nanmean(offset_allep(find(hypno2 == '1')));
m_offset_n2(s) = nanmean(offset_allep(find(hypno2 == '2')));
m_offset_n3(s) = nanmean(offset_allep(find(hypno2 == '3')));
m_offset_rem(s) = nanmean(offset_allep(find(hypno2 == 'R')));

m_exponent_wake(s) = nanmean(exponent_allep(find(hypno2 == 'W')));
m_exponent_n1(s) = nanmean(exponent_allep(find(hypno2 == '1')));
m_exponent_n2(s) = nanmean(exponent_allep(find(hypno2 == '2')));
m_exponent_n3(s) = nanmean(exponent_allep(find(hypno2 == '3')));
m_exponent_rem(s) = nanmean(exponent_allep(find(hypno2 == 'R')));


% fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
% % subplot(1,4,1)
% histogram(waves_wake.frequency,'FaceColor','b')
% hold on
% xline(mdn_wake_alpha(s),'Color','b','LineWidth',2)
% hold on
% xline(mdn_wake_theta(s),'Color','b','LineWidth',2)
% % title('Wake');

%% Stim frequency alpha and theta

mdn_ISI_all = NaN(18,1,1);
mdn_ISI_alpha = NaN(18,1,1);
mdn_ISI_theta = NaN(18,1,1);

% for s = 1:19
    
    ISI_alphastim_all_allcon = [];
    ISI_alphastim_alpha_allcon = [];
    
    for con = 1:4
        
        ISI_subcon = ISI_allsub(s).con{con};
        ISI_freq = 1./ISI_subcon;
        
        ISI_all = ISI_freq;
        ISI_alpha = ISI_freq(find(ISI_freq > lower_freq_alpha & ISI_freq < higher_freq_alpha));
        
        ISI_alphastim_all_allcon = vertcat(ISI_alphastim_all_allcon,ISI_all');
        ISI_alphastim_alpha_allcon = vertcat(ISI_alphastim_alpha_allcon,ISI_alpha');
        
        clear ISI_all ISI_alpha ISI_subcon ISI_freq
      
    end
    
    ISI_thetastim_all_allcon = [];
    ISI_thetastim_theta_allcon = [];
    
      for con = 5:8
        
        ISI_subcon = ISI_allsub(s).con{con};
        ISI_freq = 1./ISI_subcon;
        
        ISI_all = ISI_freq;
        ISI_theta = ISI_freq(find(ISI_freq > lower_freq_theta & ISI_freq < higher_freq_theta));
        
        ISI_thetastim_all_allcon = vertcat(ISI_thetastim_all_allcon,ISI_all');
        ISI_thetastim_theta_allcon = vertcat(ISI_thetastim_theta_allcon,ISI_theta');
        
        clear ISI_all ISI_theta ISI_subcon ISI_freq
      
    end
    
     m_ISI_alphastim_all(s) = nanmean(ISI_alphastim_all_allcon);
     m_ISI_alphastim_alpha(s) = nanmean(ISI_alphastim_alpha_allcon);
     m_ISI_thetastim_all(s) = nanmean(ISI_thetastim_all_allcon);
     m_ISI_thetastim_theta(s) = nanmean(ISI_thetastim_theta_allcon);
     
     mdn_ISI_alphastim_all(s) = nanmedian(ISI_alphastim_all_allcon);
     mdn_ISI_alphastim_alpha(s) = nanmedian(ISI_alphastim_alpha_allcon);
     mdn_ISI_thetastim_all(s) = nanmedian(ISI_thetastim_all_allcon);
     mdn_ISI_thetastim_theta(s) = nanmedian(ISI_thetastim_theta_allcon);

      
%     clear ISI_all_allcon ISI_alpha_allcon ISI_theta_allcon

% end

%% Oscillation frequency - findpeaks

colors = linspecer(8);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

% subplot(2,1,1)
histogram(waves_rem.frequency,30,'FaceColor',colors(5,:),'FaceAlpha',1)
[f,xi] = ksdensity(waves_rem.frequency); 
hold on
plot(xi,f*2500,'Color',[0.5216 0.3686 0.0627],'LineWidth',5)
hold on

[pks,locs] = findpeaks(f);


% Find the half max value.
halfMax = (min(f) + max(f)) / 2;
% Find where the data first drops below half the max.
index1 = find(f >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(f >= halfMax, 1, 'last');
% fwhm = index2-index1 + 1; % FWHM in indexes.
% OR, if you have an x vector
fwhmx_osc(s) = xi(index2) - xi(index1);


der = diff(f);
der_ndx = find(der < 0.001 & der > -0.001);
% plot(xi(der_ndx),f(der_ndx),'o');

alpha_peak_findpeaks(s) = NaN;
% theta_peak_findpeaks(s) = NaN;


 if find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha)
            
     [row,col] = find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha);
        
            if length(col)>1
                [r,c] = find(pks(col)==max(pks(col)));
                alpha_peak_findpeaks(s) = xi(locs(col(r)));
                scatter(alpha_peak_findpeaks(s),pks(col(r))*2500,700,'o','MarkerFaceColor',[0.5216 0.3686 0.0627],'MarkerEdgeColor','k','LineWidth',3)
            else
                alpha_peak_findpeaks(s) = xi(locs(col));
                scatter(alpha_peak_findpeaks(s),pks(col)*2500,700,'o','MarkerFaceColor',[0.5216 0.3686 0.0627],'MarkerEdgeColor','k','LineWidth',3)
            end
 end
  
  
%  saveas(fig,[Savefolder,sub_Folderpath(s).name,'_oscillation_freq_findpeaks.svg']);


% Alpha stim frequency - findpeaks

% fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

% subplot(2,1,2)
% uisetcolor()
histogram(ISI_alphastim_all_allcon,30,'FaceColor',colors(7,:),'FaceAlpha',0.5)
[f,xi] = ksdensity(ISI_alphastim_all_allcon); 
hold on
plot(xi,f*4000,'Color',colors(7,:),'LineWidth',5)
hold on

[pks,locs] = findpeaks(f);


% Find the half max value.
halfMax = (min(f) + max(f)) / 2;
% Find where the data first drops below half the max.
index1 = find(f >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(f >= halfMax, 1, 'last');
% fwhm = index2-index1 + 1; % FWHM in indexes.
% OR, if you have an x vector
fwhmx_alphastim(s) = xi(index2) - xi(index1);


der = diff(f);
der_ndx = find(der < 0.001 & der > -0.001)
% plot(xi(der_ndx),f(der_ndx),'o');

alphastim_alpha_peak_findpeaks(s) = NaN;
% alphastim_theta_peak_findpeaks(s) = NaN;


 if find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha)
            
     [row,col] = find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha);
        
            if length(col)>1
                [r,c] = find(pks(col)==max(pks(col)));
                alphastim_alpha_peak_findpeaks(s) = xi(locs(col(r)));
                scatter(alphastim_alpha_peak_findpeaks(s),pks(col(r))*4000,700,'o','MarkerFaceColor',colors(7,:),'MarkerEdgeColor','k','LineWidth',3)
            else
                alphastim_alpha_peak_findpeaks(s) = xi(locs(col));
                scatter(alphastim_alpha_peak_findpeaks(s),pks(col)*4000,700,'o','MarkerFaceColor',colors(7,:),'MarkerEdgeColor','k','LineWidth',3)
            end
 end
 
     
axis square
box off
xlim([0 14]);
xticks(0:2:14);
set(gca,'FontSize',35);
xlabel('Frequency (Hz)');
ylabel('Number')

% saveas(fig,[Savefolder,sub_Folderpath(s).name,'_alphastim_freq_findpeaks.svg']);


%% Theta stim frequency - findpeaks

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

colors = linspecer(8);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])


% subplot(2,1,1)
histogram(waves_rem.frequency,30,'FaceColor',colors(5,:),'FaceAlpha',1)
[f,xi] = ksdensity(waves_rem.frequency); 
hold on
plot(xi,f*2500,'Color',[0.5216 0.3686 0.0627],'LineWidth',5)
hold on

[pks,locs] = findpeaks(f);

der = diff(f);
der_ndx = find(der < 0.001 & der > -0.001);
% plot(xi(der_ndx),f(der_ndx),'o');

% alpha_peak_findpeaks(s) = NaN;
theta_peak_findpeaks(s) = NaN;


%  if find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha)
%             
%      [row,col] = find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha);
%         
%             if length(col)>1
%                 [r,c] = find(pks(col)==max(pks(col)));
%                 alpha_peak_findpeaks(s) = xi(locs(col(r)));
%                 scatter(alpha_peak_findpeaks(s),pks(col(r))*2500,700,'o','MarkerFaceColor',[0.5216 0.3686 0.0627],'MarkerEdgeColor','k','LineWidth',3)
%             else
%                 alpha_peak_findpeaks(s) = xi(locs(col));
%                 scatter(alpha_peak_findpeaks(s),pks(col)*2500,700,'o','MarkerFaceColor',[0.5216 0.3686 0.0627],'MarkerEdgeColor','k','LineWidth',3)
%             end
%  end

 if find(xi(locs)>=lower_freq_theta & xi(locs)<=higher_freq_theta)
            
     [row,col] = find(xi(locs)>=lower_freq_theta & xi(locs)<=higher_freq_theta);
        
            if length(col)>1
                [r,c] = find(pks(col)==max(pks(col)));
                theta_peak_findpeaks(s) = xi(locs(col(r)));
                scatter(theta_peak_findpeaks(s),pks(col(r))*2500,700,'o','MarkerFaceColor',[0.5216 0.3686 0.0627],'MarkerEdgeColor','k','LineWidth',3)
            else
                theta_peak_findpeaks(s) = xi(locs(col));
                scatter(theta_peak_findpeaks(s),pks(col)*2500,700,'o','MarkerFaceColor',[0.5216 0.3686 0.0627],'MarkerEdgeColor','k','LineWidth',3)
            end
 end


histogram(ISI_thetastim_all_allcon,30,'FaceColor',colors(7,:),'FaceAlpha',0.5)
[f,xi] = ksdensity(ISI_thetastim_all_allcon); 
hold on
plot(xi,f*2000,'Color',colors(7,:),'LineWidth',5)
hold on

[pks,locs] = findpeaks(f);


% Find the half max value.
halfMax = (min(f) + max(f)) / 2;
% Find where the data first drops below half the max.
index1 = find(f >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(f >= halfMax, 1, 'last');
% fwhm = index2-index1 + 1; % FWHM in indexes.
% OR, if you have an x vector
fwhmx_thetastim(s) = xi(index2) - xi(index1);


der = diff(f);
der_ndx = find(der < 0.001 & der > -0.001)
% plot(xi(der_ndx),f(der_ndx),'o');

% thetastim_alpha_peak_findpeaks(s) = NaN;
thetastim_theta_peak_findpeaks(s) = NaN;


%  if find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha)
%             
%      [row,col] = find(xi(locs)>=lower_freq_alpha & xi(locs)<=higher_freq_alpha);
%         
%             if length(col)>1
%                 [r,c] = find(pks(col)==max(pks(col)));
%                 thetastim_alpha_peak_findpeaks(s) = xi(locs(col(r)));
%                 scatter(thetastim_alpha_peak_findpeaks(s),pks(col(r)),50,'r','filled')
%             else
%                 thetastim_alpha_peak_findpeaks(s) = xi(locs(col));
%                 scatter(thetastim_alpha_peak_findpeaks(s),pks(col),50,'r','filled')
%             end
%  end
%  
  
 if find(xi(locs)>=lower_freq_theta & xi(locs)<=higher_freq_theta)
            
     [row,col] = find(xi(locs)>=lower_freq_theta & xi(locs)<=higher_freq_theta);
        
            if length(col)>1
                [r,c] = find(pks(col)==max(pks(col)));
                thetastim_theta_peak_findpeaks(s) = xi(locs(col(r)));
                scatter(thetastim_theta_peak_findpeaks(s),pks(col(r))*2000,700,'o','MarkerFaceColor',colors(7,:),'MarkerEdgeColor','k','LineWidth',3)
            else
                thetastim_theta_peak_findpeaks(s) = xi(locs(col));
                scatter(thetastim_theta_peak_findpeaks(s),pks(col)*2000,700,'o','MarkerFaceColor',colors(7,:),'MarkerEdgeColor','k','LineWidth',3)
            end
 end
   
axis square
box off
xlim([0 14]);
xticks(0:2:14);
set(gca,'FontSize',35);
xlabel('Frequency (Hz)');
ylabel('Number')
 
%    saveas(fig,[Savefolder,sub_Folderpath(s).name,'_thetastim_freq_findpeaks.svg']);



end

%%
save([Savefolder,'Oscillation_Freq_ISI_allsub_d03_',date,'.mat']);
% save([Savefolder,'Oscillation_Freq_ISI_allsub_30_45_Hz_',date,'.mat']);

