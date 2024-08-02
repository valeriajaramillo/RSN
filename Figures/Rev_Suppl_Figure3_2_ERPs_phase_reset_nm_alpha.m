clear all;
close all;

addpath(genpath('/users/nemo/software/eeglab'));
addpath(genpath('/users/nemo/projects/RSN'));
addpath(genpath('/users/nemo/software/Henry/useful_functions'));
% addpath(genpath('/user/HS301/m17462/matlab/colorGradient'));

Savefolder = '/parallel_scratch/nemo/RSN/analysis/analysis/Figures/';

% incl_sub = setdiff(1:19,[12]); % 14 excluded because no phasic trials, 9 because no wake eve trials
incl_sub = setdiff(1:19,[12]); % 14 excluded because no phasic trials, 9 because no wake eve trials


%% ERP phase bins - REM

% load('/parallel_scratch/nemo/RSN/analysis/analysis/erp_allsub/ERP_nm_allsub_REM_mICA_avref09-Feb-2024.mat');
load('/parallel_scratch/nemo/RSN/analysis/analysis/erp_allsub/ERP_nm_broadband_allsub_REM_mICA_avref07-Jun-2024.mat');

%% plot alpha ERP - averaged across triggers

t = (-dt*fs:dt*fs-1)/fs*1000; % time in ms

colors = linspecer(4);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
fig.WindowState = 'maximized';

% tileplot = tiledlayout(1,4);
% 
% nexttile(1)
% for trig = 1:10
trig = 1:45;



% trig = 1;

% subplot(2,5,trig)
for con = 1:4
    
%     for s = 1:size(ERP_nm_all.trial_data,1)
%         nonan_trig = find(~isnan(squeeze(ERP_nm_all.trial_data(s,con,:,1))) == 1);
%         last_trig(s) = nonan_trig(end);
%     end

plot(t,squeeze(nanmean(nanmean(ERP_nm_all.trial_data(incl_sub,con,trig,:),1),3)),'Color',colors(con,:),'LineWidth',2)
hold on

% m_ERP_nm_con = squeeze(nanmean(nanmean(ERP_nm_all.trial_data(incl_sub,con,trig,:),1),3));
% sem_ERP_nm_con = squeeze(nanmean(nanstd(ERP_nm_all.trial_data(incl_sub,con,trig,:),1),3))./sqrt(length(incl_sub));
% 
% shadedErrorBar(t,m_ERP_nm_con,sem_ERP_nm_con,'lineProps',{'Color',colors(con,:),'LineWidth',3},'patchSaturation',.3);
% hold on
% 
% clear m_ERP_nm_con sem_ERP_nm_con

end


hold on
plot(t,squeeze(nanmean(nanmean(nanmean(ERP_nm_all.trial_data(incl_sub,1:4,trig,:),3),2),1)),'Color','k','LineWidth',3)

xline(0,'LineStyle','--','LineWidth',5);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'Peak' 'Falling' 'Trough' 'Rising'});
legend off
ylim([-3 4]);
xticks(-200:200:1000);
% title(['Stimulus ',num2str(trig)]);
xtickangle(0)
set(groot,'defaultAxesXTickLabelRotationMode','manual')

% end

saveas(fig,[Savefolder,'Suppl_Figure4_ERP_nm.svg']);

%% plot alpha ERP - averaged across first 10 triggers - alphafilt

t = (-dt*fs:dt*fs-1)/fs*1000; % time in ms

colors = linspecer(4);

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
fig.WindowState = 'maximized';

% tileplot = tiledlayout(1,4);
% 
% nexttile(1)
% for trig = 1:10
trig = 1:45;
% subplot(2,5,trig)
for con = 1:4
plot(t,squeeze(nanmean(nanmean(ERP_nm_all.trial_data_alphafilt_real(incl_sub,con,trig,:),1),3)),'Color',colors(con,:),'LineWidth',2)
hold on
end

hold on
plot(t,squeeze(nanmean(nanmean(nanmean(ERP_nm_all.trial_data_alphafilt_real(incl_sub,1:4,:,:),3),2),1)),'Color','k','LineWidth',3)

xline(0,'LineStyle','--','LineWidth',2);

xlim([-200 1000]);
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
box off
axis square
% legend({'Peak' 'Falling' 'Trough' 'Rising'});
legend off
ylim([-3 3]);
xticks(-200:200:1000);
% title(['Stimulus ',num2str(trig)]);
xtickangle(0)
set(groot,'defaultAxesXTickLabelRotationMode','manual')

% end
saveas(fig,[Savefolder,'Suppl_Figure4_ERP_nm_alphafilt.svg']);


%% check non-uniformity across circle - alphafilt

fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])

trig = 1;

for samp = 1:size(ERP_nm_all.trial_data_alphafilt_phase,4)
        
    mean_phase_allcon_REM = reshape(ERP_nm_all.trial_data_alphafilt_phase(incl_sub,1:4,trig,samp),[length(incl_sub)*4 1]);
    resultant_phase_allcon_REM = reshape(ERP_nm_all.trial_data_alphafilt_phase_r(incl_sub,1:4,trig,samp),[length(incl_sub)*4 1]);
    
    mean_phase_allcon_REM(isnan(mean_phase_allcon_REM))= [];
    resultant_phase_allcon_REM(resultant_phase_allcon_REM == 0)= [];
    
    [p_val_REM(samp) z_val_REM(samp)] = circ_rtest(mean_phase_allcon_REM,resultant_phase_allcon_REM);
    
end



% plot alpha ERP phase

t = (-dt*fs:dt*fs-1)/fs*1000; % time in ms

colors = linspecer(4);


% tileplot = tiledlayout(1,4);
% 
% nexttile(1)

% subplot(2,5,trig)    
    
for con = 1:4
% plot(t,squeeze(nanmean(ERP_nm_all.trial_data(:,con,trig,:),1)),'Color',colors(con,:))
incl_sub2 = find(~isnan(squeeze(ERP_nm_all.trial_data_alphafilt_phase(:,con,trig,1))) == 1);
incl_sub3 = intersect(incl_sub,incl_sub2);
plot(t,squeeze(circ_mean(ERP_nm_all.trial_data_alphafilt_phase(incl_sub3,con,trig,:))),'Color',colors(con,:))
hold on
% sig_samps = find(p_val_REM < 0.05);
% plot(t(sig_samps),ones(length(sig_samps),1)*3.5,'*','Color','k');
% hold on
end

hold on
% plot(t,squeeze(circ_mean(nanmean(ERP_nm_all.trial_data_alphafilt_phase(incl_sub,1:4,trig,:),2))),'Color','k','LineWidth',3)
xline(0,'LineStyle','--','LineWidth',2);

xlim([-500 500]);
xlabel('Time (ms)')
ylabel('Phase (radians)')
set(gca,'Fontsize',15,'TickDir','out','LineWidth',2);
box off
axis square
% legend({'Peak' 'Falling' 'Trough' 'Rising'});
legend off
ylim([-4 4]);
yticks([-pi:pi:pi]);
yticklabels({'-\pi' '0' '\pi'})
xticks(-200:200:1000);
title(['Stimulus ',num2str(trig)]);

clear p_val_REM z_val_REM
% end

% saveas(fig,[Savefolder,'Figure6_ERP_alpha_nm_alphafilt_phase.svg']);







