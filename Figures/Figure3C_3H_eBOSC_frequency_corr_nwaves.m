clear all;
close all;

%%

load('D:\Valeria\RSN\data\for_sharing\data_to_make_figures\Oscillation_Freq_ISI_allsub_d03_08-Feb-2024.mat');
close

Savefolder = 'D:\Valeria\RSN\data\for_sharing\data_to_make_figures\Figures\';

%%
incl_sub = setdiff(1:19,[12]);

nonan_sub_alpha = intersect(find(~isnan(alpha_peak_findpeaks)==1),incl_sub);
nonan_sub_alphastim_alpha = intersect(find(~isnan(alphastim_alpha_peak_findpeaks)==1),incl_sub);
nonan_sub_alpha = intersect(nonan_sub_alpha,nonan_sub_alphastim_alpha);
length(nonan_sub_alpha)

[r p] = corr(alpha_peak_findpeaks(nonan_sub_alpha)',alphastim_alpha_peak_findpeaks(nonan_sub_alpha)')

   fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
%     scatter(mch_allcon_ifq_alpha_off(incl_sub),nanmean(m_ISI_theta(incl_sub,5:8),2),100,'k','LineWidth',3);
    scatter(alpha_peak_findpeaks(nonan_sub_alpha),alphastim_alpha_peak_findpeaks(nonan_sub_alpha),300,'MarkerEdgeColor','k','LineWidth',3);
    l = lsline;
    l.LineWidth = 3;
    l.Color = 'k';
    xlabel('REM Alpha Freq (Hz)');
    ylabel('Alpha ISI Freq (Hz)');
%     title(['OFF, r = ', num2str(r), ', p = ',num2str(p)]);
    set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
    box off
    axis square
    xlim([7 11]);
    ylim([7 11]);
    xticks(7.0:1.0:11.0);
    xticklabels(7.0:1.0:11.0);
    xtickformat('%.1f')
    yticks(7.0:1.0:11.0);
    yticklabels(7.0:1.0:11.0);
    ytickformat('%.1f')

    saveas(fig,[Savefolder,'Figure3C_corr_REM_alpha_oscillations_ISI_Alpha_',num2str(lower_freq_alpha),'_',num2str(higher_freq_alpha),'.svg']);

%%


incl_sub = setdiff(1:19,[12]);

nonan_sub_theta = intersect(find(~isnan(theta_peak_findpeaks)==1),incl_sub);
nonan_sub_thetastim_theta = intersect(find(~isnan(thetastim_theta_peak_findpeaks)==1),incl_sub);
nonan_sub_theta = intersect(nonan_sub_theta,nonan_sub_thetastim_theta);
length(nonan_sub_theta)

[r p] = corr(theta_peak_findpeaks(nonan_sub_theta)',thetastim_theta_peak_findpeaks(nonan_sub_theta)')

   fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 1 1])
%     scatter(mch_allcon_ifq_alpha_off(incl_sub),nanmean(m_ISI_theta(incl_sub,5:8),2),100,'k','LineWidth',3);
    scatter(theta_peak_findpeaks(nonan_sub_theta)',thetastim_theta_peak_findpeaks(nonan_sub_theta)',300,'MarkerEdgeColor','k','LineWidth',3);
    l = lsline;
    l.LineWidth = 3;
    l.Color = 'k';
    xlabel('REM Theta Freq (Hz)');
    ylabel('Theta ISI Freq (Hz)');
%     title(['OFF, r = ', num2str(r), ', p = ',num2str(p)]);
    set(gca,'Fontsize',35,'TickDir','out','LineWidth',3);
    box off
    axis square
    xlim([4.5 8]);
    ylim([4.5 8]);
    xticks(4.5:1:8);
    xticklabels(4.5:1:8);
    xtickformat('%.1f')
    yticks(4.5:1:8);
    yticklabels(4.5:1:8);
    ytickformat('%.1f')


    saveas(fig,[Savefolder,'Figure3H_corr_REM_alpha_oscillations_ISI_Alpha_',num2str(lower_freq_alpha),'_',num2str(higher_freq_alpha),'.svg']);

%% number of waves

m_n_rem_alpha = round(nanmean(n_rem_alpha(incl_sub)))
% sem_n_rem_alpha = round(nanstd(n_rem_alpha(incl_sub))/sqrt(length(incl_sub)))
sd_n_rem_alpha = round(nanstd(n_rem_alpha(incl_sub)))

m_n_rem_theta = round(nanmean(n_rem_theta(incl_sub)))
% sem_n_rem_theta = round(nanstd(n_rem_theta(incl_sub))/sqrt(length(incl_sub)))
sd_n_rem_theta = round(nanstd(n_rem_theta(incl_sub)))


m_abundance_rem_alpha = nanmean(abundance_rem_alpha(incl_sub))
% sem_abundance_rem_alpha = nanstd(abundance_rem_alpha(incl_sub))/sqrt(length(incl_sub))
sd_abundance_rem_alpha = nanstd(abundance_rem_alpha(incl_sub))

m_abundance_rem_theta = nanmean(abundance_rem_theta(incl_sub))
% sem_abundance_rem_theta = nanstd(abundance_rem_theta(incl_sub))/sqrt(length(incl_sub))
sd_abundance_rem_theta = nanstd(abundance_rem_theta(incl_sub))


% m_density_rem_alpha = nanmean(density_rem_alpha(incl_sub))
% sem_density_rem_alpha = nanstd(density_rem_alpha(incl_sub))/sqrt(length(incl_sub))
% 
% m_density_rem_theta = nanmean(density_rem_theta(incl_sub))
% sem_density_rem_theta = nanstd(density_rem_theta(incl_sub))/sqrt(length(incl_sub))
% 
%% FWHM

m_fwhm_rem_osc = nanmean(fwhmx_osc(incl_sub))
sd_fwhm_rem_osc = nanstd(fwhmx_osc(incl_sub))

m_fwhm_rem_alphastim = nanmean(fwhmx_alphastim(incl_sub))
sd_fwhm_rem_alphastim = nanstd(fwhmx_alphastim(incl_sub))

m_fwhm_rem_thetastim = nanmean(fwhmx_thetastim(incl_sub))
sd_fwhm_rem_thetastim = nanstd(fwhmx_thetastim(incl_sub))


