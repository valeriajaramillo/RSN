# RSN: REM sleep neuromodulation


## Overview

This repository contains scripts for the paper:


Raw data and data to create the figures can be found here:
[Zenodo DOI](10.5281/zenodo.10663994)


Required toolboxes to run the scripts: 

- [eeglab v2021.1](https://sccn.ucsd.edu/eeglab/download.php)
- [fieldtrip revision ebc13229d](https://www.fieldtriptoolbox.org/) 
- [circular statistics](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
- [eBOSC](https://github.com/jkosciessa/eBOSC)


## Preprocessing

### sleep

**A1_writeedf_for_scoring**: write .edf file that is used for sleep scoring (re-referencing 8 EEG channels used for sleep scoring to mastoids), you need to export .eeg file as .edf file using BrainVisionRecorder before running this script (.edf file is needed)  

**A2_extract_triggers**: extract AEP and CLAS triggers from .eeg file (raw data, provided), save in _nm.mat file  

**B1_load_filter_data**: load and filter sleep EEG data from .eeg file (raw data, provided), save as _fil.fdt and .set file  

**C1_artcorr_extract_REM**: load filtered data, sleep scoring, and phasic/tonic scoring data, visualize REM sleep EEG using eeglab, visually identify and interpolate bad channels, visually mark segments that contain artefacts, save _goodREM.fdt and .set files that only contain 'good' (artefact-free) REM samples, save _goodREM.mat file containing information about sleep and phasic/tonic scoring, bad channels, and artefacts  

**D_ICA**: perform ICA on _goodREM data using eeglab using AMICA1.6 plugin  

**E_manual_ICA**: visualize independent components using eeglab, visually identify eye movement, cardiac, muscle activity, and channel noise components and remove those components, save as _manual_ICA.fdt and .set files  

**F_insert_REM_after_ICA_rereferencing**: put back cleaned REM sleep data (after ICA) into the whole-night data and re-reference to the average across all channels  

**G_goodtrigs**: save triggers that occurred during a 'good' REM sample in _nm_good.mat  


### wake

**A_load_filter_data_wake**: load and filter wake EEG data from .eeg file (raw data, provided), save as _fil.fdt and .set file  

**B_extract_triggers_wake**: extract AEP triggers from .eeg file (raw data, provided), save in _nm.mat file  

**C_artcorr_wake**: load filtered data, visualize EEG using eeglab, visually identify and interpolate bad channels, visually mark segments that contain artefacts, save _goodwake.fdt and .set files that only contain 'good' (artefact-free) wake samples, save _goodwake.mat file containing information about bad channels, and artefacts  

**D_ICA_wake**: perform ICA on _goodwake data using eeglab using AMICA1.6 plugin  

**E_manual_ICA**: visualize independent components using eeglab, visually identify eye movement, cardiac, muscle activity, and channel noise components and remove those components, save as _manual_ICA.fdt and .set files  

**F_insert_wake_after_ICA_rereferencing**: re-reference wake data to the average across all channels  

**G_goodtrigs**: save triggers that occurred during a 'good' wake sample in _nm_good.mat  


## Analyses


### connectivity: perform connectivity analyses, run in order A-D

**A_phase_connectivity**: calculate the phase using Hilbert transform for trial data  
**B_PLV_sub**: calculate PLV and PLI for all channel pairs for each participant  
**C_PLV_allsub**: save PLV and PLI for all participants  
**D_cluster_analysis_SnPM_connectivity**: calculate lme for each electrode with cluster correction for alpha CLAS  
**D_cluster_analysis_SnPM_connectivity_theta**: calculate lme for each electrode with cluster correction for theta CLAS  

### eBOSC: detect oscillations to determine individual peak frequency, run in order A-B

**A_eBOSC**: detect oscillations for each epoch  
**B_eBOSC**: extract oscillations with at least 3 cycles and 300 ms duration, plot histogram and detect peak  

### ERP: perform AEP analyses for REM and wake data, run in order A-B**

**A_ERP_sub_REM**: extract AEP REM data for each participant  
**A_ERP_sub_wake**: extract AEP wake data for each participant  
**B_ERP_allsub_REM**: save AEP REM data for all participants  
**B_ERP_allsub_wake**: save AEP wake data for all participants  

### frequency: perform frequency analyses, run in order A-C

**A_instantaneous_frequency**: calculate instantaneous frequency for REM and wake data for each participant  
**B_freq_allsub**: calculate frequency for ON and OFF windows and save for all participants  
**C_cluster_analysis_SnPM_frequency**:  calculate lme for each electrode with cluster correction for alpha CLAS  
**C_cluster_analysis_SnPM_frequency_theta**: calculate lme for each electrode with cluster correction for theta CLAS  

### ISI: calculate inter-stimulus-intervals (ISI)

**A_ISI_allsub_allch**: calculate ISI's and save for all participants  

### psd: perform power analyses

**A1_Power_ON_OFF_goodREM**: perform power and phase analysis for all trials of each condition and for each participant  
**B_Power_ON_OFF_allsub_allch**: average across all good trials, all good phasic and all good tonic trials and save for all participants  
**C_cluster_analysis_SnPM_power**: calculate lme for each electrode with cluster correction for alpha CLAS  
**C_cluster_analysis_SnPM_power_theta**: calculate lme for each electrode with cluster correction for theta CLAS  

### resultant: perform phase-locking accuracy analyses

**A_resultant_allsub_allch**: calculate resultant and mean phase for each participant, condition and channel and save for all participants  

### sleep_parameters: calculate sleep stage and phasic/tonic percentages

**A_Table1_Sleep_parameters**: calculate classical sleep parameters  
**B_Table1_phasic_tonic_conditions**: calculate percentage of phasic and tonic REM sleep for each participant and for ON and OFF windows for each condition, save for all participants  


## Figures

**Figure1C_phasic_tonic_epoch**: plot example 10-s data for phasic and tonic REM sleep  
**Figure1D_Suppl_Figure1_ntrials_across_night_erps**: make boxplots for number of AEP and CLAS trials (all and phasic)  
**Figure2_Suppl_Figure2_phasic_tonic_psd_AEPs**: plot power spectrum and AEPs for phasic, tonic and wake EO (eyes open) and EC (eyes closed)  
**Figure3A_3F_ind_psd**: plot power spectrum for each individual in alpha and theta range  
**Figure3B_3G_eBOSC_frequency_histogram**: plot histogram of detected oscillations and stimulation frequencies for example participant  
**Figure3C_3H_eBOSC_frequency_corr_nwaves**: plot correlation scatter plot between individual oscillation peak frequency and peak stimulation frequency  
**Figure3D_3I_polarhistogram**: plot polar histogram of mean phase and resultant across participants  
**Figure3E_3J_resultant_topo**: plot resultant for each electrode across scalp, test non-uniformity across circle  
**Figure4A_4G_Suppl_Figure5-9and12-17_Power_ON_OFF_allsub_topo_boxplots**: plot topography of power change lme F-values and significant electrode clusters  
**Figure4B_4C_Power_ON_OFF_change_psd_phasictonic_alpha**: plot power changes across spectrum for alpha CLAS  
**Figure4D_4J_Suppl_Figure10and18_Frequency_ON_OFF_allsub_topo**: plot topography of frequency change lme F-values and significant electrode clusters  
**Figure4E-F_4K-L_Frequency_ON_OFF_change_boxplots_time**: make boxplots for frequency change and plot frequency change over time across different conditions  
**Figure4H_4I_Power_ON_OFF_change_psd_phasictonic_theta**: plot power changes across spectrum for theta CLAS  
**Suppl_Figure3_ISI_autocorr**: plot autocorrelation across ISI's  
**Suppl_Figure11_connectivity_alpha**: plot topography of connectivity change lme F-values and significant electrode clusters for alpha CLAS  
**Suppl_Figure19_connectivity_theta**: plot topography of connectivity change lme F-values and significant electrode clusters for theta CLAS  
**Suppl_Table1_questionnaires**: make table for KSS and VAMS changes from evening to morning  








