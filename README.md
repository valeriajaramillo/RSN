# RSN

[Code from](https://github.com/valeriajaramillo/RSN/)

required toolboxes: 

- [eeglab v2021.1](https://sccn.ucsd.edu/eeglab/download.php)
- [fieldtrip revision ebc13229d](https://www.fieldtriptoolbox.org/) 


# preprocessing

sleep: preprocessing scripts for sleep data, run in order A-G

1. A1_writeedf_for_scoring: write .edf file that is used for sleep scoring (re-referencing 8 EEG channels used for sleep scoring to mastoids), you need to export .eeg file as .edf file using BrainVisionRecorder before running this script (.edf file is needed)

2. A2_extract_triggers: extract AEP and CLAS triggers from .eeg file (raw data, provided), save in _nm.mat file

3. B1_load_filter_data: load and filter sleep EEG data from .eeg file (raw data, provided), save as _fil.fdt and .set file

4. C1_artcorr_extract_REM: load filtered data, sleep scoring, and phasic/tonic scoring data, visualize REM sleep EEG using eeglab, visually identify and interpolate bad channels, visually mark segments that contain artefacts, save _goodREM.fdt and .set files that only contain 'good' (artefact-free) REM samples, save _goodREM.mat file containing information about sleep and phasic/tonic scoring, bad channels, and artefacts

5. D_ICA: perform ICA on _goodREM data using eeglab using AMICA1.6 plugin

6. E_manual_ICA: visualize independent components using eeglab, visually identify eye movement, cardiac, muscle activity, and channel noise components and remove those components, save as _manual_ICA.fdt and .set files

7. F_insert_REM_after_ICA_rereferencing: put back cleaned REM sleep data (after ICA) into the whole-night data and re-reference to the average across all channels

8. G_goodtrigs: save triggers that occurred during a 'good' REM sample in _nm_good.mat


wake: preprocessing scripts for wake data, run in order A-G

1. A_load_filter_data_wake: load and filter wake EEG data from .eeg file (raw data, provided), save as _fil.fdt and .set file

2. B_extract_triggers_wake: extract AEP triggers from .eeg file (raw data, provided), save in _nm.mat file

3. C_artcorr_wake: load filtered data, visualize EEG using eeglab, visually identify and interpolate bad channels, visually mark segments that contain artefacts, save _goodwake.fdt and .set files that only contain 'good' (artefact-free) wake samples, save _goodwake.mat file containing information about bad channels, and artefacts

4. D_ICA_wake: perform ICA on _goodwake data using eeglab using AMICA1.6 plugin

5. E_manual_ICA: visualize independent components using eeglab, visually identify eye movement, cardiac, muscle activity, and channel noise components and remove those components, save as _manual_ICA.fdt and .set files

6. F_insert_wake_after_ICA_rereferencing: re-reference wake data to the average across all channels

7. G_goodtrigs: save triggers that occurred during a 'good' wake sample in _nm_good.mat













