clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));


Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/RSN_024/';

czref_file = dir([Folderpath,'*sleep*czref.set']);
manual_ICA_file = dir([Folderpath,'ICA2/','*czref_goodREM_ICA_manual_ICA.set']);
goodrem_mat_file = dir([Folderpath,'*czref_goodREM.mat']);

%%

EEG = pop_loadset('filename',[czref_file(1).name],'filepath',[Folderpath]);
EEG = pop_select(EEG,'nochannel',{'Echt' 'EMG' 'lEOG' 'rEOG' 'ECG'});

EEG_mICA = pop_loadset('filename',[manual_ICA_file(1).name],'filepath',[Folderpath,'ICA2']);

load([Folderpath,goodrem_mat_file(1).name]);

EEG.data(:,rem_goodsamp2)= EEG_mICA.data;

EEG = pop_saveset(EEG, 'filename', [czref_file(1).name(1:end-4),'_mICA'], 'filepath', Folderpath);


% put back Cz
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:) = zeros(1, EEG.pnts);
EEG.chanlocs(1,EEG.nbchan).labels = 'Cz';
EEG = pop_chanedit(EEG, 'lookup','/user/HS301/m17462/matlab/Scripts/RSN/preprocessing/sleep/standard_1005.elc','eval','chans = pop_chancenter( chans, [],[]);');

originaleeg = EEG;

%% Average referencing

EEG = pop_reref(EEG, []); % average referencing
EEG = pop_saveset(EEG, 'filename', [czref_file(1).name(1:end-4),'_mICA_avref'], 'filepath', Folderpath);
clear EEG

