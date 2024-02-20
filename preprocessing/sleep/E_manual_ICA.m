close all
clear 
clc

% addpath('/users/psychology01/software/automaticanalysis/');

% clear;
% aa_ver5

addpath(genpath('/user/HS301/m17462/matlab/eeglab/'));

%%

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/RSN_006/ICA2/';
Folderpath_dir = dir([Folderpath,'*ICA.set']);
% Folderpath_dir = dir([Folderpath,'*ICs_removed.set']);
filename = Folderpath_dir(1).name;

Folderpath_auxch = '/vol/research/nemo/datasets/RSN/data/hdEEG/RSN_006/';
Folderpath_auxch_dir = dir([Folderpath_auxch,'*sleep*_auxch_all.set']);

%%

% cd('/mnt/beegfs/users/psychology01/Henry')
% aap = aarecipe('aap_tasklist_meeg.xml');
% SPM = aas_inittoolbox(aap,'spm');
% SPM.load;

% EL = aas_inittoolbox(aap,'eeglab');
% EL.load;


%% 

% cond_name = {'sham';'sham2';'sham3';'stim';'stim2';'stim3'};

%%

close all

% sub = 12;
% cond = 1;


%% View IC's

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',[filename],'filepath',[Folderpath]);
% EEG = pop_loadset('filename',['IV_',dataName,'ICA.set'],'filepath',[Folderpath,filesep,dataName]);

[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 0, 1, 1);

pop_viewprops(EEG, 0, 1:30, {}, {}, 0);%,1,[],1:3);%,{},{},0,{},{});

EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); % get component data

% EEG = pop_saveset(EEG, 'filename', [filename(1:end-4),'_autoartcorr_reref_components'], 'filepath', Folderpath);

%%

EEG.manualrejcomp = [11 18 17 21 29]; % define eye movement components

%% Plot eye movements and eye components

eeg_aux = pop_loadset('filename', [Folderpath_auxch_dir(1).name], 'filepath', Folderpath_auxch);

matfile = dir([Folderpath_auxch,'*czref_goodREM.mat']);
load([Folderpath_auxch,matfile(1).name],'rem_goodsamp2');
eeg_aux.data = eeg_aux.data(:,rem_goodsamp2);

if ~isequal(length(rem_goodsamp2),length(EEG.data))
    error('Auxch file and EEG data file do not have same length');
end

EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:); % get component data

c1 = EEG.icaact(EEG.manualrejcomp(1),:);
c2 = EEG.icaact(EEG.manualrejcomp(2),:);

lEOG = eeg_aux.data(3,:);
rEOG = eeg_aux.data(4,:);
% lEOG = eeg_aux.data(4,:);
% rEOG = eeg_aux.data(5,:);

epochl = 10;

for ep = 30:40

figure

subplot(2,1,1)
plot(lEOG((ep-1)*epochl*EEG.srate+1:ep*epochl*EEG.srate),'Color',[1 0 0])
hold on
plot(rEOG((ep-1)*epochl*EEG.srate+1:ep*epochl*EEG.srate),'Color',[0 0.4470 0.7410])
ylim([-100 100])
title('EOG')
legend({'lEOG' 'rEOG'})
    
subplot(2,1,2)
plot(c1((ep-1)*epochl*EEG.srate+1:ep*epochl*EEG.srate),'Color',[0.3010 0.7450 0.9330])
hold on
plot(c2((ep-1)*epochl*EEG.srate+1:ep*epochl*EEG.srate),'Color',[0.8500 0.3250 0.0980])
ylim([-20 20])
title('Components')
legend({['IC',num2str(EEG.manualrejcomp(1))] ['IC',num2str(EEG.manualrejcomp(2))]})

end

    %% Reject IC's

% close all
EEG = eeg_checkset(EEG);
EEG = pop_subcomp(EEG, EEG.manualrejcomp, 1); % define manually which components to reject

% pop_eegplot(EEG);

pop_newset(ALLEEG, EEG, 1,'savenew',[Folderpath,filename(1:end-4),'_manual_ICA'],'gui','off'); 
