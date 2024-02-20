clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Matlab'));

%% adapt paths and recording start

Folderpath = '/vol/research/nemo/datasets/DRI007/data/hdEEG/DRI007_025/edf/alpha1';
filename = 'alpha_Vol101_Phase0.edf';

recording_start = [2023 09 05 11 12 00.000]; % [2022 01 20 22 05 37.000] adapt recording start date and time [yyyy mm dd hh mm ss.msmsms]

%% Get header

% [hdr data] = edfread2([Folderpath,filesep,filename]);
% hdr.startdate
% hdr.starttime
% year = str2num(['20',hdr.startdate(end-1:end)]);
% clear data
% [2022 01 20 22 05 37.000]

%% Filter and resample data

fs = 256;

% EEG = pop_fileio([Folderpath,filesep,filename],'channels',[1]);
EEG = pop_fileio([Folderpath,filesep,filename]);
EEG = pop_eegfiltnew(EEG,1,40);
EEG = pop_resample(EEG,fs);

% Take the channels out that we need
F4 = 29;
C4 = 24;
P4 = 19;
O2 = 18;
F3 = 3;
C3 = 8;
P3 = 14;
O1 = 16;
A1 = 10; % TP9
A2 = 21; % TP10
LOC = 130;
ROC = 131;
EMG = 129;

F4A1 = EEG.data(F4,:)-EEG.data(A1,:); % orig
C4A1 = EEG.data(C4,:)-EEG.data(A1,:); % orig
P4A1 = EEG.data(P4,:)-EEG.data(A1,:); % orig
O2A1 = EEG.data(O2,:)-EEG.data(A1,:); % orig
ROCA1 = EEG.data(ROC,:); %EEG.data(ROC,:)-EEG.data(A1,:); % orig

F3A2 = EEG.data(F3,:)-EEG.data(A2,:); % orig
C3A2 = EEG.data(C3,:)-EEG.data(A2,:); % orig
P3A2 = EEG.data(P3,:)-EEG.data(A2,:); % orig
O1A2 = EEG.data(O1,:)-EEG.data(A2,:); % orig
LOCA2 = EEG.data(LOC,:); %EEG.data(LOC,:)-EEG.data(A2,:); % orig

EMG2EMG3 = EEG.data(EMG,:); % orig

%%%%%
% F4A1 = EEG.data(F4,:)-EEG.data(A2,:); % re-reference right to RM
% C4A1 = EEG.data(C4,:)-EEG.data(A2,:); % re-reference right to RM
% P4A1 = EEG.data(P4,:)-EEG.data(A2,:); % re-reference right to RM
% O2A1 = EEG.data(O2,:)-EEG.data(A2,:); % re-reference right to RM
% ROCA1 = EEG.data(ROC,:); %EEG.data(ROC,:)-EEG.data(A2,:);

% F3A2 = EEG.data(F3,:)-EEG.data(A1,:); % re-reference left to LM
% C3A2 = EEG.data(C3,:)-EEG.data(A1,:); % re-reference left to LM
% P3A2 = EEG.data(P3,:)-EEG.data(A1,:); % re-reference left to LM
% O1A2 = EEG.data(O1,:)-EEG.data(A1,:); % re-reference left to LM
% LOCA2 = EEG.data(LOC,:); %EEG.data(LOC,:)-EEG.data(A1,:); 


data2 = [F4A1;C4A1;P4A1;O2A1;LOCA2;ROCA1;F3A2;C3A2;P3A2;O1A2;EMG2EMG3];
data2 = double(data2);

% labels = hdr.label;
labels={'F4A1' 'C4A1' 'P4A1' 'O2A1' 'LOCA2' 'ROCA1' 'F3A2' 'C3A2' 'P3A2' 'O1A2'  'EMG1EMG2'}; % labels need to be as follows to being able to open it in EOG Scoring App: {'F4A1' 'C4A1' 'P4A1' 'O2A1' 'LOCA2' 'ROCA1' 'F3A2' 'C3A2' 'P3A2' 'O1A2'  'EMG1EMG2'}

% pop_writeeeg(EEG, filename, 'TYPE', 'EDF')
addpath(genpath('/user/HS301/m17462/matlab/Biosig3.7.9'));
eeglab % open eeglab in order to run writeeeg function, otherwise it doesn't work (path error?)
writeeeg([Folderpath,filesep,filename(1:end-4),'_edf2.edf'], data2, fs, 'TYPE', 'EDF','Label',labels,'T0',recording_start); 
rmpath(genpath('/user/HS301/m17462/matlab/Biosig3.7.9'));


% edfw = edfwrite([Folderpath,filesep,'edftest3.edf'],hdr,data2); 


%%

% Folderpath = 'C:\Users\m17462\OneDrive - University of Surrey\Airforce_Surrey_psa\REM_Scoring\Valeria\af0468\AFOSR_AF0468_SEBN2';
% filename = 'AFOSR_AF0468_SEBN2.edf';
% 
% Folderpath = 'C:\Users\m17462\OneDrive - University of Surrey\Data\REM_study\pilots\night4\hdEEG\sara_test';
% filename = 'Sleep2_edftest';



% edfw = edfwrite([Folderpath,filesep,'edftest.edf'],hdr,data); 

% [hdr data] = edfread2('C:\Users\m17462\OneDrive - University of Surrey\Data\REM_study\pilots\night4\hdEEG\sara_test\edftest\edftest2.edf');


