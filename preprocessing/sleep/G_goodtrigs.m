clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('//user/HS301/m17462/matlab/Scripts/RSN'));

addpath /mnt/beegfs/users/psychology01/useful_functions

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/RSN_024/';
goodrem_mat_file = dir([Folderpath,'*czref_goodREM.mat']);
nm_mat_file = dir([Folderpath,'*sleep*nm.mat']);

%% Load scoring data

load([Folderpath,goodrem_mat_file(1).name]);

%% Load nm file

load([Folderpath,nm_mat_file(1).name]);

%% Exclude bad ERP triggers

[nm.trigs_ECHT_ERP_good nm.trigs_ECHT_ERP_goodndx] = intersect(nm.trigs_ECHT_ERP,rem_goodsamp2);
[nm.trigs_StimTrak_ERP_good nm.trigs_StimTrak_ERP_goodndx] = intersect(nm.trigs_StimTrak_ERP,rem_goodsamp2);

%% Exclude bad windows and triggers

for con = 1:8
        
    [nm.OFF_start_good{con} OFF_start_goodndx] = intersect(nm.OFF_start{con},rem_goodsamp2);
    nm.OFF_end_good{con} = nm.OFF_end{con}(OFF_start_goodndx);
    
    [nm.ON_start_good{con} ON_start_goodndx] = intersect(nm.ON_start{con},rem_goodsamp2);
    nm.ON_end_good{con} = nm.ON_end{con}(ON_start_goodndx);
    
    [nm.ON_trigs_ECHT_good{con} nm.ON_trigs_ECHT_goodndx] = intersect(nm.ON_trigs_ECHT{con},rem_goodsamp2);
%   nm.vol_trigs_ECHT_good{con} =    nm.vol_trigs_ECHT{con}(nm.ON_trigs_ECHT_goodndx); 
%   voltrigs are wrong atm, need to correct script B and rerun

    [nm.ON_trigs_StimTrak_good{con} nm.ON_trigs_StimTrak_goodndx] = intersect(nm.ON_trigs_StimTrak{con},rem_goodsamp2);
%   nm.vol_trigs_StimTrak_good{con} =    nm.vol_trigs_StimTrak{con}(nm.ON_trigs_StimTrak_goodndx); 
%   voltrigs are wrong atm, need to correct script B and rerun    
   
end

%%

save([Folderpath,nm_mat_file(1).name(1:end-4),'_good.mat'],'nm');




