clear all;
close all;

addpath(genpath('/user/HS301/m17462/matlab/eeglab'));
addpath(genpath('/user/HS301/m17462/matlab/Scripts/RSN'));
addpath(genpath('/user/HS301/m17462/matlab/Henry/useful_functions'));
addpath(genpath('/user/HS301/m17462/matlab/colorGradient'));

Savefolder = '/vol/research/nemo/datasets/RSN/data/analysis/Figures/';

incl_sub = setdiff(1:19,12);

Folderpath = '/vol/research/nemo/datasets/RSN/data/hdEEG/';
sub_Folderpath = dir([Folderpath,'RSN*']);

%% load all trial data - REM

for s = 1:length(sub_Folderpath)
    
    erp_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_sleep*_erp_good.mat']);
    
    load([Folderpath,sub_Folderpath(s).name,filesep,erp_good_file(1).name]);
    
    ERP_trials_allsub_REM{s} = ERP; 
    
end


%% load all trial data - wake eve

for s = 1:length(sub_Folderpath)
    
    erp_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_wake_e*_erp_good.mat']);
    
    load([Folderpath,sub_Folderpath(s).name,filesep,erp_good_file(1).name]);
    
    ERP_trials_allsub_wake_e{s} = ERP; 
    
end

%% load all trial data - wake mor

for s = 1:length(sub_Folderpath)
    
    erp_good_file = dir([Folderpath,sub_Folderpath(s).name,filesep,'*_wake_m*_erp_good.mat']);
    
    load([Folderpath,sub_Folderpath(s).name,filesep,erp_good_file(1).name]);
    
    ERP_trials_allsub_wake_m{s} = ERP; 
    
end

%% phase response curve - tonic REM

colors = linspecer(4);
             
%compute for many onset phases speeding up/slowing down
      
state = 1; % 1 = eyes open, 2 = eyes closed
chan = 2;
             
clear angle_start_end phase_after R_after angs clear sem_end pred_angle R_end

trigsamp = 501;
nbins = 10;
fs = 500;
osc_freq = 9;

% angs(1,:) = wrapToPi(deg2rad([0:360/10:360]));
angs(1,:) = wrapToPi(deg2rad([0:360/nbins:360]));
angs(2,:) = angs-pi/nbins;
angs(3,:) = angs(1,:)+pi/nbins;
             
fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 .4 1])
hold on
% for chan = [12]
%                  clear phase_after phase_start_end
G = 1;

cycs_all = [20/2 100/2 150/2 200/2 250/2];

for angle = 1:length(cycs_all) % [10 75 85 95 105];%75 100 125 150]%(25:25:200)
    cycs = cycs_all(angle);
    A = 1;
    for ang = 1:size(angs,2) %0:18:360;%[36:10:360-36]
             
         for cond = 1
                 clear y
                 for s = 1:length(incl_sub)
                     
                     sub = incl_sub(s);
                     
%                         tmp = [];
%                         for condo = 1:4
%                          tmp = [tmp;group_phase.sub{sub}.chan{chan}.cond{condo}];
                        tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.tonic_ndx,:); % trials x samps
%                         tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.goodtrial_ndx,:); % trials x samps
%                           tmp = ERP_trials_allsub_wake_m{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_m{sub}.good_trial_ndx{state},:);
%                           tmp = ERP_trials_allsub_wake_e{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_e{sub}.good_trial_ndx{state},:);
%                         end
                                              
                        tmp = tmp(tmp(:,trigsamp)>angs(2,ang)&tmp(:,trigsamp)<angs(3,ang),:);
                         
                        if ~isempty(tmp)
                                      
                        Mn = circ_mean(tmp); % average across all trials in that bin
                        R = circ_r(tmp);
                        to_find = Mn(trigsamp);
                     
                        cosed = cos(circ_dist(Mn,to_find));
                                          
                        [pks,locs] = findpeaks(cosed);
                                       
                     
                        phase_before(s,cond) =  Mn(trigsamp);
                        cycle_length(s,cond) = abs(trigsamp-locs(find(locs==trigsamp)-1)); % duration from previous peak to stim peak in samps
                        phase_after(s,cond,1) = Mn(trigsamp);
                        R_after(s,cond,1) = R(trigsamp);
                                        
                    
                        D = (diff(locs(locs>=trigsamp-100 & locs<=trigsamp)));
                        %D = D(D>30);
                        cycles(s,cond,chan) = mean(D); % mean duration of all cycles between -200 ms and 0
                    
                    
                            for cyc = 1:round(cycs)
                                   
                                phase_after(s,cond,cyc+1) = Mn(trigsamp+cyc); % mean phase for all samples from trigsamp to endsamp
                                phase_change(s,cond,cyc) = circ_dist(phase_after(s,cond,cyc),phase_before(s,cond)); % phase distance from each samp after to trigsamp
                                R_after(s,cond,cyc+1) = R(trigsamp+cyc);%  resultant for all samples from trigsamp to endsamp      
                        
                            end
                    
                                                               
                        else
                 
                            phase_after(s,cond,:) = NaN(1,1,cycs+1);
                            R_after(s,cond,:) = NaN(1,1,cycs+1);
                        end
            
                 end
           
         end
             
             
%              if cycs == trigsamp
%              phase_series(ang,:) = squeeze(circ_mean(phase_after));
%              end
%              clear cumul_phase_change
             
             cumul_phase_change = unwrap(squeeze(phase_after(:,1,:))')';

             angle_start_end(1,ang) = (wrapToPi(angs(1,ang))); % angle at center of bin
             tmp = phase_after(:,cond,end); % phase of endsamp
             tmp2 = R_after(:,cond,end); %  resultant of endsamp
             tmp2 = tmp2(~isnan(tmp2));
             tmp = tmp(~isnan(tmp));
                 
             angle_start_end(2,ang) = (wrapToPi(circ_mean(tmp))); % mean of endsamp phase across participants
             yy = num2str(cycs);
             yy = str2num(yy(end-1:end))/(fs/osc_freq); % calculate how much of a cycle this amount of time is
             pred_angle(ang) = wrapToPi(angs(1,ang)+ 2*yy*pi);
             R_end(1,ang) = circ_r(tmp,tmp2); % resultant of endsamp across participants
                
             angle_start_end(3,ang) = (circ_dist(circ_mean(tmp),pred_angle(ang))); % phase distance from endsamp to predicted angle

             sem_end(1,ang) = circ_std(tmp)./sqrt(length(incl_sub)); % sem of end phase
             sem_end(2,ang) = circ_std(circ_dist(angs(1,ang),tmp))./sqrt(length(incl_sub)); % sem of difference between angle at center of bin and end phase

             A = A+1;
             
             end
             
             [B I] = sort((angle_start_end(2,:)));
             
             colorz = linspecer(length(angle_start_end(3,:)));

             subplot(5,4,(angle-1)*4+1)
             hold on
             errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color','k','LineStyle','none') % angle at center of bin vs mean phase at end, errorbar
             scatter(angle_start_end(1,:),pred_angle,2.5,'k','filled') % angle at center of bin vs predicted angle
             
             for i = 1:length(angle_start_end(3,:))

                if chan == 2

                    scatter(angle_start_end(1,i),angle_start_end(2,i),25,colorz(i,:),'filled','MarkerEdgeColor','k') % angle at center of bins vs mean at end 
                else
                    scatter(angle_start_end(1,:),angle_start_end(2,:),25,colors(2,:))
                    errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color',colors(2,:),'LineStyle','none') 
                end
             end
             
%              G = G+1;

             ylim([-pi pi])
             xlim([-pi pi])
             title(['Delay: ',num2str(cycs/fs*1000),' ms'])     

             xlabel('Start Phase')
             ylabel('End Phase')
             
             subplot(5,4,(angle-1)*4+2)
             hold on
             errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color','k','LineStyle','none') % angle at center of bin vs phase distance from endsamp to predicted angle, errorbar
             
             for i = 1:length(angle_start_end(3,:))
                if chan == 2
                    scatter(angle_start_end(1,i),angle_start_end(3,i),25,colorz(i,:),'filled','MarkerEdgeColor','k')                         
                else
                    scatter(angle_start_end(1,:),angle_start_end(3,:),25,colors(2,:))
                    errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color',colors(2,:),'LineStyle','none')    
                end
             end
             
%             G = G+1;
            a = 1;
            for A = [0 90 180 270]
%             for A = [340 60 160 240]
                xline(wrapToPi(deg2rad(A)),'Color',colors(a,:)) 
                a = a+1;
            end
             %xlim([-180 180])
             xlabel('Start Phase')
             ylabel('Phase Change')
             yline(0)
%              text(.95*-pi,-.9*pi,'slower')
             %ylim([-200 200])
             %ylim([.1 .14])
             ylim([-pi pi])
             xlim([-pi pi])
%              text(.95*-pi,.9*pi,'faster')
                         
                          
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(2,:)')],[0 (circ_r(angle_start_end(2,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%             if chan == 2
%                 B = 0.5;
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(2,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     B = B+.025;
%                     hold on
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%   
%                         hold on
%                     end
%                 end
%                 
%             else
%                 polarplot(angle_start_end(1,:),(circ_dist(circ_mean(angle_start_end(3,:)'),angle_start_end(3,:)))+pi,'Color',colors(2,:))
% 
%             end
%              hold on
% 
%              rlim([0 1.1])
%              G = G+1;
%              
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(3,:)')],[0 (circ_r(angle_start_end(3,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%              if chan == 2
%                 B = .5;
%                 
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(3,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     hold on
%                     B = B+.025;
%     
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%                         hold on
%                     end
%                 end
%              end
%              
%              hold on
%              rlim([0 1.1])
%              G = G+1;
             
end
             
             
saveas(fig,[Savefolder,'Figure6_phase_response_curve_freq',num2str(osc_freq),'_alphafilt_tonic_REM.svg']);

%% phase response curve - phasic REM

% colors = linspecer(4);
%              
% %compute for many onset phases speeding up/slowing down
%       
% state = 1; % 1 = eyes open, 2 = eyes closed
% chan = 2;
%              
% clear angle_start_end phase_after R_after angs clear sem_end pred_angle R_end
% 
% trigsamp = 501;
% nbins = 10;
% fs = 500;
% osc_freq = 9;
% 
% % angs(1,:) = wrapToPi(deg2rad([0:360/10:360]));
% angs(1,:) = wrapToPi(deg2rad([0:360/nbins:360]));
% angs(2,:) = angs-pi/nbins;
% angs(3,:) = angs(1,:)+pi/nbins;
%              
% fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 .4 1])
% hold on
% % for chan = [12]
% %                  clear phase_after phase_start_end
% G = 1;
% 
% for cycs = [20/2 100/2 150/2 200/2 250/2]% [10 75 85 95 105];%75 100 125 150]%(25:25:200)
%     A = 1;
%     for ang = 1:size(angs,2) %0:18:360;%[36:10:360-36]
%              
%          for cond = 1
%                  clear y
%                  for s = 1:length(incl_sub)
%                      
%                      sub = incl_sub(s);
%                      
% %                         tmp = [];
% %                         for condo = 1:4
% %                          tmp = [tmp;group_phase.sub{sub}.chan{chan}.cond{condo}];
%                         tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.phasic_ndx,:); % trials x samps
% %                         tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.goodtrial_ndx,:); % trials x samps
% %                           tmp = ERP_trials_allsub_wake_m{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_m{sub}.good_trial_ndx{state},:);
% %                           tmp = ERP_trials_allsub_wake_e{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_e{sub}.good_trial_ndx{state},:);
% %                         end
%                                               
%                         tmp = tmp(tmp(:,trigsamp)>angs(2,ang)&tmp(:,trigsamp)<angs(3,ang),:);
%                          
%                         if ~isempty(tmp)
%                                       
%                         Mn = circ_mean(tmp); % average across all trials in that bin
%                         R = circ_r(tmp);
%                         to_find = Mn(trigsamp);
%                      
%                         cosed = cos(circ_dist(Mn,to_find));
%                                           
%                         [pks,locs] = findpeaks(cosed);
%                                        
%                      
%                         phase_before(s,cond) =  Mn(trigsamp);
%                         cycle_length(s,cond) = abs(trigsamp-locs(find(locs==trigsamp)-1)); % duration from previous peak to stim peak in samps
%                         phase_after(s,cond,1) = Mn(trigsamp);
%                         R_after(s,cond,1) = R(trigsamp);
%                                         
%                     
%                         D = (diff(locs(locs>=trigsamp-100 & locs<=trigsamp)));
%                         %D = D(D>30);
%                         cycles(s,cond,chan) = mean(D); % mean duration of all cycles between -200 ms and 0
%                     
%                     
%                             for cyc = 1:round(cycs)
%                                    
%                                 phase_after(s,cond,cyc+1) = Mn(trigsamp+cyc); % mean phase for all samples from trigsamp to endsamp
%                                 phase_change(s,cond,cyc) = circ_dist(phase_after(s,cond,cyc),phase_before(s,cond)); % phase distance from each samp after to trigsamp
%                                 R_after(s,cond,cyc+1) = R(trigsamp+cyc);%  resultant for all samples from trigsamp to endsamp      
%                         
%                             end
%                     
%                                                                
%                         else
%                  
%                             phase_after(s,cond,:) = NaN(1,1,cycs+1);
%                             R_after(s,cond,:) = NaN(1,1,cycs+1);
%                         end
%             
%                  end
%            
%          end
%              
%              
% %              if cycs == trigsamp
% %              phase_series(ang,:) = squeeze(circ_mean(phase_after));
% %              end
% %              clear cumul_phase_change
%              
%              cumul_phase_change = unwrap(squeeze(phase_after(:,1,:))')';
% 
%              angle_start_end(1,ang) = (wrapToPi(angs(1,ang))); % angle at center of bin
%              tmp = phase_after(:,cond,end); % phase of endsamp
%              tmp2 = R_after(:,cond,end); %  resultant of endsamp
%              tmp2 = tmp2(~isnan(tmp2));
%              tmp = tmp(~isnan(tmp));
%                  
%              angle_start_end(2,ang) = (wrapToPi(circ_mean(tmp))); % mean of endsamp phase across participants
%              yy = num2str(cycs);
%              yy = str2num(yy(end-1:end))/(fs/osc_freq); % calculate how much of a cycle this amount of time is
%              pred_angle(ang) = wrapToPi(angs(1,ang)+ 2*yy*pi);
%              R_end(1,ang) = circ_r(tmp,tmp2); % resultant of endsamp across participants
%                 
%              angle_start_end(3,ang) = (circ_dist(circ_mean(tmp),pred_angle(ang))); % phase distance from endsamp to predicted angle
% 
%              sem_end(1,ang) = circ_std(tmp)./sqrt(length(incl_sub)); % sem of end phase
%              sem_end(2,ang) = circ_std(circ_dist(angs(1,ang),tmp))./sqrt(length(incl_sub)); % sem of difference between angle at center of bin and end phase
% 
%              A = A+1;
%              
%              end
%              
%              [B I] = sort((angle_start_end(2,:)));
%              
%              colorz = linspecer(length(angle_start_end(3,:)));
% 
%              subplot(5,4,G)
%              hold on
%              errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color','k','LineStyle','none') % angle at center of bin vs mean phase at end, errorbar
%              scatter(angle_start_end(1,:),pred_angle,2.5,'k','filled') % angle at center of bin vs predicted angle
%              
%              for i = 1:length(angle_start_end(3,:))
% 
%                 if chan == 2
% 
%                     scatter(angle_start_end(1,i),angle_start_end(2,i),25,colorz(i,:),'filled','MarkerEdgeColor','k') % angle at center of bins vs mean at end 
%                 else
%                     scatter(angle_start_end(1,:),angle_start_end(2,:),25,colors(2,:))
%                     errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color',colors(2,:),'LineStyle','none') 
%                 end
%              end
%              
%              G = G+1;
% 
%              ylim([-pi pi])
%              xlim([-pi pi])
%              title(['Delay: ',num2str(cycs/fs*1000),' ms'])     
% 
%              xlabel('Start Phase')
%              ylabel('End Phase')
%              
%              subplot(5,4,G)
%              hold on
%              errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color','k','LineStyle','none') % angle at center of bin vs phase distance from endsamp to predicted angle, errorbar
%              
%              for i = 1:length(angle_start_end(3,:))
%                 if chan == 2
%                     scatter(angle_start_end(1,i),angle_start_end(3,i),25,colorz(i,:),'filled','MarkerEdgeColor','k')                         
%                 else
%                     scatter(angle_start_end(1,:),angle_start_end(3,:),25,colors(2,:))
%                     errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color',colors(2,:),'LineStyle','none')    
%                 end
%              end
%              
%             G = G+1;
%             a = 1;
%             for A = [0 90 180 270]
% %             for A = [340 60 160 240]
%                 xline(wrapToPi(deg2rad(A)),'Color',colors(a,:)) 
%                 a = a+1;
%             end
%              %xlim([-180 180])
%              xlabel('Start Phase')
%              ylabel('Phase Change')
%              yline(0)
%              text(.95*-pi,-.9*pi,'slower')
%              %ylim([-200 200])
%              %ylim([.1 .14])
%              ylim([-pi pi])
%              xlim([-pi pi])
%              text(.95*-pi,.9*pi,'faster')
%                          
%                           
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(2,:)')],[0 (circ_r(angle_start_end(2,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%             if chan == 2
%                 B = 0.5;
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(2,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     B = B+.025;
%                     hold on
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%   
%                         hold on
%                     end
%                 end
%                 
%             else
%                 polarplot(angle_start_end(1,:),(circ_dist(circ_mean(angle_start_end(3,:)'),angle_start_end(3,:)))+pi,'Color',colors(2,:))
% 
%             end
%              hold on
% 
%              rlim([0 1.1])
%              G = G+1;
%              
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(3,:)')],[0 (circ_r(angle_start_end(3,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%              if chan == 2
%                 B = .5;
%                 
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(3,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     hold on
%                     B = B+.025;
%     
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%                         hold on
%                     end
%                 end
%              end
%              
%              hold on
%              rlim([0 1.1])
%              G = G+1;
%              
% end
%              
             
% saveas(fig,[Savefolder,'Figure6_phase_response_curve_freq',num2str(osc_freq),'_phasic_REM.svg']);

%% phase response curve - wake eve (eyes closed)

% colors = linspecer(4);
%              
% %compute for many onset phases speeding up/slowing down
%       
% state = 2; % 1 = eyes open, 2 = eyes closed
% chan = 2;
%              
% clear angle_start_end phase_after R_after angs clear sem_end pred_angle R_end
% 
% trigsamp = 501;
% nbins = 10;
% fs = 500;
% osc_freq = 10;
% 
% % angs(1,:) = wrapToPi(deg2rad([0:360/10:360]));
% angs(1,:) = wrapToPi(deg2rad([0:360/nbins:360]));
% angs(2,:) = angs-pi/nbins;
% angs(3,:) = angs(1,:)+pi/nbins;
%              
% fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 .4 1])
% hold on
% % for chan = [12]
% %                  clear phase_after phase_start_end
% G = 1;
% 
% for cycs = [20/2 100/2 150/2 200/2 250/2]% [10 75 85 95 105];%75 100 125 150]%(25:25:200)
%     A = 1;
%     for ang = 1:size(angs,2) %0:18:360;%[36:10:360-36]
%              
%          for cond = 1
%                  clear y
%                  for s = 1:length(incl_sub)
%                      
%                      sub = incl_sub(s);
%                      
% %                         tmp = [];
% %                         for condo = 1:4
% %                          tmp = [tmp;group_phase.sub{sub}.chan{chan}.cond{condo}];
% %                         tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.phasic_ndx,:); % trials x samps
% %                         tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.goodtrial_ndx,:); % trials x samps
% %                           tmp = ERP_trials_allsub_wake_m{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_m{sub}.good_trial_ndx{state},:);
%                           tmp = ERP_trials_allsub_wake_e{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_e{sub}.good_trial_ndx{state},:);
% %                         end
%                                               
%                         tmp = tmp(tmp(:,trigsamp)>angs(2,ang)&tmp(:,trigsamp)<angs(3,ang),:);
%                          
%                         if ~isempty(tmp)
%                                       
%                         Mn = circ_mean(tmp); % average across all trials in that bin
%                         R = circ_r(tmp);
%                         to_find = Mn(trigsamp);
%                      
%                         cosed = cos(circ_dist(Mn,to_find));
%                                           
%                         [pks,locs] = findpeaks(cosed);
%                                        
%                      
%                         phase_before(s,cond) =  Mn(trigsamp);
%                         cycle_length(s,cond) = abs(trigsamp-locs(find(locs==trigsamp)-1)); % duration from previous peak to stim peak in samps
%                         phase_after(s,cond,1) = Mn(trigsamp);
%                         R_after(s,cond,1) = R(trigsamp);
%                                         
%                     
%                         D = (diff(locs(locs>=trigsamp-100 & locs<=trigsamp)));
%                         %D = D(D>30);
%                         cycles(s,cond,chan) = mean(D); % mean duration of all cycles between -200 ms and 0
%                     
%                     
%                             for cyc = 1:round(cycs)
%                                    
%                                 phase_after(s,cond,cyc+1) = Mn(trigsamp+cyc); % mean phase for all samples from trigsamp to endsamp
%                                 phase_change(s,cond,cyc) = circ_dist(phase_after(s,cond,cyc),phase_before(s,cond)); % phase distance from each samp after to trigsamp
%                                 R_after(s,cond,cyc+1) = R(trigsamp+cyc);%  resultant for all samples from trigsamp to endsamp      
%                         
%                             end
%                     
%                                                                
%                         else
%                  
%                             phase_after(s,cond,:) = NaN(1,1,cycs+1);
%                             R_after(s,cond,:) = NaN(1,1,cycs+1);
%                         end
%             
%                  end
%            
%          end
%              
%              
% %              if cycs == trigsamp
% %              phase_series(ang,:) = squeeze(circ_mean(phase_after));
% %              end
% %              clear cumul_phase_change
%              
%              cumul_phase_change = unwrap(squeeze(phase_after(:,1,:))')';
% 
%              angle_start_end(1,ang) = (wrapToPi(angs(1,ang))); % angle at center of bin
%              tmp = phase_after(:,cond,end); % phase of endsamp
%              tmp2 = R_after(:,cond,end); %  resultant of endsamp
%              tmp2 = tmp2(~isnan(tmp2));
%              tmp = tmp(~isnan(tmp));
%                  
%              angle_start_end(2,ang) = (wrapToPi(circ_mean(tmp))); % mean of endsamp phase across participants
%              yy = num2str(cycs);
%              yy = str2num(yy(end-1:end))/(fs/osc_freq); % calculate how much of a cycle this amount of time is
%              pred_angle(ang) = wrapToPi(angs(1,ang)+ 2*yy*pi);
%              R_end(1,ang) = circ_r(tmp,tmp2); % resultant of endsamp across participants
%                 
%              angle_start_end(3,ang) = (circ_dist(circ_mean(tmp),pred_angle(ang))); % phase distance from endsamp to predicted angle
% 
%              sem_end(1,ang) = circ_std(tmp)./sqrt(length(incl_sub)); % sem of end phase
%              sem_end(2,ang) = circ_std(circ_dist(angs(1,ang),tmp))./sqrt(length(incl_sub)); % sem of difference between angle at center of bin and end phase
% 
%              A = A+1;
%              
%              end
%              
%              [B I] = sort((angle_start_end(2,:)));
%              
%              colorz = linspecer(length(angle_start_end(3,:)));
% 
%              subplot(5,4,G)
%              hold on
%              errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color','k','LineStyle','none') % angle at center of bin vs mean phase at end, errorbar
%              scatter(angle_start_end(1,:),pred_angle,2.5,'k','filled') % angle at center of bin vs predicted angle
%              
%              for i = 1:length(angle_start_end(3,:))
% 
%                 if chan == 2
% 
%                     scatter(angle_start_end(1,i),angle_start_end(2,i),25,colorz(i,:),'filled','MarkerEdgeColor','k') % angle at center of bins vs mean at end 
%                 else
%                     scatter(angle_start_end(1,:),angle_start_end(2,:),25,colors(2,:))
%                     errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color',colors(2,:),'LineStyle','none') 
%                 end
%              end
%              
%              G = G+1;
% 
%              ylim([-pi pi])
%              xlim([-pi pi])
%              title(['Delay: ',num2str(cycs/fs*1000),' ms'])     
% 
%              xlabel('Start Phase')
%              ylabel('End Phase')
%              
%              subplot(5,4,G)
%              hold on
%              errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color','k','LineStyle','none') % angle at center of bin vs phase distance from endsamp to predicted angle, errorbar
%              
%              for i = 1:length(angle_start_end(3,:))
%                 if chan == 2
%                     scatter(angle_start_end(1,i),angle_start_end(3,i),25,colorz(i,:),'filled','MarkerEdgeColor','k')                         
%                 else
%                     scatter(angle_start_end(1,:),angle_start_end(3,:),25,colors(2,:))
%                     errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color',colors(2,:),'LineStyle','none')    
%                 end
%              end
%              
%             G = G+1;
%             a = 1;
%             for A = [0 90 180 270]
% %             for A = [340 60 160 240]
%                 xline(wrapToPi(deg2rad(A)),'Color',colors(a,:)) 
%                 a = a+1;
%             end
%              %xlim([-180 180])
%              xlabel('Start Phase')
%              ylabel('Phase Change')
%              yline(0)
%              text(.95*-pi,-.9*pi,'slower')
%              %ylim([-200 200])
%              %ylim([.1 .14])
%              ylim([-pi pi])
%              xlim([-pi pi])
%              text(.95*-pi,.9*pi,'faster')
%                          
%                           
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(2,:)')],[0 (circ_r(angle_start_end(2,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%             if chan == 2
%                 B = 0.5;
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(2,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     B = B+.025;
%                     hold on
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%   
%                         hold on
%                     end
%                 end
%                 
%             else
%                 polarplot(angle_start_end(1,:),(circ_dist(circ_mean(angle_start_end(3,:)'),angle_start_end(3,:)))+pi,'Color',colors(2,:))
% 
%             end
%              hold on
% 
%              rlim([0 1.1])
%              G = G+1;
%              
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(3,:)')],[0 (circ_r(angle_start_end(3,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%              if chan == 2
%                 B = .5;
%                 
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(3,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     hold on
%                     B = B+.025;
%     
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%                         hold on
%                     end
%                 end
%              end
%              
%              hold on
%              rlim([0 1.1])
%              G = G+1;
%              
% end
%              
%              
% % saveas(fig,[Savefolder,'Figure6_phase_response_curve_freq',num2str(osc_freq),'_wake_eve_EC.svg']);
% 
% %% phase response curve - wake mor (eyes closed)
% 
% colors = linspecer(4);
%              
% %compute for many onset phases speeding up/slowing down
%       
% state = 2; % 1 = eyes open, 2 = eyes closed
% chan = 2;
%              
% clear angle_start_end phase_after R_after angs clear sem_end pred_angle R_end
% 
% trigsamp = 501;
% nbins = 10;
% fs = 500;
% osc_freq = 10;
% 
% % angs(1,:) = wrapToPi(deg2rad([0:360/10:360]));
% angs(1,:) = wrapToPi(deg2rad([0:360/nbins:360]));
% angs(2,:) = angs-pi/nbins;
% angs(3,:) = angs(1,:)+pi/nbins;
%              
% fig = figure('Renderer','painters','units','normalized','outerposition',[0 0 .4 1])
% hold on
% % for chan = [12]
% %                  clear phase_after phase_start_end
% G = 1;
% 
% for cycs = [20/2 100/2 150/2 200/2 250/2]% [10 75 85 95 105];%75 100 125 150]%(25:25:200)
%     A = 1;
%     for ang = 1:size(angs,2) %0:18:360;%[36:10:360-36]
%              
%          for cond = 1
%                  clear y
%                  for s = 1:length(incl_sub)
%                      
%                      sub = incl_sub(s);
%                      
% %                         tmp = [];
% %                         for condo = 1:4
% %                          tmp = [tmp;group_phase.sub{sub}.chan{chan}.cond{condo}];
% %                         tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.phasic_ndx,:); % trials x samps
% %                         tmp = ERP_trials_allsub_REM{sub}.trial_data_alphafilt_phase(ERP_trials_allsub_REM{sub}.goodtrial_ndx,:); % trials x samps
%                           tmp = ERP_trials_allsub_wake_m{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_m{sub}.good_trial_ndx{state},:);
% %                           tmp = ERP_trials_allsub_wake_e{sub}.trial_data_alphafilt_phase{state}(ERP_trials_allsub_wake_e{sub}.good_trial_ndx{state},:);
% %                         end
%                                               
%                         tmp = tmp(tmp(:,trigsamp)>angs(2,ang)&tmp(:,trigsamp)<angs(3,ang),:);
%                          
%                         if ~isempty(tmp)
%                                       
%                         Mn = circ_mean(tmp); % average across all trials in that bin
%                         R = circ_r(tmp);
%                         to_find = Mn(trigsamp);
%                      
%                         cosed = cos(circ_dist(Mn,to_find));
%                                           
%                         [pks,locs] = findpeaks(cosed);
%                                        
%                      
%                         phase_before(s,cond) =  Mn(trigsamp);
%                         cycle_length(s,cond) = abs(trigsamp-locs(find(locs==trigsamp)-1)); % duration from previous peak to stim peak in samps
%                         phase_after(s,cond,1) = Mn(trigsamp);
%                         R_after(s,cond,1) = R(trigsamp);
%                                         
%                     
%                         D = (diff(locs(locs>=trigsamp-100 & locs<=trigsamp)));
%                         %D = D(D>30);
%                         cycles(s,cond,chan) = mean(D); % mean duration of all cycles between -200 ms and 0
%                     
%                     
%                             for cyc = 1:round(cycs)
%                                    
%                                 phase_after(s,cond,cyc+1) = Mn(trigsamp+cyc); % mean phase for all samples from trigsamp to endsamp
%                                 phase_change(s,cond,cyc) = circ_dist(phase_after(s,cond,cyc),phase_before(s,cond)); % phase distance from each samp after to trigsamp
%                                 R_after(s,cond,cyc+1) = R(trigsamp+cyc);%  resultant for all samples from trigsamp to endsamp      
%                         
%                             end
%                     
%                                                                
%                         else
%                  
%                             phase_after(s,cond,:) = NaN(1,1,cycs+1);
%                             R_after(s,cond,:) = NaN(1,1,cycs+1);
%                         end
%             
%                  end
%            
%          end
%              
%              
% %              if cycs == trigsamp
% %              phase_series(ang,:) = squeeze(circ_mean(phase_after));
% %              end
% %              clear cumul_phase_change
%              
%              cumul_phase_change = unwrap(squeeze(phase_after(:,1,:))')';
% 
%              angle_start_end(1,ang) = (wrapToPi(angs(1,ang))); % angle at center of bin
%              tmp = phase_after(:,cond,end); % phase of endsamp
%              tmp2 = R_after(:,cond,end); %  resultant of endsamp
%              tmp2 = tmp2(~isnan(tmp2));
%              tmp = tmp(~isnan(tmp));
%                  
%              angle_start_end(2,ang) = (wrapToPi(circ_mean(tmp))); % mean of endsamp phase across participants
%              yy = num2str(cycs);
%              yy = str2num(yy(end-1:end))/(fs/osc_freq); % calculate how much of a cycle this amount of time is
%              pred_angle(ang) = wrapToPi(angs(1,ang)+ 2*yy*pi);
%              R_end(1,ang) = circ_r(tmp,tmp2); % resultant of endsamp across participants
%                 
%              angle_start_end(3,ang) = (circ_dist(circ_mean(tmp),pred_angle(ang))); % phase distance from endsamp to predicted angle
% 
%              sem_end(1,ang) = circ_std(tmp)./sqrt(length(incl_sub)); % sem of end phase
%              sem_end(2,ang) = circ_std(circ_dist(angs(1,ang),tmp))./sqrt(length(incl_sub)); % sem of difference between angle at center of bin and end phase
% 
%              A = A+1;
%              
%              end
%              
%              [B I] = sort((angle_start_end(2,:)));
%              
%              colorz = linspecer(length(angle_start_end(3,:)));
% 
%              subplot(5,4,G)
%              hold on
%              errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color','k','LineStyle','none') % angle at center of bin vs mean phase at end, errorbar
%              scatter(angle_start_end(1,:),pred_angle,2.5,'k','filled') % angle at center of bin vs predicted angle
%              
%              for i = 1:length(angle_start_end(3,:))
% 
%                 if chan == 2
% 
%                     scatter(angle_start_end(1,i),angle_start_end(2,i),25,colorz(i,:),'filled','MarkerEdgeColor','k') % angle at center of bins vs mean at end 
%                 else
%                     scatter(angle_start_end(1,:),angle_start_end(2,:),25,colors(2,:))
%                     errorbar(angle_start_end(1,:),angle_start_end(2,:),sem_end(1,:),'Color',colors(2,:),'LineStyle','none') 
%                 end
%              end
%              
%              G = G+1;
% 
%              ylim([-pi pi])
%              xlim([-pi pi])
%              title(['Delay: ',num2str(cycs/fs*1000),' ms'])     
% 
%              xlabel('Start Phase')
%              ylabel('End Phase')
%              
%              subplot(5,4,G)
%              hold on
%              errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color','k','LineStyle','none') % angle at center of bin vs phase distance from endsamp to predicted angle, errorbar
%              
%              for i = 1:length(angle_start_end(3,:))
%                 if chan == 2
%                     scatter(angle_start_end(1,i),angle_start_end(3,i),25,colorz(i,:),'filled','MarkerEdgeColor','k')                         
%                 else
%                     scatter(angle_start_end(1,:),angle_start_end(3,:),25,colors(2,:))
%                     errorbar(angle_start_end(1,:),angle_start_end(3,:),sem_end(2,:),'Color',colors(2,:),'LineStyle','none')    
%                 end
%              end
%              
%             G = G+1;
%             a = 1;
%             for A = [0 90 180 270]
% %             for A = [340 60 160 240]
%                 xline(wrapToPi(deg2rad(A)),'Color',colors(a,:)) 
%                 a = a+1;
%             end
%              %xlim([-180 180])
%              xlabel('Start Phase')
%              ylabel('Phase Change')
%              yline(0)
%              text(.95*-pi,-.9*pi,'slower')
%              %ylim([-200 200])
%              %ylim([.1 .14])
%              ylim([-pi pi])
%              xlim([-pi pi])
%              text(.95*-pi,.9*pi,'faster')
%                          
%                           
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(2,:)')],[0 (circ_r(angle_start_end(2,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%             if chan == 2
%                 B = 0.5;
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(2,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     B = B+.025;
%                     hold on
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%   
%                         hold on
%                     end
%                 end
%                 
%             else
%                 polarplot(angle_start_end(1,:),(circ_dist(circ_mean(angle_start_end(3,:)'),angle_start_end(3,:)))+pi,'Color',colors(2,:))
% 
%             end
%              hold on
% 
%              rlim([0 1.1])
%              G = G+1;
%              
%              subplot(5,4,G)
%              polarplot([0 circ_mean(angle_start_end(3,:)')],[0 (circ_r(angle_start_end(3,:)',R_end'))],'k','LineWidth',2)
%              hold on
% 
%              if chan == 2
%                 B = .5;
%                 
%                 for i = 1:length(angle_start_end(3,:))
%                     polarscatter(angle_start_end(3,i),R_end(i),35,colorz(i,:),'filled','MarkerEdgeColor','k')
%                     hold on
%                     B = B+.025;
%     
%                     if G == 3
%                         before(i) = angle_start_end(2,i);
%                     else
%                         hold on
%                     end
%                 end
%              end
%              
%              hold on
%              rlim([0 1.1])
%              G = G+1;
%              
% end
%              
%              
% % saveas(fig,[Savefolder,'Figure6_phase_response_curve_freq',num2str(osc_freq),'_wake_mor_EC.svg']);
% 
%              
%              
% 
% 
%              
%              
% 
