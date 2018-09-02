function stats = MS_gamma_stats(cfg_in, Events)
%%
%
%
%
%
%



%% set up defaults

cfg_def = [];
cfg_def.pot_trk = 'pot';

cfg = ProcessConfig2(cfg_def, cfg_in)

global PARAMS

%% collect all the gamma event counts
out.low = [];  out.high = [];
out_norm.low = []; out_norm.high= [];
rate.low =[]; rate.high = [];

sites = {'PL'    'IL'    'OFC'    'Piri_O'    'NAc'    'Piri_N'    'CG'};
bands = {'low', 'high'};
phases = [PARAMS.Phases, 'control',];
sub_list = fieldnames(Events);

for iBand = 1:2
    for iSite = 1:length(sites)
        for iPhase = 1:length(phases)
            all_events_len.(bands{iBand}).(sites{iSite}).(phases{iPhase}) =[];
        end
    end
end

for iBand = 1:length(bands)
    for iSub = 1:length(sub_list)
        sess_list = fieldnames(Events.(sub_list{iSub}));
        % pre allocate temp fields
        temp = NaN(length(PARAMS.Phases)+1, length(sites),4);
        temp_norm = NaN(length(PARAMS.Phases)+1, length(sites),4);
        single_subject.rate = NaN(length(PARAMS.Phases)+1, length(sites),4);

        for iSess = 1:length(sess_list);
            these_sites = fieldnames(Events.(sub_list{iSub}).(sess_list{iSess}));
            for iSite = 1:length(these_sites)
                       site_idx = strcmp(sites, these_sites{iSite}(1:end-4)); 
%                        if sum(site_idx) ==1
                           site_idx = find(site_idx ==1);
                           site_field = strcmp(these_sites, [sites{site_idx} '_' cfg.pot_trk]); % find the correct subfield;
                         if  sum(site_field) ==1
                           site_field = find(site_field==1);
                           
                           for iPhase = 1:4
                               temp(iPhase,site_idx,iSess) = length(Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).tstart);
                               single_subject.rate(iPhase,site_idx , iSess) = Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).rate;
                                
                               evt_lens = Events.(Subjects{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).tend - Events.(Subjects{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).tstart;
                               all_events_len.(bands{iBand}).(these_sites{site_field}).(phases{iPhase}) = cat(1,all_events_len.(bands{iBand}).(these_sites{site_field}).(phases{iPhase}),evt_lens);
                           end
                           temp(5, site_idx, iSess) = nanmean([temp(1,site_idx,iSess),temp(4,site_idx,iSess)]); 
                           single_subject.rate(5,site_idx , iSess) = Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{5}).(bands{iBand}).rate;
                           evt_lens = Events.(Subjects{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{5}).(bands{iBand}).tend - Events.(Subjects{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{5}).(bands{iBand}).tstart;
                           all_events_len.(bands{iBand}).(these_sites{site_field}).(phases{5}) = cat(1,all_events_len.(bands{iBand}).(these_sites{site_field}).(phases{5}),evt_lens);
                        
                           temp_norm(:,site_idx,iSess) = temp(:,site_idx,iSess)./temp(5,site_idx,iSess);
                          
                           single_subject.num = temp;
                           single_subject.num_norm = temp_norm;

                        
                         end
%                     if iSess == 1
%                         out.(bands{iBand})(iSite, :,iSess) = temp.(bands{iBand})(iSite, :,iSess);
%                         norm.(bands{iBand})(iSite, :,iSess) = temp_norm.(bands{iBand})(iSite, :,iSess);
%                     else
%                         out.(bands{iBand}) = [out.(bands{iBand}); temp.(bands{iBand}).(site_list{iSite})(iSess, :)];
%                         norm.(bands{iBand}) = [norm.(bands{iBand}); temp_norm.(bands{iBand}).(site_list{iSite})(iSess, :)];
%                     end
            end
%             temp = []; temp_norm;
        end
        all_subs.(Subjects{iSub}).(bands{iBand}) = single_subject;

        out.(bands{iBand}) = cat(3,out.(bands{iBand}), temp);
        out_norm.(bands{iBand}) = cat(3,out_norm.(bands{iBand}), temp_norm);
        rate.(bands{iBand}) = cat(3,rate.(bands{iBand}), single_subject.rate);
        temp = []; temp_norm =[];

    end
end

%% make sure there are no infs in the normalized version

out_norm.low(out_norm.low == inf) = NaN;
out_norm.high(out_norm.high == inf) = NaN;


%% collect the output 
% out_norm.low = circshift(

stats_file = fopen([PARAMS.stats_dir 'Count_stats.txt'], 'w');
types= {'Four', 'Piri'}; 
for iTypes = 1:length(types)
    if strcmp(types{iTypes}, 'Four')
        s_idx = [1 2 3 5 7];
    else strcmp(types{iTypes}, 'Piri')
        s_idx = [3 4 5 6];
    end
cfg_stats = [];
        cfg_stats.title = ['Count low gamma ' types{iTypes}] ;
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out.low)
close all

cfg_stats = [];
        cfg_stats.title = ['Count high gamma ' types{iTypes}];
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out.high)

%normalized
cfg_stats = [];
        cfg_stats.title = ['Count low gamma norm ' types{iTypes}];
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out_norm.low)
close all

cfg_stats = [];
        cfg_stats.title = ['Count high gamma norm ' types{iTypes}];
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out_norm.high)
close all
end


%% print all the stats
fprintf('\n\n\nGamma Event Stats\n')
fileID = fopen(cat(2,PARAMS.stats_dir,'\Naris_stats_events.txt'),'w');
fprintf(fileID,'Gamma Event Stats\n')
fprintf(fileID, ['_________________________________________\n'])
fprintf(fileID, [date '\n'])

for iband = 1:length(bands)
%     for iSite =
        all_stats.(bands{iband}).avg_len = nanmean(all_stats.(bands{iband}).length);
        all_stats.(bands{iband}).std_len = nanstd(all_stats.(bands{iband}).length);
        all_stats.(bands{iband}).total = length(all_stats.(bands{iband}).length);
        all_stats.(bands{iband}).avg_rate = nanmean(all_stats.(bands{iband}).rate); % divide by the number of minutes in the session.
        all_stats.(bands{iband}).avg_nCycle = nanmean(all_stats.(bands{iband}).nCycles);
        fprintf(fileID, [bands{iband}  ':       mean length: ' num2str(all_stats.(bands{iband}).avg_len*1000)...
            ' std: ' num2str(all_stats.(bands{iband}).std_len*1000) ...
            ' nCycle: ' num2str(all_stats.(bands{iband}).avg_nCycle)...
            ' rate: ' num2str(all_stats.(bands{iband}).avg_rate)...
            ' total: ' num2str(all_stats.(bands{iband}).total) '\n']);
%     end
    end

    fclose(fileID);
end
