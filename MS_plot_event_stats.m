function MS_plot_event_stats(cfg_in, Events)
%% counts the number of events for each phase of the Naris experiment.  Runs stats on the events as well.




%% defaults and parameters
global PARAMS

cfg_def = [];


cfg = ProcessConfig2(cfg_def, cfg_in);

%check for sites
if ~isfield(cfg, 'sites') || isempty(cfg.sites)
    error('No sites specified in cfg')
end

Subjects = fieldnames(Events); % get the list of subjects
%% cycle through the sites specified in the input cfg
bands = {'low', 'high'};
phases = PARAMS.Phases;
phases{5} = 'control';

% predefine the event length structure
for iBand = 1:2
    for iSite = 1:length(cfg.sites)
        for iPhase = 1:length(phases)
            all_events_len.(bands{iBand}).(cfg.sites{iSite}).(phases{iPhase}) =[];
        end
    end
end
                    
for iBand = 1:2
    all_out.(bands{iBand}).rate = []; all_out.(bands{iBand}).num = [];
    for iSub = 1:length(Subjects)
        % cycle through sessions per subject
        sess_list = fieldnames(Events.(Subjects{iSub}));
        for iSess = 1:length(sess_list)
            for iSite = 1:length(cfg.sites)
                for iPhase = 1:length(phases)
                    label{iSite, iPhase} = [cfg.sites{iSite} '-' phases{iPhase}];
                    % collect the events as [Sites, Phase, session]
                    if  isfield(Events.(Subjects{iSub}).(sess_list{iSess}), cfg.sites{iSite})
                        single_subject.rate(iSite, iPhase, iSess) = Events.(Subjects{iSub}).(sess_list{iSess}).(cfg.sites{iSite}).(phases{iPhase}).(bands{iBand}).rate;
                        single_subject.num(iSite, iPhase, iSess) = length(Events.(Subjects{iSub}).(sess_list{iSess}).(cfg.sites{iSite}).(phases{iPhase}).(bands{iBand}).tstart);
                        
                        % get the length of every event 
                        evt_lens = Events.(Subjects{iSub}).(sess_list{iSess}).(cfg.sites{iSite}).(phases{iPhase}).(bands{iBand}).tend - Events.(Subjects{iSub}).(sess_list{iSess}).(cfg.sites{iSite}).(phases{iPhase}).(bands{iBand}).tstart;
                        all_events_len.(bands{iBand}).(cfg.sites{iSite}).(phases{iPhase}) = cat(1,all_events_len.(bands{iBand}).(cfg.sites{iSite}).(phases{iPhase}),evt_lens);
                        if sum(evt_lens>1) >=1
                            disp([(Subjects{iSub}) '-' (sess_list{iSess}) '-' (bands{iBand}) '-' (cfg.sites{iSite}) '-' (phases{iPhase})])
                            pause
                        end
                    else
                        single_subject.rate(iSite, iPhase, iSess) = NaN;
                        single_subject.num(iSite, iPhase, iSess) = NaN;
                    end
                end
            end
        end
        all_subs.(Subjects{iSub}).(bands{iBand}) = single_subject;
        all_out.(bands{iBand}).rate = cat(3,all_out.(bands{iBand}).rate, single_subject.rate);
        all_out.(bands{iBand}).num = cat(3,all_out.(bands{iBand}).num, single_subject.num);
    end
end


%% compile everything

% for iBands = 1:2
%     Mean_rate.(bands{iBand}) = nanmean(all_out.(bands{iBand}).rate, 3);
%
%     Mean_num.(bands{iBand}) = nanmean(all_out.(bands{iBand}).num, 3);
%
%


%% print all the stats
fprintf('\n\n\nGamma Event Stats\n')
fileID = fopen(cat(2,PARAMS.stats_dir,'\Naris_stats_events.txt'),'w');
fprintf(fileID,'Gamma Event Stats\n')
fprintf(fileID, ['_________________________________________\n'])
fprintf(fileID, [date '\n'])

for iband = 1:length(bands)
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
end

fclose(fileID);