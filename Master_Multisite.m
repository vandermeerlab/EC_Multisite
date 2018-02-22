% % Multi_Master:
% %   Master control script for multi-site loading, analysis, and statistics.
% %   Requirements:
% %        - van der Meer lab codebase
% %        - EC_Multisite functions (vandermeerlab/EC_Multisite on github)
% %
% 
% 
% % terminology:
% %    - subject: each rat used.  Two experiments "multisite" and xpiri
% %
% %    - session: recording day which included four phases
% %        - phase : recording phases of the naris protocol as "pre", "ipsi",
% %                  "contra", "post"
% %
% 
%% make a log
fprintf(PARAMS.log, date);
% Extract the data from each recroding phase within each session and separate pot vs track sections

for iSub = 1:length(PARAMS.Subjects)
    if isunix
        cd([PARAMS.data_dir '/' PARAMS.Subjects{iSub}])
    else
        cd([PARAMS.data_dir '\' PARAMS.Subjects{iSub}])
    end

    dir_files = dir(); % get all the sessions for the current subject
    dir_files(1:2) = [];
    sess_list = [];
    for iDir = 1:length(dir_files)
        if dir_files(iDir).isdir == 1 && strcmp(dir_files(iDir).name(1:4), PARAMS.Subjects{iSub}) % ensure that the session folders have the correct names
            sess_list = [sess_list;dir_files(iDir).name];  % extract only the folders for the seesions
        end
    end
    sess_list = cellstr(sess_list);
    % load the data for each session within the current subject.
    for iSess = 1:length(sess_list)
        cfg_loading = [];
        cfg_loading.fname = sess_list{iSess};
        fprintf(['\n' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        fprintf(PARAMS.log,['\nLoading ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);

        [data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), cfg_loading] = MS_load_data_fast(cfg_loading);
        fprintf(PARAMS.log, '...complete');
    end
    % ensure the correct number of sessions exist per rat
    if length(fieldnames(data.(PARAMS.Subjects{iSub}))) ~=4
        error('too many or too few sessions for multisite experiment.  Should only contain 4 per rat')
    end
end

%% get the gamma event counts per recording phase
fprintf(PARAMS.log,'\n\nCollecting Events');

for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess  = 1:length(sess_list)
        fprintf(PARAMS.log,['\nEvents ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))] = MS_extract_gamma([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
end
% summary of naris events
% 
% %     stats = MS_gamma_stats([], Events);
%% generate PSDs 
fprintf(PARAMS.log,'\n\nExtracting Power Metrics');
for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess  = 1:length(sess_list)
        fprintf(['Session ' sess_list{iSess} '\n'])
        fprintf(PARAMS.log,['\nGetting Power ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_collect_psd([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
end

% %% save the intermediate files
fprintf(PARAMS.log,'\n\nSaving intermediates');
% mkdir(PARAMS.data_dir, 'temp');
save([PARAMS.data_dir 'MS_data.mat'], 'data', '-v7.3')
save([PARAMS.data_dir 'MS_naris.mat'], 'Naris', '-v7.3')
save([PARAMS.data_dir 'MS_events.mat'], 'Events', '-v7.3')
% 
% % fclose(PARAMS.log);

%% get the ratio of the power in multiple bands relative to the exponential f curve
fprintf(PARAMS.log,'\n\nExtracting Power Ratio');
for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(Naris.(PARAMS.Subjects{iSub}));
    for iSess  = 1:length(sess_list)
        fprintf(PARAMS.log,['\nGetting ratio ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        cfg_pow_ratio.id = sess_list{iSess};
        [Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), cfg_p_ratio] = MS_get_power_ratio(cfg_pow_ratio,Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%         Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_get_power_ratio(cfg_pow_ratio,Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
end
save([PARAMS.data_dir 'MS_naris.mat'], 'Naris', '-v7.3')


%% load the intermediate files
% PARAMS.log = fopen([PARAMS.data_dir '/PARAMS.log_2.txt'], 'w');
% fprintf(PARAMS.log, date);
% % fprintf(PARAMS.log,'\n\nLoading intermediates');
% load([PARAMS.data_dir 'MS_data.mat'])
% load([PARAMS.data_dir 'MS_naris.mat'])
% load([PARAMS.data_dir 'MS_events.mat'])

%% split pot vs trk
% fprintf(PARAMS.log,'\n\nSplitting the data into pot and trk');
% [Naris_pot, Naris_trk]  = MS_pot_trk_split(Naris);
% [data_pot, data_trk]  = MS_pot_trk_split(data);
% % [Events_pot, Events_trk]  = MS_pot_trk_split(Events);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the PSDs
cfg_psd.type = 'white'; 
MS_plot_psd(cfg_psd, Naris);

% %% count the events
% cfg_evt_plot =[];
% cfg_evt_plot.sites = {'PL_pot', 'IL_pot', 'OFC_pot', 'NAc_pot', 'CG_pot'};
% 
% MS_plot_event_stats(cfg_evt_plot, Events)

%% get an example event from each session and plot all sites together for the same event.
for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Events ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        MS_event_fig([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
end


%% generate a spectrogram across each session for each site.
for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(data_pot.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        MS_spec_fig([], data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
end
%% plot the gamma band power ratios
% cfg_pow_ratio_plot.ylims = [-100 100];
% cfg_pow_ratio_plot.plot_type = 'raw';
% cfg_pow_ratio_plot.ylims_norm = [-2 2];
% % temporary
% cfg_pow_ratio_plot.power_ratio.contrast = [25 45; 90 110];
% cfg_pow_ratio_plot.power_ratio.gamma_freq = [45 65; 70 90];
% 
% cfg_pow_ratio_plot.pot_trk = 'pot'; 
% MS_plot_power_ratio(cfg_pow_ratio_plot, Naris)

%% get the phase slope values across all subjects, sessions, pairs, events
% for iSub = 1:length(PARAMS.Subjects)
%     sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
%     for iSess = 1:length(sess_list)
%         fprintf(PARAMS.log,['\nExtracting phase slope ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
%         mat_out.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_get_phase_metrics([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%         fprintf(PARAMS.log, '...complete');
%     end
% end
% 
% 
% save([PARAMS.data_dir 'MS_mat.mat'], 'mat_out', '-v7.3')


%%
% MS_plot_power([], Naris);

% %% Get the phase coherence metrics
% % create pairs of channels for detected events.
% 
% for iSub = 1:length(PARAMS.Subjects)
%     sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
%     for iSess = 1:length(sess_list)
%         [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), Coh_mat.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))]  = MS_event_pairs([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%     end
% end
% 
% %% plot the COH metrics
% 
% stats_coh =  MS_Coh_plot_stats(Coh_mat);
%
%
%
% %% get the coordinates from the Expkeys
%
% % stats_subjects = MS_get_subject_info(data);
% fclose(PARAMS.log);
