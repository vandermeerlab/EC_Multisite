% Multi_Master:
%   Master control script for multi-site loading, analysis, and statistics.
%   Requirements:
%        - van der Meer lab codebase
%        - EC_Multisite functions (vandermeerlab/EC_Multisite on github)
%


% terminology:
%    - subject: each rat used.  Two experiments "multisite" and xpiri
%
%    - session: recording day which included four phases
%        - phase : recording phases of the naris protocol as "pre", "ipsi",
%                  "contra", "post"
%

%% make a log
MS_log = fopen([PARAMS.data_dir '/MS_log.txt'], 'w');
fprintf(MS_log, date);
%% Extract the data from each recroding phase within each session and separate pot vs track sections
%
% for iSub = 1:length(PARAMS.Subjects)
%     if isunix
%         cd([PARAMS.data_dir '/' PARAMS.Subjects{iSub}])
%     else
%         cd([PARAMS.data_dir '\' PARAMS.Subjects{iSub}])
%     end
%
%     dir_files = dir(); % get all the sessions for the current subject
%     dir_files(1:2) = [];
%     sess_list = [];
%     for iDir = 1:length(dir_files)
%         if dir_files(iDir).isdir == 1 && strcmp(dir_files(iDir).name(1:4), PARAMS.Subjects{iSub}) % ensure that the session folders have the correct names
%             sess_list = [sess_list;dir_files(iDir).name];  % extract only the folders for the seesions
%         end
%     end
%     sess_list = cellstr(sess_list);
%     % load the data for each session within the current subject.
%     for iSess = 1:length(sess_list)
%         cfg_loading = [];
%         cfg_loading.fname = sess_list{iSess};
%         fprintf(['\n' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
%         fprintf(MS_log,['\nLoading ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
%
%         [data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), cfg_loading] = MS_load_data_fast(cfg_loading);
%         fprintf(MS_log, '...complete');
%     end
%     % ensure the correct number of sessions exist per rat
%     if length(fieldnames(data.(PARAMS.Subjects{iSub}))) ~=4
%         error('too many or too few sessions for multisite experiment.  Should only contain 4 per rat')
%     end
% end
%
% %% get the gamma event counts per recording phase
% fprintf(MS_log,'\n\nCollecting Events');
%
% for iSub = 1:length(PARAMS.Subjects)
%     sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
%     for iSess  = 1:length(sess_list)
%         fprintf(MS_log,['\nEvents ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
%         [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))] = MS_extract_gamma([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%         fprintf(MS_log, '...complete');
%     end
% end
% % summary of naris events
%
% %     stats = MS_gamma_stats([], Events);
% %% generate PSDs & get the relative power ratios
% fprintf(MS_log,'\n\nExtracting Power Metrics');
% for iSub = 1:length(PARAMS.Subjects)
%     sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
%     for iSess  = 1:length(sess_list)
%         fprintf(['Session ' sess_list{iSess} '\n'])
%         fprintf(MS_log,['\nGetting Power ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
%         Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_collect_psd([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%         fprintf(MS_log, '...complete');
%     end
% end
%
% %% save the intermediate files
% fprintf(MS_log,'\n\nSaving intermediates');
% mkdir(PARAMS.data_dir, 'temp');
% save([PARAMS.data_dir 'MS_data.mat'], 'data', '-v7.3')
% save([PARAMS.data_dir 'MS_naris.mat'], 'Naris', '-v7.3')
% save([PARAMS.data_dir 'MS_events.mat'], 'Events', '-v7.3')

% fclose(MS_log);

%% load the intermediate files
MS_log = fopen([PARAMS.data_dir '/MS_log_2.txt'], 'w');
fprintf(MS_log, date);
fprintf(MS_log,'\n\nLoading intermediates');
mkdir(PARAMS.data_dir, 'temp');
load([PARAMS.data_dir '/MS_data.mat'])
load([PARAMS.data_dir '/MS_naris.mat'])
load([PARAMS.data_dir '/MS_events.mat'])

%% split pot vs trk

[Naris_pot, Naris_trk]  = MS_pot_trk_split(Naris);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% get the ratio of the power in multiple bands relative to the exponential f curve
% % fprintf(MS_log,'\n\nExtracting Power Ratio');
% for iSub = 1:length(PARAMS.Subjects)
%     sess_list = fieldnames(Naris.(PARAMS.Subjects{iSub}));
%     for iSess  = 1:length(sess_list)
%         %         fprintf(MS_log,['\nGetting ratio ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
%         cfg_pow_ratio.id = sess_list{iSess};
%         Naris_pot.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_get_power_ratio(cfg_pow_ratio,Naris_pot.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%         Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_get_power_ratio(cfg_pow_ratio,Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%         %         fprintf(MS_log, '...complete');
%     end
% end
%% plot the PSDs
MS_plot_psd([], Naris);
%% get an example event from each session and plot all sites together for the same event.
for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        fprintf(MS_log,['\nPlotting Events ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        MS_event_fig([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(MS_log, '...complete');
    end
end


%% generate a spectrogram across each session for each site.
for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        fprintf(MS_log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        MS_spec_fig([], data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(MS_log, '...complete');
    end
end
% %% plot the gamma band power ratios
% cfg_pow_ratio_plot.ylims = [-75 75];
% cfg_pow_ratio_plot.plot_type = 'raw';
% cfg_pow_ratio_plot.ylims_norm = [0 2];
%
% MS_plot_power_ratio(cfg_pow_ratio_plot, Naris_pot)
% MS_plot_power_ratio(cfg_pow_ratio_plot, Naris_trk)

% MS_plot_power([], Naris);

%% Get the phase coherence metrics
% create pairs of channels for detected events.

for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), Coh_mat.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))]  = MS_event_pairs([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    end
end

%% plot the COH metrics

stats_coh =  MS_Coh_plot_stats(Coh_mat);
%
%
%
% %% get the coordinates from the Expkeys
%
% % stats_subjects = MS_get_subject_info(data);
fclose(MS_log);
