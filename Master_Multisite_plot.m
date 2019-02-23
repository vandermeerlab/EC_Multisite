% Multi_Master_plot:
%   Master control script for multi-site plotting.
%   Requirements:
%        - van der Meer lab codebase
%        - EC_Multisite functions (vandermeerlab/EC_Multisite on github)
%


% terminology:
%
%    - session: recording day which included four phases
%        - phase : recording phases of the naris protocol as "pre", "ipsi",
%                  "contra", "post"
%

%% make a log
global PARAMS
fprintf(PARAMS.log, date);
PARAMS.inter_dir = '/Volumes/Fenrir/MS_temp/';

% Extract the data from each recroding phase within each session and separate pot vs track sections
close all
for iSub = 1:length(PARAMS.Subjects)
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'])
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris.mat'])
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Events.mat'])

    %% generate sample events for each session for each site.

        sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Events ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        % once for low events
        cfg_in = [];
        cfg_in.type = 'low';
        MS_event_fig(cfg_in, Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(strrep(sess_list{iSess}, '-', '_')));
        %once for high events
        cfg_in = [];
        cfg_in.type = 'high';
        MS_event_fig(cfg_in, Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
    
    %% generate a spectrogram across each session for each site.

        sess_list = fieldnames(data);
        for iSess = 1:length(sess_list)
            fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            MS_spec_fig([], data.(strrep(sess_list{iSess}, '-', '_')));
            fprintf(PARAMS.log, '...complete');
        end
    
    %% generate a coher-o-gram for each session used in Fig 6 (very time consuming.  Run if needed) Session R107-2018-08-04 was used in paper
    sess_list = fieldnames(data);
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        MS_coherogram_fig([], data.(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
    
        %% generate a amplitude cross corr-o-gram for each session used in Fig 7 (very time consuming.  Run if needed) Session R107-2018-08-04 was used in paper
    sess_list = fieldnames(data);
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Ampxcorr_o_gram ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        MS_plot_amp_xcorrogram([], data.(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
    
    
    %% plot the PSD for each session and each site
    cfg_psd = [];
    %         cfg_psd.type = 'white';
    MS_plot_psd(cfg_psd, Naris.(PARAMS.Subjects{iSub}))
    
    
end

%% %%%%%%% CROSS SESSION ANALYSES  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for iSub = 1:length(PARAMS.Subjects)    
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'])
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Events.mat'])

     all_Naris.(PARAMS.Subjects{iSub}) = Naris;
     all_Events.(PARAMS.Subjects{iSub}) = Events;

    clearvars -except iSub PARAMS all_Naris all_Events
    close all
end
    
%% generate the averages for each site across sessions

% plot the average psd across sessions
MS_plot_psd_avg([], all_Naris)
close all

% plot all the power ratio statistics
cfg_pow_ratio = [];
cfg_pow_ratio.method = 'median';
cfg_pow_ratio.stats_method = 'lme';

MS_plot_power_ratio(cfg_pow_ratio, all_Naris)
close all

% plot all the gamma event statistics
cfg_count = [];
cfg_count.method = 'median';
cfg_count.stats_method = 'lme';
MS_plot_gamma_stats(cfg_count, all_Events)
close all
%% get the session wide coherence plots

    % %% generate a Coherogram across each session for each site.
    % fprintf(PARAMS.log,['\nPlotting Coh Sess ' PARAMS.Subjects{iSub}]);
    mkdir(PARAMS.inter_dir, 'sess')
    cfg_coh = [];
    cfg_coh.measure = 'coh';
    cfg_coh.plot_type = 'no_piri';
    if iSub ==1
        cfg_coh.legend = 'on';
    end
    MS_plot_session_phase(cfg_coh, all_Naris);
    fprintf(PARAMS.log, '...complete');
    
    close all
%% get the session wide amplitude plots
    
    fprintf(PARAMS.log,['\nPlotting Amp Sess ' PARAMS.Subjects{iSub}]);
    cfg_coh = [];
    cfg_coh.measure = 'amp';
    MS_plot_session_phase(cfg_coh, all_Naris);
    fprintf(PARAMS.log, '...complete');
        
    

%% plot the phase measures for all events
close all
cfg_event = [];
% cfg_event.Subjects = {'all'}; % for speed
cfg_event.Subjects = {'all', 'R102', 'R104','R107', 'R108', 'R112','R122','R123'}; % each subject by themselves.  Used for subpanels on F3
MS_plot_event_phase(cfg_event)








    %% generate a Coherogram across each session for each site.
    %     load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'])
    %
    %     sess_list = fieldnames(data);
    %     for iSess = 1:length(sess_list)
    %         fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
    %         MS_coherogram_fig([], data.(strrep(sess_list{iSess}, '-', '_')));
    %         fprintf(PARAMS.log, '...complete');
    %     end
    
    
    
   
%% Naris power by distance.  Used for stats only.  
cfg_dist_stats = [];
MS_get_naris_dist(cfg_dist_stats, all_Naris);


%% Figure S2 all coherence for each site pair (no piri)
cfg_S2 = [];
cfg_S2.measure = 'coh';

MS_Figure_S2(cfg_S2, all_Naris);
close all

cfg_S2 = [];
cfg_S2.measure = 'amp';

MS_Figure_S2(cfg_S2, all_Naris);

close all
% just for debugging
cfg_S2 = [];
cfg_S2.measure = 'lag';

MS_Figure_S2(cfg_S2, all_Naris);

