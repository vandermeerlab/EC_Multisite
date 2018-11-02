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

     all_Naris.(PARAMS.Subjects{iSub}) = Naris.(PARAMS.Subjects{iSub});
     all_Events.(PARAMS.Subjects{iSub}) = Events.(PARAMS.Subjects{iSub});

    clearvars -except iSub PARAMS all_Naris all_Events
    close all
end
    
%% generate the averages for each site across sessions

% plot the average psd across sessions
MS_plot_psd_avg([], all_Naris)
close all

% plot all the power ratio statistics
MS_plot_power_ratio([], all_Naris)
close all

% plot all the gamma event statistics
MS_plot_gamma_stats([], all_Events)
close all
%% get the session wide coherence and amplitude plots

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
    %% Get the coherence across each session.
    
    fprintf(PARAMS.log,['\nPlotting Amp Sess ' PARAMS.Subjects{iSub}]);
    cfg_coh = [];
    cfg_coh.measure = 'amp';
    MS_plot_session_phase(cfg_coh, all_Naris);
    fprintf(PARAMS.log, '...complete');
        
    
    
    %% generate a Coherogram across each session for each site.
    %     load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'])
    %
    %     sess_list = fieldnames(data);
    %     for iSess = 1:length(sess_list)
    %         fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
    %         MS_coherogram_fig([], data.(strrep(sess_list{iSess}, '-', '_')));
    %         fprintf(PARAMS.log, '...complete');
    %     end
    
    
    
   



%% Figure S2 all coherence for each site pair (no piri)
cfg_S2 = [];
cfg_S2.measure = 'coh';

MS_Figure_S2(cfg_S2, all_Naris);

cfg_S2 = [];
cfg_S2.measure = 'amp';

MS_Figure_S2(cfg_S2, all_Naris);


cfg_S2 = [];
cfg_S2.measure = 'lag';

MS_Figure_S2(cfg_S2, all_Naris);

