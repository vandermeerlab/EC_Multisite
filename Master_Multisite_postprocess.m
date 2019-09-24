% Multi_Master_postprocess:
%   Master control script for multi-site loading, analysis, and statistics.
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
% Extract the data from each recroding phase within each session and separate pot vs track sections
cd(PARAMS.inter_dir)

%% loop through subjects
for iSub = 1%1:length(PARAMS.Subjects)

% iSub = 4;
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'])
    
    %%%%%%%%%%%%%%% Power Measures %%%%%%%%%%%%%%%%%
    
    
    % generate the PSD for each session for each site
    fprintf(PARAMS.log,'\n\nExtracting Power Metrics');
    sess_list = fieldnames(data);
    for iSess  = 1:length(sess_list)
        fprintf(['Session ' sess_list{iSess} '\n'])
        fprintf(PARAMS.log,['\nGetting Power ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_collect_psd([],data.(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end

    % get the ratio of the power in multiple bands relative to the exponential f curve
    fprintf(PARAMS.log,'\n\nExtracting Power Ratio');
    sess_list = fieldnames(Naris.(PARAMS.Subjects{iSub}));
    for iSess  = 1:length(sess_list)
        fprintf(PARAMS.log,['\nGetting ratio ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        cfg_pow_ratio = [];
        if iSub == 6 || iSub == 7
            cfg_pow_ratio.Fs = 1875;
        else
            cfg_pow_ratio.Fs = 2000;
        end
        cfg_pow_ratio.id = sess_list{iSess};
        [Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), cfg_p_ratio] = MS_get_power_ratio(cfg_pow_ratio,Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        %         Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_get_power_ratio(cfg_pow_ratio,Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
        close all
    end
% end


%     d_t = data;
%     clear data
%     data.(PARAMS.Subjects{iSub}) = d_t;
%     clear d_t
save([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris.mat'], 'Naris', '-v7.3')
%     if exist([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat']) ==2
%         load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'])
%     end

%% get amplitude xcorr for each session .
%for iSub = length(PARAMS.Subjects):-1:5
    sess_list = fieldnames(data);
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).amp = MS_amp_xcorr_session_2([], data.(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
%    end

%% Get the coherence across each session.
% for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(data);
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).coh = MS_coh_session([], data.(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
% end

%% save the intermediate files
fprintf(PARAMS.log,'\n\nSaving intermediates');
% mkdir(PARAMS.data_dir, 'temp');
% save([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'], 'data', '-v7.3')
save([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'], 'Naris', '-v7.3')

clearvars -except iSub PARAMS

end

%% Collect the phase information in an 'All' structure
MS_collect_phase()



    
    % PPC for cut sessions (did not use)
    
    

% end

