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
% Extract the data from each recroding phase within each session and separate pot vs track sections

for iSub = 1:length(PARAMS.Subjects)
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'])


%% generate a Coherogram across each session for each site.
fprintf(PARAMS.log,['\nPlotting Coh Sess ' PARAMS.Subjects{iSub}]);
cfg_coh = [];
cfg_coh.measure = 'coh';
MS_plot_session_phase(cfg_coh, Naris);
fprintf(PARAMS.log, '...complete');

close all
%% Get the coherence across each session.

% fprintf(PARAMS.log,['\nPlotting Amp Sess ' PARAMS.Subjects{iSub}]);
% cfg_coh = [];
% cfg_coh.measure = 'amp';
% MS_plot_session_phase(cfg_coh, Naris);
% fprintf(PARAMS.log, '...complete');

clearvars -except iSub PARAMS
end

