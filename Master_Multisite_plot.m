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
close all

for iSub = 1:length(PARAMS.Subjects)
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'])

% %% generate a Coherogram across each session for each site.
% fprintf(PARAMS.log,['\nPlotting Coh Sess ' PARAMS.Subjects{iSub}]);
% mkdir(PARAMS.inter_dir, 'sess')
% cfg_coh = [];
% cfg_coh.measure = 'coh';
% cfg_coh.plot_type = 'no_piri';
% if iSub ==1
%     cfg_coh.legend = 'on';
% end
% MS_plot_session_phase(cfg_coh, Naris);
% fprintf(PARAMS.log, '...complete');
% 
% close all
% %% Get the coherence across each session.
% 
% fprintf(PARAMS.log,['\nPlotting Amp Sess ' PARAMS.Subjects{iSub}]);
% cfg_coh = [];
% cfg_coh.measure = 'amp';
% MS_plot_session_phase(cfg_coh, Naris);
% fprintf(PARAMS.log, '...complete');

%% plot example events



%% generate a Coherogram across each session for each site.
%     load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'])
% 
%     sess_list = fieldnames(data);
%     for iSess = 1:length(sess_list)
%         fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
%         MS_coherogram_fig([], data.(strrep(sess_list{iSess}, '-', '_')));
%         fprintf(PARAMS.log, '...complete');
%     end



all_Naris.(PARAMS.Subjects{iSub}) = Naris.(PARAMS.Subjects{iSub});

clearvars -except iSub PARAMS all_Naris
close all
end


%% Figure S2 all coherence for each site pair (no piri)
cfg_S2 = [];
cfg_S2.measure = 'coh';

MS_Figure_S2(cfg_S2, all_Naris);

cfg_S2 = [];
cfg_S2.measure = 'amp';

MS_Figure_S2(cfg_S2, all_Naris);

