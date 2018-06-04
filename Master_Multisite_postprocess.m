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

for iSub = 1:length(PARAMS.Subjects)
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'])
	d_t = data;
clear data
data.(PARAMS.Subjects{iSub}) = d_t;
clear d_t


%% generate a Coherogram across each session for each site.
for iSub = 1:length(PARAMS.Subjects)
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).amp = MS_amp_xcorr_session_2([], data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
end



%% save the intermediate files
fprintf(PARAMS.log,'\n\nSaving intermediates');
% mkdir(PARAMS.data_dir, 'temp');
% save([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'], 'data', '-v7.3')
save([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Naris_amp.mat'], 'Naris', '-v7.3')

clearvars -except iSub PARAMS
end

