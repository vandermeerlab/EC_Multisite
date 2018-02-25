% for iSub = 1:length(PARAMS.Subjects)
restoredefaultpath
global PARAMS

if isunix
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Multisite/'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Multisite/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Multisite/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = '/Users/jericcarmichael/Documents/GitHub/EC_Multisite'; % where the multisite repo can be found

else
    PARAMS.data_dir = 'G:\JK_recordings\Naris\'; % where to find the raw data
    PARAMS.inter_dir = 'G:\JK_recordings\Naris\Multisite\temp\'; % where to put intermediate files
    PARAMS.stats_dir = 'G:\JK_recordings\Naris\Multisite\temp\Stats\'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite'; % where the multisite repo can be found
end
% log the progress
PARAMS.log = fopen([PARAMS.data_dir 'MS_log.txt'], 'w');
% define subjects, phases, 
PARAMS.Phases = {'pre', 'ipsi', 'contra', 'post'}; % recording phases within each session
PARAMS.Subjects = {'R102', 'R104','R107', 'R108', 'R112', 'R122', 'R123'}; %list of subjects
PARAMS.Sub_xpiri = {'R108', 'R112', 'R122', 'R123'};  % subjects with electrodes spanning the piriform cortex
PARAMS.all_sites = {'PL', 'IL', 'OFC', 'Piri_OFC', 'NAc', 'Piri_NAc', 'CG'}; 
% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_MS_dir));
cd(PARAMS.data_dir) % move to the data folder

% formatting
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'DefaultLegendInterpreter', 'none')
set(groot,'defaulttextinterpreter','none');  

load([PARAMS.inter_dir 'MS_data_R104.mat'])
% load([PARAMS.inter_dir 'MS_naris.mat'])
load([PARAMS.inter_dir 'MS_events.mat'])

iSub = 2;
sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
parfor iSess = 1:length(sess_list)
    %             fprintf(PARAMS.log,['\nExtracting phase slope ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
    mat_all{iSub,iSess} = MS_get_phase_metrics([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    
    %             fprintf(PARAMS.log, '...complete');
end
%     end


% save([PARAMS.inter_dir 'MS_mat2.mat'], 'mat_all', '-v7.3')