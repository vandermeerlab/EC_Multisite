% Multi_Master():
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

%% initialize default PARAMeters
clear all; close all
restoredefaultpath
global PARAMS

if isunix
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Multisite'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Multisite/Temp/'; % where to put intermediate files
    PARAMS.stats_out = '/Users/jericcarmichael/Documents/Multisite/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = '/Users/jericcarmichael/Documents/GitHub/EC_Multisite'; % where the multisite repo can be found
    
else
    PARAMS.data_dir = 'G:\JK_recordings\Naris'; % where to find the raw data
    PARAMS.inter_dir = 'G:\JK_recordings\Naris\Multisite\temp\'; % where to put intermediate files
    PARAMS.stats_out = 'G:\JK_recordings\Naris\Multisite\temp\Stats\'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite'; % where the multisite repo can be found
end

PARAMS.Phases = {'pre', 'ipsi', 'contra', 'post'}; % recording phases within each session
PARAMS.Subjects = {'R102', 'R104','R107', 'R108', 'R112', 'R122', 'R123'}; %list of subjects
PARAMS.Sub_xpiri = {'R108', 'R112', 'R122', 'R123'};  % subjects with electrodes spanning the piriform cortex

% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_MS_dir));
cd(PARAMS.data_dir) % move to the data folder

%% Extract the data from each recroding phase within each session and separate pot vs track sections

for iSub = 3:5%:length(PARAMS.Subjects)
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
        [data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), cfg_loading] = MS_load_data_fast(cfg_loading);
    end
    % ensure the correct number of sessions exist per rat
    if length(fieldnames(data.(PARAMS.Subjects{iSub}))) ~=4
        error('too many or too few sessions for multisite experiment.  Should only contain 4 per rat')
    end
end

%% get the gamma event counts per recording phase
for iSub = 4:5%1:length(PARAMS.Subjects)
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess  = 1:length(sess_list)
        [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))] = MS_extract_gamma([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    end
end
    % summary of naris events
    
%     stats = MS_gamma_stats([], Events);
%% generate PSDs & get the relative power ratios
for iSub = 3:5%1:length(PARAMS.Subjects)
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess  = 1:length(sess_list)
        fprintf(['Session ' sess_list{iSess} '\n'])
        Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_collect_psd([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    
%         Naris_pot.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_get_power_ratio([],Naris_pot.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    end
end

%% plot the gamma band power ratios
MS_plot_power([], Naris);


%% Get the phase coherence metrics
%%%%%%%%%%%%%%%% you are here %%%%%%%%%%%%%%%%%%%%


% create pairs of channels for detected events.  

for iSub = 3%:length(PARAMS.Subjects)
    sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), Coh_mat.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))]  = MS_event_pairs([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    end    
end

%% plot the COH metrics 

stats_coh =  MS_Coh_plot_stats(Coh_mat);

%% plot the PSDs
MS_plot_psd([], Naris);



%% get an example event from each session and plot all sites together for the same event. 

for iSub = 3:5%1:length(PARAMS.Subjects)
% iSess = 1; iSub = 1;     
sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        MS_event_fig([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    end    
end

%% generate a spectrogram across each session for each site. 

for iSub = 3:5%1:length(PARAMS.Subjects)
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        MS_spec_fig([], data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
    end    
end


%% get the coordinates from the Expkeys

stats_subjects = MS_get_subject_info(data);
