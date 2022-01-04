%% MS_initialize: initialize the default parameters and

% clear all;
close all
restoredefaultpath
global PARAMS

if isunix
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Multisite/'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Multisite/Temp/'; % where to put intermediate files
    %     PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Multisite/Stats/'; % where to put the statistical output .txt
    PARAMS.stats_dir = '/Volumes/Fenrir/MS_temp/Stats2/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = '/Users/jericcarmichael/Documents/GitHub/EC_Multisite'; % where the multisite repo can be found
    PARAMS.Chronux_dir = '/Users/jericcarmichael/Documents/chronux_2_11';
elseif strcmp(computer, 'GLNXA64')
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Multisite/'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Multisite/Temp/'; % where to put intermediate files
    %     PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Multisite/Stats/'; % where to put the statistical output .txt
    PARAMS.stats_dir = '/Volumes/Fenrir/MS_temp/Stats2/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = '/home/ecarmichael/Documents/GitHub/EC_Multisite'; % where the multisite repo can be found
    
else
    PARAMS.data_dir = 'G:\Multisite\'; % where to find the raw data
    PARAMS.inter_dir = 'G:\Multisite\temp\'; % where to put intermediate files
    PARAMS.stats_dir = 'G:\Multisite\temp\Stats\'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite'; % where the multisite repo can be found
end
% log the progress
PARAMS.log = fopen([PARAMS.data_dir 'MS_log.txt'], 'w');
% define subjects, phases,
PARAMS.Phases = {'pre', 'ipsi', 'contra', 'post'}; % recording phases within each session
PARAMS.Subjects = {'R102', 'R104','R107', 'R108', 'R112', 'R122', 'R123'}; %list of subjects
PARAMS.Sub_xpiri = {'R108', 'R112', 'R122', 'R123'};  % subjects with electrodes spanning the piriform cortex
PARAMS.all_sites = {'PL', 'IL', 'OFC', 'PiriO', 'NAc', 'PiriN', 'CG'};
rng(10,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_MS_dir));
cd(PARAMS.data_dir) % move to the data folder

% Set of colours for consistency between plots
PARAMS.all_pairs = {'IL_PL' 'IL_OFC', 'IL_NAc', 'IL_CG', 'IL_PiriO', 'IL_PiriN','PL_OFC', 'PL_NAc', 'PL_CG', 'PL_PiriO', 'PL_PiriN', ...
    'OFC_NAc', 'OFC_CG', 'OFC_PiriO', 'OFC_PiriN', 'NAc_CG', 'NAc_PiriO', 'NAc_PiriN','CG_PiriO', 'CG_PiriN', 'PiriO_PiriN'};
PARAMS.pair_c_ord = linspecer(length(PARAMS.all_pairs)+1);
PARAMS.pair_c_ord(5,:) = [];

PARAMS.fig_blue = [0,173,216]./255;
PARAMS.fig_green = [123,225,160]./255;
PARAMS.fig_pink  = [255,168,213]./255;

% formatting
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'DefaultLegendInterpreter', 'none')
set(groot,'defaulttextinterpreter','none');

%MS_amp_xcorr_session
%Master_Multisite_preprocess;

