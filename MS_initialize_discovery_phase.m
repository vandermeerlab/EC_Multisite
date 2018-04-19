%% MS_initialize: initialize the default parameters and 
tic
clear all; close all
restoredefaultpath
global PARAMS

if isunix
    PARAMS.data_dir = '/global/scratch/ecarmichael/Multisite/'; % where to find the raw data
    PARAMS.inter_dir = '/dartfs-hpc/rc/lab/M/MeerM/EC/temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/dartfs-hpc/rc/home/r/f00287r/MS/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/dartfs-hpc/rc/home/r/f00287r/Code/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = '/dartfs-hpc/rc/home/r/f00287r/Code/EC_Multisite'; % where the multisite repo can be found
    
else
    PARAMS.data_dir = 'G:\JK_recordings\Naris\'; % where to find the raw data
    PARAMS.inter_dir = 'G:\JK_recordings\Naris\Multisite\temp\'; % where to put intermediate files
    PARAMS.stats_dir = 'G:\JK_recordings\Naris\Multisite\temp\Stats\'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite'; % where the multisite repo can be found
end
% log the progress
PARAMS.log = fopen([PARAMS.inter_dir 'MS_log.txt'], 'w');
% define subjects, phases, 
PARAMS.Phases = {'pre', 'ipsi', 'contra', 'post'}; % recording phases within each session
PARAMS.Subjects = {'R102', 'R104', 'R107', 'R108', 'R112', 'R122', 'R123'}; %list of subjects
PARAMS.Sub_xpiri = {'R108', 'R112', 'R122', 'R123'};  % subjects with electrodes spanning the piriform cortex
PARAMS.all_sites = {'PL', 'IL', 'OFC', 'Piri_O', 'NAc', 'Piri_N', 'CG'}; 

% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_MS_dir));
cd(PARAMS.data_dir) % move to the data folder

% formatting
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
set(groot, 'DefaultLegendInterpreter', 'none')
set(groot,'defaulttextinterpreter','none');  
Master_Phase