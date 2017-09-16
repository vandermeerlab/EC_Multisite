%% MS_initialize: initialize the default parameters and 
tic
clear all; close all
restoredefaultpath
global PARAMS

if isunix
    PARAMS.data_dir = '/global/scratch/ecarmichael//Multisite'; % where to find the raw data
    PARAMS.inter_dir = '/ihome/jcarmich/MS/Temp/'; % where to put intermediate files
    PARAMS.stats_out = '/ihome/jcarmich/MS/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/ihome/jcarmich/Code/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = '/ihome/jcarmich/Code/EC_Multisite'; % where the multisite repo can be found
    
else
    PARAMS.data_dir = 'G:\JK_recordings\Naris'; % where to find the raw data
    PARAMS.inter_dir = 'G:\JK_recordings\Naris\Multisite\temp\'; % where to put intermediate files
    PARAMS.stats_out = 'G:\JK_recordings\Naris\Multisite\temp\Stats\'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_MS_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite'; % where the multisite repo can be found
end

PARAMS.Phases = {'pre', 'ipsi', 'contra', 'post'}; % recording phases within each session
PARAMS.Subjects = {'R102', 'R104', 'R108', 'R112', 'R122', 'R123'}; %list of subjects
PARAMS.Sub_xpiri = {'R108', 'R112', 'R122', 'R123'};  % subjects with electrodes spanning the piriform cortex

% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_MS_dir));
cd(PARAMS.data_dir) % move to the data folder

%% run the analysis
Master_Multisite
toc
disp('this worked')