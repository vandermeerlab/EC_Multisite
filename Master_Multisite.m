% Multi_Master():
%   Master control script for multi-site loading, analysis, and statistics.
%   Requirements: 
%        - van der Meer lab codebase
%        - EC_Multisite functions (vandermeerlab/EC_Multisite on github)
%



%% initialize default PARAMeters
global PARAMS

if isunix
PARAMS.data_dir = '/Users/jericcarmichael/Documents/Data'; % where to find the raw data
PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Data/Temp'; % where to put intermediate files
PARAMS.stats_out = '/Users/jericcarmichael/Documents/Data/Stats'; % where to put the statistical output .txt

else
PARAMS.data_dir = '/Users/jericcarmichael/Documents/Data'; % where to find the raw data
PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Data/Temp'; % where to put intermediate files
PARAMS.stats_out = '/Users/jericcarmichael/Documents/Data/Stats'; % where to put the statistical output .txt   
end

PARAMS.Phases = {'pre', 'ipsi', 'contra', 'post'}; % recording phases within each session

PARAMS.Subjects = {'102', '104', '122', '123'}; %list of subjects 
PARAMS.Sub_xpiri = {'122', '123'};  % subjects with electrodes spanning the piriform cortex 
%% Extract the data from each recroding phase within each session and separate pot vs track sections

cfg_loading = [];
[data, cfg_loading] = MS_load_data(cfg_loading); 
