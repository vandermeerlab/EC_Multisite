% PPC_init

%% used for cluster runs of the PPC_comparison script.  Paths and parameters are put into a cfg struct here.  

% addpath('/dartfs-hpc/rc/home/r/f00287r/Code/EC_NVHL')
% addpath(genpath('/dartfs-hpc/rc/home/r/f00287r/Code/vandermeerlab/code-matlab/shared'))
% addpath('/dartfs-hpc/rc/home/r/f00287r/Code/fieldtrip')
addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'))

addpath('/Users/jericcarmichael/Documents/GitHub/fieldtrip')
ft_defaults
addpath(genpath('/Users/jericcarmichael/Documents/GitHub/EC_Multisite'))
ft_defaults

cfg_in = [];
cfg_in.dataset = {'CSC12.ncs'};
cfg_in.data_dir = '/Users/jericcarmichael/Documents/Multisite/R112-2017-08-01_post';
cfg_in.inter_dir = '/Users/jericcarmichael/Documents/Multisite/R112-2017-08-01_post';
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST6_cut'; 

PPC_comparison_CG(cfg_in); 

%%
disp('108 ST3 ...')
cfg_in = [];
cfg_in.dataset = {'CSC5.ncs'};
cfg_in.data_dir = '/Users/jericcarmichael/Documents/Multisite/R108-2017-08-04_pre';
cfg_in.inter_dir = '/Users/jericcarmichael/Documents/MultisiteR108-2017-08-04_pre';
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST3_cut'; 
PPC_comparison_CG(cfg_in); 
disp('complete')
