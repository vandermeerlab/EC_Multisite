% PPC_init

%% used for cluster runs of the PPC_comparison script.  Paths and parameters are put into a cfg struct here.  

addpath(genpath('/dartfs-hpc/rc/home/r/f00287r/Code/EC_Multisite'))
addpath(genpath('/dartfs-hpc/rc/home/r/f00287r/Code/vandermeerlab/code-matlab/shared'))
addpath('/dartfs-hpc/rc/home/r/f00287r/Code/fieldtrip')
% addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'))

% addpath('/Users/jericcarmichael/Documents/GitHub/fieldtrip')
% ft_defaults
% addpath(genpath('/Users/jericcarmichael/Documents/GitHub/EC_Multisite'))
ft_defaults
data_root = '/dartfs-hpc/rc/lab/M/MeerM/EC';
mkdir(data_root, 'PPC_figs')
nShuffle = 100; 
%%

disp('112 ST1 ...')
cfg_in = [];
cfg_in.dataset = {'CSC1.ncs'};
cfg_in.data_dir = [data_root '/R112-2017-07-31_post'];
cfg_in.inter_dir = [data_root '/PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = nShuffle; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST1_SS'; 
PPC_comparison_CG(cfg_in); 
disp('complete')
%%
disp('112 St1 ...')
cfg_in = [];
cfg_in.dataset = {'CSC1.ncs'};
cfg_in.data_dir = [data_root '/R112-2017-08-01_post'];
cfg_in.inter_dir = [data_root '/PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = nShuffle; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST1_cut'; 
PPC_comparison_CG(cfg_in); 
disp('complete')

%%
disp('112 ST6 ...')
cfg_in = [];
cfg_in.dataset = {'CSC12.ncs'};
cfg_in.data_dir = [data_root '/R112-2017-08-01_post'];
cfg_in.inter_dir = [data_root '/PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = nShuffle; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST6_cut'; 
%PPC_comparison_CG(cfg_in); 
%disp('complete')


%%
disp('108 ST3 ...')
%cfg_in = [];
%cfg_in.dataset = {'CSC5.ncs'};
%cfg_in.data_dir = [data_root '/R108-2017-08-04_pre'];
%cfg_in.inter_dir = [data_root '/PPC_figs'];
%cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
%cfg_in.shuffle = nShuffle; 
%cfg_in.plot = 1;
%cfg_in.min_nSpikes = 500; 
%cfg_in.spike_id = 'ST3_cut'; 
%PPC_comparison_CG(cfg_in); 
%disp('complete')

%%
disp('108 TT7 ...')
cfg_in = [];
cfg_in.dataset = {'CSC16.ncs'};
cfg_in.data_dir = [data_root '/R108-2017-08-04_pre'];
cfg_in.inter_dir = [data_root '/PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = nShuffle; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT7_cut'; 
PPC_comparison_CG(cfg_in); 
disp('complete')

%%
disp('123 TT1 ...')
cfg_in = [];
cfg_in.dataset = {'CSC1.ncs'};
cfg_in.data_dir = [data_root '/R123-2017-02-27_pre2'];
cfg_in.inter_dir = [data_root '/PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = nShuffle; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT1_cut'; 
PPC_comparison_CG(cfg_in); 
disp('complete')

%%
disp('123 TT4 ...')
cfg_in = [];
cfg_in.dataset = {'CSC16.ncs'};
cfg_in.data_dir = [data_root '/R123-2017-02-27_pre2'];
cfg_in.inter_dir = [data_root '/PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = nShuffle; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT4_cut'; 
PPC_comparison_CG(cfg_in); 
disp('complete')
