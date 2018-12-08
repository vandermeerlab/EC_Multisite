% PPC_init

%% used for cluster runs of the PPC_comparison script.  Paths and parameters are put into a cfg struct here.  

% addpath('/dartfs-hpc/rc/home/r/f00287r/Code/EC_NVHL')
% addpath(genpath('/dartfs-hpc/rc/home/r/f00287r/Code/vandermeerlab/code-matlab/shared'))
% addpath('/dartfs-hpc/rc/home/r/f00287r/Code/fieldtrip')
if isunix
    addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'))
    addpath('/Users/jericcarmichael/Documents/GitHub/fieldtrip')
    ft_defaults
    addpath(genpath('/Users/jericcarmichael/Documents/GitHub/EC_Multisite'))
    ft_defaults
%     data_root = 

else
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'))
    addpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip')
    ft_defaults
    addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite'))
    ft_defaults
    data_root = 'G:\Multisite\temp\'; 
    mkdir(data_root, 'PPC_figs')

end



%% R108
cfg_in = [];
cfg_in.dataset = {'CSC4.ncs'};
cfg_in.data_dir = 'R108-2017-08-04_post';
cfg_in.inter_dir = 'R108-2017-08-04_post';
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST3_cut'; 

PPC_comparison_CG(cfg_in); 

%%
% disp('108 ST3 ...')
% cfg_in = [];
% cfg_in.dataset = {'CSC5.ncs'};
% cfg_in.data_dir = '/Users/jericcarmichael/Documents/Multisite/R108-2017-08-04_pre';
% cfg_in.inter_dir = '/Users/jericcarmichael/Documents/MultisiteR108-2017-08-04_pre';
% cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
% cfg_in.shuffle = 100; 
% cfg_in.plot = 1;
% cfg_in.min_nSpikes = 500; 
% cfg_in.spike_id = 'ST3_cut'; 
% PPC_comparison_CG(cfg_in); 
% disp('complete')
% %% R112
% 
% disp('112 ST1 ...')
% cfg_in = [];
% cfg_in.dataset = {'CSC1.ncs'};
% cfg_in.data_dir = [data_root '/R112-2017-07-31_post'];
% cfg_in.inter_dir = [data_root '/PPC_figs'];
% cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
% cfg_in.shuffle = nShuffle; 
% cfg_in.plot = 1;
% cfg_in.min_nSpikes = 500; 
% cfg_in.spike_id = 'ST1_SS'; 
% PPC_comparison_CG(cfg_in); 
% disp('complete')
% %%
% disp('112 St1 ...')
% cfg_in = [];
% cfg_in.dataset = {'CSC1.ncs'};
% cfg_in.data_dir = [data_root 'R112-2017-08-01_post'];
% cfg_in.inter_dir = [data_root 'PPC_figs'];
% cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
% cfg_in.shuffle = nShuffle; 
% cfg_in.plot = 1;
% cfg_in.min_nSpikes = 500; 
% cfg_in.spike_id = 'ST1_cut'; 
% PPC_comparison_CG(cfg_in); 
% disp('complete')
% 
% %%
% disp('112 ST6 ...')
% cfg_in = [];
% cfg_in.dataset = {'CSC12.ncs'};
% cfg_in.data_dir = [data_root '/R112-2017-08-01_post'];
% cfg_in.inter_dir = [data_root '/PPC_figs'];
% cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
% cfg_in.shuffle = nShuffle; 
% cfg_in.plot = 1;
% cfg_in.min_nSpikes = 500; 
% cfg_in.spike_id = 'ST6_cut'; 

%% Complete
disp('112 ST1 ...') %IL
cfg_in = [];
cfg_in.dataset = {'CSC1_2k.ncs'};
cfg_in.data_dir = [data_root 'R112-2017-07-31_post'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST1_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
close all
%% Complete
disp('112 ST2 ...') %PL
cfg_in = [];
cfg_in.dataset = {'CSC3_2k.ncs'};
cfg_in.data_dir = [data_root 'R112-2017-07-31_post'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST2_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
close all
%% Complete
disp('112 TT7 ...') %CG
cfg_in = [];
cfg_in.dataset = {'CSC14_2k.ncs'};
cfg_in.data_dir = [data_root 'R112-2017-07-31_post'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT7_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
close all
%% 
disp('112 ST1 ...') % IL
cfg_in = [];
cfg_in.dataset = {'CSC1.ncs'};
cfg_in.data_dir = [data_root 'R112-2017-08-04_post'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST1_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
close all

%% Complete
disp('107 TT2 ...') %OFC
cfg_in = [];
cfg_in.dataset = {'CSC7_2k.ncs'};
cfg_in.data_dir = [data_root 'R107-2017-07-31_pre'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT2_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
close all
%% Complete
disp('122 TT1 ...') % Piri_OFC
cfg_in = [];
cfg_in.dataset = {'CSC1.ncs'};
cfg_in.data_dir = [data_root 'Spikes_R122-2017-02-26_SPIKES'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT1_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
close all
%% Complete
disp('40262 TT2 ...') % OFC
cfg_in = [];
cfg_in.dataset = {'CSC6_2k.ncs'};
cfg_in.data_dir = [data_root 'R40262-2017-12-02_trk'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT2_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')

%% 
disp('40264 TT2 ...') % OFC
cfg_in = [];
cfg_in.dataset = {'CSC7_2k.ncs'};
cfg_in.data_dir = [data_root 'R40264-2017-12-02_trk'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT2_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
%%
disp('112 ST6 ...') %Piri_NAc*
cfg_in = [];
cfg_in.dataset = {'CSC12.ncs'};
cfg_in.data_dir = [data_root 'R112-2017-08-01_post'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST6_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')

%% complete
disp('112 TT7 ...') %CG
cfg_in = [];
cfg_in.dataset = {'CSC14.ncs'};
cfg_in.data_dir = [data_root 'R112-2017-07-31_post'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'TT7_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')
%% Complete
disp('108 ST3 ...') % OFC
cfg_in = [];
cfg_in.dataset = {'CSC5.ncs'};
cfg_in.data_dir = [data_root 'R108-2017-08-04_post'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.spike_id = 'ST3_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')

%% complete
disp('123 TT1 ...') %Piri_OFC
cfg_in = [];
cfg_in.dataset = {'CSC1.ncs'};
cfg_in.data_dir = [data_root 'R123-2017-02-26_post_trk'];
cfg_in.inter_dir = [data_root 'PPC_figs'];
cfg_in.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_in.shuffle = 100; 
cfg_in.plot = 1;
cfg_in.min_nSpikes = 500; 
cfg_in.tstart = 78*2000; % saturations will yield NaNs.  This trims the first half of the sessions since it contains saturations.   
cfg_in.tstop = (578*2000)-cfg_in.tstart; % same as above
cfg_in.spike_id = 'TT1_SS*.t'; 

PPC_comparison_CG(cfg_in); 
disp('complete')