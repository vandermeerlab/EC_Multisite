%% PPC across cells for Jimmie data

clear all
close all
restoredefaultpath
cd('/Users/jericcarmichael/Documents/Multisite/R123-2017-02-27_pre_pot')

%% prepare for fieldtrip
addpath('/Users/jericcarmichael/Documents/GitHub/fieldtrip')
ft_defaults
addpath('/Users/jericcarmichael/Documents/GitHub/EC_Multisite/Basic_functions/')
%% Load everything
cfg        = [];
cfg.dataset = {'CSC16.ncs'};

data = ft_read_neuralynx_interp({cfg.dataset});

%get the Ts files for all cells. and appended them to the csc file
% t_id = FindFiles('*.t');
% for iS = 1:length(t_id)
t_id = 'TT4_01.t'
    spike = ft_read_spike(t_id); % needs fixed read_mclust_t.m
    data = ft_appendspike([],data, spike);

% end
% %% redefine trials ninto pre, main, post
% evt = LoadEvents([])
% 
% data_trel = ft_redefinetrial(cfg_redef, data)
% 
% if isfield(cfg_in, 'plot')
%     cfg_pliot = []; 
%     cfg_plot.spikechannel = data.label(2:end);
%    ft_spike_plot_raster(cfg_plot, data) 
%     
%     
%     
%     
% end

%% STA
cfg              = [];
cfg.timwin       = [-0.5 0.5]; %
cfg.spikechannel = data.label{2};
cfg.channel      = data.label{1}; 
staAll           = ft_spiketriggeredaverage(cfg, data);
 
% plot
figure
plot(staAll.time, staAll.avg(:,:)');
legend(data.label); h = title(cfg.spikechannel); set(h,'Interpreter','none');
set(gca,'FontSize',14,'XLim',cfg.timwin,'XTick',cfg.timwin(1):0.1:cfg.timwin(2)); 
xlabel('time (s)'); grid on;

%% ppc etc
cfg           = [];
cfg.method    = 'mtmconvol';
cfg.foi       = 5:1:100;
cfg.t_ftimwin = 11./cfg.foi; % cycles per frequency
cfg.taper     = 'hanning';
cfg.spikechannel = data.label{2};
cfg.channel      = data.label{1};
stsConvol     = ft_spiketriggeredspectrum(cfg, data); % note, use raw or interpolated version

% plot
% plot(stsConvol.freq,nanmean(sq(abs(stsConvol.fourierspctrm{1}))))

%%
cfg               = [];
cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
cfg.spikechannel  = stsConvol.label;
cfg.channel       = stsConvol.lfplabel; % selected LFP channels
cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
cfg.timwin        = 'all'; % compute over all available spikes in the window
%cfg.latency       = [-2.5 0]; % sustained visual stimulation period
statSts           = ft_spiketriggeredspectrum_stat(cfg,stsConvol);

% plot the results
figure;
plot(statSts.freq,statSts.ppc0')
set(0,'DefaultTextInterpreter','none');
set(gca,'FontSize',18);
xlabel('frequency')
ylabel('PPC')
title(cfg.spikechannel);