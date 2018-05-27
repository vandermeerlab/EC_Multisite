function PPC_comparison(cfg_in)
%% PPC_comparison

%%
cfg_def = [];
cfg_def.dataset = {'CSC8.ncs','CSC24.ncs','CSC30.ncs'};
cfg_def.data_dir = '/Users/jericcarmichael/Documents/R111-2017-06-20-Rec_auto for Eroc';
cfg_def.inter_dir = '/Users/jericcarmichael/Documents/R111-2017-06-20-Rec_auto for Eroc';
cfg_def.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_def.shuffle = 10;
cfg_def.min_nSpikes = 500;
cfg = ProcessConfig(cfg_def, cfg_in);

cd(cfg.data_dir)

data = ft_read_neuralynx_interp(cfg.dataset);

LFP_list = cfg.dataset;

%get the Ts files for all cells. and appended them to the csc file
t_id = FindFiles('*.t');
S_list = {};
for iS = 1:length(t_id)
    spike = ft_read_spike(t_id{iS}); % needs fixed read_mclust_t.m
    disp([spike.label{1} 'Contained: ' num2str(length(spike.timestamp{1})) ' spikes'])
    if length(spike.timestamp{1}) < cfg.min_nSpikes
        continue
    end
    S_list{end+1} = spike.label{1};
    data = ft_appendspike([],data, spike);
end
%% redefine trials as pre, task, post
evt = LoadEvents([]);
cfg_csc.fc = LFP_list(1);
d_temp = LoadCSC(cfg_csc); % get the data in the TSD format to get the trial index

idx = strfind(evt.label, 'Starting Recording');
start_idx = find(not(cellfun('isempty', idx)));
idx = strfind(evt.label, 'Stopping Recording');
stop_idx = find(not(cellfun('isempty', idx)));

tstart = nearest_idx(evt.t{start_idx}(cfg.phase), d_temp.tvec);
tstop = nearest_idx(evt.t{stop_idx}(cfg.phase), d_temp.tvec);

cfg_trl = [];
cfg_trl.begsample = tstart;
cfg_trl.endsample = tstop;
data_trl = ft_redefinetrial(cfg_trl, data);

if isfield(cfg, 'plot')
    figure(99)
    plot(data_trl.time{1}(1:2000), data_trl.trial{1}(1:length(LFP_list),1:2000))
    
end

%%
for  iLFP = 1:length(LFP_list)
    for iS = 1:length(S_list)
        spk_chan = S_list{iS};
        lfp_chan = LFP_list{iLFP};
        
        cfg_i              = [];
        cfg_i.timwin       = [-0.002 0.006]; % remove 4 ms around every spike
        cfg_i.spikechannel = spk_chan;
        cfg_i.channel      = lfp_chan(1:end-4);
        cfg_i.method       = 'linear'; % remove the replaced segment with interpolation
        data_i        = ft_spiketriggeredinterpolation(cfg_i, data_trl);
        
        %% STA
        cfg_sta              = [];
        cfg_sta.timwin       = [-0.5 0.5]; %
        cfg_sta.spikechannel = spk_chan;
        cfg_sta.channel      = lfp_chan(1:end-4);
        staAll           = ft_spiketriggeredaverage(cfg_sta, data_i);
        
        % plot
        if isfield(cfg, 'plot')
            
            figure
            plot(staAll.time, staAll.avg(:,:)');
            legend(data.label); h = title(cfg_sta.spikechannel); set(h,'Interpreter','none');
            set(gca,'FontSize',14,'XLim',cfg_sta.timwin,'XTick',cfg_sta.timwin(1):0.1:cfg_sta.timwin(2));
            xlabel('time (s)'); grid on;
        end
        %% ppc etc
        cfg_ppc            = [];
        cfg_ppc.method    = 'mtmconvol';
        cfg_ppc.foi       = 1:1:100;
        cfg_ppc.t_ftimwin = 5./cfg_ppc.foi; % cycles per frequency
        cfg_ppc.taper     = 'hanning';
        cfg_ppc.spikechannel = spk_chan;
        cfg_ppc.channel      = lfp_chan(1:end-4);
        stsConvol     = ft_spiketriggeredspectrum(cfg_ppc, data_i); % note, use raw or interpolated version
        
        % plot
        if isfield(cfg, 'plot')
            
            plot(stsConvol.freq,nanmean(sq(abs(stsConvol.fourierspctrm{1}))))
        end
        %%
        cfg_ppc                = [];
        cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
        cfg_ppc.spikechannel = spk_chan;
        cfg_ppc.channel      = lfp_chan(1:end-4);
        %cfg.dojack = 1;
        cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
        cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
        %cfg.latency       = [-2.5 0]; % sustained visual stimulation period
        statSts           = ft_spiketriggeredspectrum_stat(cfg_ppc,stsConvol);
        
        % plot the results
        figure;
        plot(statSts.freq,statSts.ppc0')
        set(0,'DefaultTextInterpreter','none');
        set(gca,'FontSize',18);
        xlabel('frequency')
        ylabel('PPC')
        title(cfg_ppc.spikechannel);
        
        obs_freq = statSts.freq;
        obs_ppc = statSts.ppc0';
        
        
        %%
        nShuf = cfg.shuffle;
        iChan = length(data_i.label)+1;
        data_i.label{iChan} = 'temp_shuf';
        %%
        t_idx = strfind(data_i.label, S_list{iS});
        spk_idx = find(not(cellfun('isempty', t_idx)));
        shuf_ppc = zeros(nShuf,length(obs_freq));
        parfor iShuf = 1:nShuf
            fprintf('Shuffle %d...\n',iShuf);
           shuf_ppc(iShuf,:) =  Shuffle_PPC(data_i, spk_idx, lfp_chan, iChan);
        end

        id = strrep(spk_chan, '-', '_');
        
        %% plot
        close all
        if isfield(cfg, 'plot')
            figure(111)
            hold on;
            h(1) = plot(obs_freq,obs_ppc,'k','LineWidth',2);
            plot(obs_freq,nanmean(shuf_ppc,1),'r');
            plot(obs_freq,nanmean(shuf_ppc,1)+nanstd(shuf_ppc,1),'r:');
            plot(obs_freq,nanmean(shuf_ppc,1)-nanstd(shuf_ppc,1),'r:');
            
            
            set(0,'DefaultTextInterpreter','none');
            legend(h,{'observed','shuffled'},'Location','Northeast'); legend boxoff;
            set(gca,'FontSize',18);
            xlabel('frequency')
            ylabel('PPC')
            title([spk_chan '_' lfp_chan]);
            
            mkdir(cfg.inter_dir, 'PPC')
            if isunix
                saveas(gcf, [cfg.inter_dir '/PPC/' id '_csc_' lfp_chan(1:end-4)], 'fig');
                saveas(gcf, [cfg.inter_dir '/PPC/' id '_csc_' lfp_chan(1:end-4)], 'png');
            else
                saveas(gcf, [cfg.inter_dir '\PPC\' id '_csc_' lfp_chan(1:end-4)], 'fig');
                saveas(gcf, [cfg.inter_dir '\PPC\' id '_csc_' lfp_chan(1:end-4)], 'png');
            end
        end
        %% collect variables for export
        PPC.(lfp_chan(1:end-4)).(id).obs_freq = obs_freq;
        PPC.(lfp_chan(1:end-4)).(id).obs_ppc = obs_ppc;
        PPC.(lfp_chan(1:end-4)).(id).shuf_freq = shuf_ppc;
        PPC.(lfp_chan(1:end-4)).(id).staAll = staAll;
        
    end
end
[~,dir_id] = fileparts(pwd);
dir_id = strrep(dir_id, '-', '_');
dir_id = strrep(dir_id, ' ', '_');
save(['PPC_' dir_id], 'PPC',  '-v7.3');
end

%%
function shuf_ppc = Shuffle_PPC(data_i, spk_idx, lfp_chan, iChan)
% shuffle once
for iT = 1:length(data_i.trial) % shuffle each trial separately
    orig_data = data_i.trial{iT}(spk_idx,:);
    data_i.trial{iT}(iChan,:) = orig_data(randperm(length(orig_data)));

end

%% ppc etc
cfg_ppc            = [];
cfg_ppc.method    = 'mtmconvol';
cfg_ppc.foi       = 1:1:100;
cfg_ppc.t_ftimwin = 5./cfg_ppc.foi; % cycles per frequency
cfg_ppc.taper     = 'hanning';
cfg_ppc.spikechannel = 'temp_shuf';
cfg_ppc.channel      = lfp_chan(1:end-4);
stsConvol     = ft_spiketriggeredspectrum(cfg_ppc , data_i); % note, use raw or interpolated version

% plot
%plot(stsConvol.freq,nanmean(sq(abs(stsConvol.fourierspctrm{1}))))

%%
cfg_ppc               = [];
cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
cfg_ppc.spikechannel = 'temp_shuf';
cfg_ppc.channel      = lfp_chan(1:end-4);
%cfg.dojack = 1;
cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
%cfg.latency       = [-2.5 0]; % sustained visual stimulation period
statSts           = ft_spiketriggeredspectrum_stat(cfg_ppc ,stsConvol);

shuf_ppc = statSts.ppc0';


end % of shuffles

