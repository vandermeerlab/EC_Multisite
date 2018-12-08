function PPC_comparison(cfg_in)
%% PPC_comparison
% clear all
% close all
% restoredefaultpath
% cfg_in = [];
% %
% if isunix
%     addpath(genpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'))
%     %
%     addpath('/Users/jericcarmichael/Documents/GitHub/fieldtrip')
%     ft_defaults
%     addpath('/Users/jericcarmichael/Documents/GitHub/EC_Multisite/Basic_functions/')
% else
%     addpath(genpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'))
%     %
%     addpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip')
%     ft_defaults
%     addpath('D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite\Basic_functions\')
% end




%%
cfg_def = [];
cfg_def.dataset = {'CSC7.ncs'};  % [testing] Use CSC8 for MClust and CSC7 for SS3D
if isunix
cfg_def.data_dir = '/Users/jericcarmichael/Documents/Multisite/R112-2017-08-01_post';
else
cfg_def.data_dir = 'G:\Multisite\temp\R112-2017-08-04_post';    % SS3D debug
% cfg_def.data_dir = 'G:\R111-2017-06-20-Rec_auto for Eroc'; % MClust data debug
end
cfg_def.inter_dir = cfg_def.data_dir; % can be changed to some other place
cfg_def.phase = 1; % corresponds to the first recording phase ('pre').  2 = 'task', 3 = 'post'.
cfg_def.shuffle = 10;
cfg_def.min_nSpikes = 300;
cfg_def.spike_id = '.t';
cfg_def.check = 1;

cfg = ProcessConfig(cfg_def, cfg_in);
cd(cfg.data_dir)

data = ft_read_neuralynx_interp(cfg.dataset);
% resample if LFP Fs is greater than 2k
if data.fsample > 2000
    cfg_resample = [];
    cfg_resample.resamplefs = 2000;
    cfg_resample.sampleindex = 'yes';
    data = ft_resampledata(cfg_resample, data); %
end
LFP_list = cfg.dataset;

%get the Ts files for all cells. and appended them to the csc file
t_id = FindFile_str(cd, cfg.spike_id);
% hack around above finding '.txt files
keep_idx = [];
for ii =1:length(t_id)
    if strcmp(t_id{ii}(end-1:end), '.t')
        keep_idx(ii) = 1;
    else
        keep_idx(ii) = 0;
    end
end
t_id(keep_idx==0) = [];
t_id
%% load all the .t files and append them to the data struct
S_list = {};
for iS = 1:length(t_id)
    spike = ft_read_spike(t_id{iS}); % needs fixed read_mclust_t.m
    spike.label{1} = t_id{iS}(1:end-2);
    disp([spike.label{1} ' Contained: ' num2str(length(spike.timestamp{1})) ' spikes'])
    if length(spike.timestamp{1}) < cfg.min_nSpikes
        continue
    end
    S_list{end+1} = spike.label{1};
    data = ft_appendspike([],data, spike);
end

% %% TEMP visualization
% if cfg.check
%    figure(100) 
%    plot(data.time{1}, data.trial{1}(1,:), 'b')
%    hold on
%    iC = 2; % spike index
%    nSpikes = length(data.trial{1}(iC,:));
%    xvals = [ data.trial{1}(iC,:)' ; data.trial{1}(iC,:)' ; nan(size(data.trial{1}(iC,:)')) ];
%    yvals = [ iC.*ones(1,nSpikes);...
%        iC.*ones(1,nSpikes) ; nan(1,nSpikes) ];
%    
%    xvals = xvals(:);
%    yvals = yvals(:);
%    
%    plot(xvals,yvals,'Color','k','LineWidth',2);
%    
%    plot(data.time{1}, data.trial{1}(2,:), '.k')
% %    
% %    for ii = 1:10000%length(data.trial{1}):-1:1
% %    line([data.trial{1}(2,ii), data.trial{1}(2,ii)], [-800 -700])
% %    end
% %     
%     
%     
% end
%% redefine trials as pre, task, post
evt = LoadEvents([]);
cfg_csc.fc = LFP_list(1);
d_temp = LoadCSC(cfg_csc); % get the data in the TSD format to get the trial index

idx = strfind(evt.label, 'Starting Recording');
start_idx = find(not(cellfun('isempty', idx)));
idx = strfind(evt.label, 'Stopping Recording');
stop_idx = find(not(cellfun('isempty', idx)));


% usefule for dealing with saturations in the LFP which kill PPC
if isfield(cfg, 'tstart')
    tstart = cfg.tstart;
else
    tstart = nearest_idx(evt.t{start_idx}(cfg.phase), d_temp.tvec);
end
% usefule for dealing with saturations in the LFP which kill PPC
if isfield(cfg, 'tstop')
    tstop = cfg.tstop;
else
    tstop = nearest_idx(evt.t{stop_idx}(cfg.phase), d_temp.tvec);
end



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
            
            %             mkdir(cfg.inter_dir, 'PPC')
            if isunix
                sess_id = strsplit(pwd, '/');
            else
                sess_id = strsplit(pwd, '\');
            end
            sess_id = strrep(sess_id{end}, '-', '_');
            if isunix
                saveas(gcf, [cfg.inter_dir '/' sess_id '_' id '_' lfp_chan(1:end-4)], 'fig');
                saveas(gcf, [cfg.inter_dir '/' sess_id '_' id '_' lfp_chan(1:end-4)], 'png');
                saveas_eps([ sess_id '_' id '_' lfp_chan(1:end-4)], cfg.inter_dir)
            else
                saveas(gcf, [cfg.inter_dir '\' sess_id '_' id '_' lfp_chan(1:end-4)], 'fig');
                saveas(gcf, [cfg.inter_dir '\' sess_id '_' id '_' lfp_chan(1:end-4)], 'png');
                saveas_eps([ sess_id '_' id '_' lfp_chan(1:end-4)], cfg.inter_dir)
            end
        end
        %% collect variables for export
        PPC.(lfp_chan(1:end-4)).(id).obs_freq = obs_freq;
        PPC.(lfp_chan(1:end-4)).(id).obs_ppc = obs_ppc;
        PPC.(lfp_chan(1:end-4)).(id).shuf_freq = shuf_ppc;
        PPC.(lfp_chan(1:end-4)).(id).staAll = staAll;
        
    end
end
if ~isempty(S_list)
    [~,dir_id] = fileparts(pwd);
    dir_id = strrep(dir_id, '-', '_');
    dir_id = strrep(dir_id, ' ', '_');
    save(['PPC_MS' dir_id], 'PPC',  '-v7.3');
end
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

