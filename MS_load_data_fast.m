function [data, cfg] = MS_load_data(cfg_in)
%% MS_load_data:
%  used to load the raw LFP recodings across multiple subject and segragate
%  the data into the different phases ('pre', 'ipsi', 'contra', 'post') as
%  well as pot vs track segments
%
%  inputs:
%      global PARAMS from MASTER_Multisite
%    - cfg_in: input configuration
%
%
%  outputs:
%    - data [struct]: contains data for each subjects separated into
%    session(day), phases(pre, ipsi, contra, post), segments(pot, trk)
%    -cfg_out [struct]: contains the configurations use to process the data

%% Initial paramters
global PARAMS

mfun = mfilename;
cfg_def = [];

cfg= ProcessConfig2(cfg_def, cfg_in);


%% load the data fpr each phase within the session.
if isunix
    cd([PARAMS.data_dir '/' cfg.fname(1:4) '/' cfg.fname ])
else
    cd([PARAMS.data_dir '\' cfg.fname(1:4) '\' cfg.fname ])
end

ExpKeys = MS_Load_Keys();
%%
for iPhase = 1:length(PARAMS.Phases)
    if isunix
        cd([PARAMS.data_dir '/' cfg.fname(1:4) '/' cfg.fname '/' cfg.fname '_' PARAMS.Phases{iPhase}])
    else
        cd([PARAMS.data_dir '\' cfg.fname(1:4) '\' cfg.fname '\' cfg.fname '_' PARAMS.Phases{iPhase}])
    end
    
    evt= LoadEvents([]);
    % check to see if I was dumb and didnt stop recording between the pot
    % and track.  If so, then split the session in half shave off 1/10 of
    % off the end for the first half and the start of the second.
    idx = strfind(evt.label, 'Starting Recording');
    start_idx = find(not(cellfun('isempty', idx)));
    idx = strfind(evt.label, 'Stopping Recording');
    stop_idx = find(not(cellfun('isempty', idx)));
    
    % check to make sure NLX didn't mess up the start stop events (this
    % happened for R104-2016-09-26_ipsi. It has 19 start times for some
    % reason)
    if length(evt.t{start_idx}) >= 3 % should only have 2
        [~, trk_idx]  = max(evt.t{stop_idx}(2:end)-evt.t{start_idx}(2:end)); % find the largest gap
        trk_idx = trk_idx+1; % offset by one to compensate for the 'trk' being the second phase
    else
        trk_idx = [];
    end
    
    if evt.t{stop_idx}(1)-evt.t{start_idx}(1) < 60*9; % check to ensure the recording is roughly 10mins long.
        error('pot session is too short');
    elseif evt.t{stop_idx}(trk_idx)-evt.t{start_idx}(trk_idx) < 60*9;
        error('trk session is too short');
    end
    
%% load all the channels of interest
    for iChan = 1:length(ExpKeys.Chan_to_use)
        cfg_load.fc = ExpKeys.Chan_to_use(iChan);
        cfg_load.resample = 2000;
        csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)) = LoadCSC(cfg_load);
        
        fprintf(['\n' ExpKeys.Chan_to_use_labels{iChan}(1:end-1) '_FS:' num2str(csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).cfg.hdr{1}.SamplingFrequency) '\n'])
        % check to see if the data has been sampled appropriotely as per
        % cfg_load.resample
        if csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).cfg.hdr{1}.SamplingFrequency > cfg_load.resample
            cfg.decimateByFactor = csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).cfg.hdr{1}.SamplingFrequency/cfg_load.resample;
            
            fprintf('%s: Decimating by factor %d...\n',mfun,cfg.decimateByFactor)
            csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).data = decimate(csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).data,cfg.decimateByFactor);
            csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).tvec = csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).tvec(1:cfg.decimateByFactor:end);
            csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).cfg.hdr{1}.SamplingFrequency = csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).cfg.hdr{1}.SamplingFrequency./cfg.decimateByFactor;
            
            if round(mode(diff(csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).tvec))*10000) ~= floor((1/(csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)).cfg.hdr{1}.SamplingFrequency)*10000));
                error('Something went wrong.  The diff in tvec samples does not match the sampling frequency following resampling')
            end
        end
        
        % split into pot and track.
        data.(PARAMS.Phases{iPhase}).([ExpKeys.Chan_to_use_labels{iChan}(1:end-1) '_pot']) = restrict(csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)),evt.t{start_idx}(1),evt.t{stop_idx}(1));
        if isempty(trk_idx)
            data.(PARAMS.Phases{iPhase}).([ExpKeys.Chan_to_use_labels{iChan}(1:end-1) '_trk']) = restrict(csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)),evt.t{start_idx}(2),evt.t{stop_idx}(2));
        else
            data.(PARAMS.Phases{iPhase}).([ExpKeys.Chan_to_use_labels{iChan}(1:end-1) '_trk']) = restrict(csc_out.(ExpKeys.Chan_to_use_labels{iChan}(1:end-1)),evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
        end
    end
    % get the position for each phase
    if isempty(FindFiles('*.nvt'))
        unzip('VT1.zip')
    end
    cfg_load = [];
    pos = LoadPos(cfg_load);
    data.(PARAMS.Phases{iPhase}).pos.pot = restrict(pos,evt.t{start_idx}(1),evt.t{stop_idx}(1));
    data.(PARAMS.Phases{iPhase}).pos.trk = restrict(pos,evt.t{start_idx}(2),evt.t{stop_idx}(2));
    
    data.(PARAMS.Phases{iPhase}).ExpKeys = ExpKeys;
    
end
