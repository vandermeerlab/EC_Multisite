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


cfg_def = [];

cfg= ProcessConfig2(cfg_def, cfg_in);


%% load the data fpr each phase within the session.
for iPhase = 1:length(PARAMS.Phases)
    if isunix
        cd([PARAMS.data_dir '/' cfg.fname(1:4) '/' cfg.fname '/' cfg.fname '_' PARAMS.Phases{iPhase}])
    else
        cd([PARAMS.data_dir '\' cfg.fname(1:4) '\' cfg.fname '\' cfg.fname '_' PARAMS.Phases{iPhase}])
    end
    LoadExpKeys()
    
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
        [~, trk_idx]  = max(diff(evt.t{start_idx}(2:end))); % find the largest gap
        trk_idx = trk_idx +1; % offset by one to compensate for the "diff"
    else
        trk_idx = [];
    end
    
    % subjects with only OFC, NAc, CG.
    if strcmp(ExpKeys.ratID, 'R547606') || strcmp(ExpKeys.ratID, 'R547576') || strcmp(ExpKeys.ratID, 'R547574')
        cfg_load = [];
        cfg_load.fc = ExpKeys.Chan_to_use(1);
        OFC_csc = LoadCSC(cfg_load);
        OFC_csc.data_label = ExpKeys.Chan_to_use_labels{1};
        cfg_load.fc = ExpKeys.Chan_to_use(2);
        NAc_csc = LoadCSC(cfg_load);
        NAc_csc.data_label = ExpKeys.Chan_to_use_labels{2};
        cfg_load.fc = ExpKeys.Chan_to_use(3);
        CG_csc = LoadCSC(cfg_load);
        CG_csc.data_label = ExpKeys.Chan_to_use_labels{3};
        % collect the data
        %         data.(PARAMS.Phases{iPhase}).OFC = OFC_csc;
        %         data.(PARAMS.Phases{iPhase}).NAc = NAc_csc;
        %         data.(PARAMS.Phases{iPhase}).CG = CG_csc;
        
        % split into pot and track.
        data.(PARAMS.Phases{iPhase}).OFC_pot = restrict(OFC_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).NAc_pot = restrict(NAc_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).CG_pot = restrict(CG_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        
        
        collect the data
        data.(PARAMS.Phases{iPhase}).OFC_trk = restrict(OFC_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
        data.(PARAMS.Phases{iPhase}).NAc_trk = restrict(NAc_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
        data.(PARAMS.Phases{iPhase}).CG_trk = restrict(CG_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
        
        % subjects with PL, OFC, NAc, CG
    elseif strcmp(ExpKeys.ratID, 'R547607') || strcmp(ExpKeys.ratID, 'R547577') || strcmp(ExpKeys.ratID, 'R102') || strcmp(ExpKeys.ratID, 'R104')
        cfg_load = [];
        % get the PL
        cfg_load.fc = ExpKeys.Chan_to_use(1);
        PL_csc = LoadCSC(cfg_load);
        PL_csc.data_label = ExpKeys.Chan_to_use_labels{1};
        % get the OFC
        cfg_load.fc = ExpKeys.Chan_to_use(2);
        OFC_csc = LoadCSC(cfg_load);
        OFC_csc.data_label = ExpKeys.Chan_to_use_labels{2};
        % get the NAc
        cfg_load.fc = ExpKeys.Chan_to_use(3);
        NAc_csc = LoadCSC(cfg_load);
        NAc_csc.data_label = ExpKeys.Chan_to_use_labels{3};
        %
        cfg_load.fc = ExpKeys.Chan_to_use(4);
        CG_csc = LoadCSC(cfg_load);
        CG_csc.data_label = ExpKeys.Chan_to_use_labels{4};
        
        
        % split into pot and track.
        data.(PARAMS.Phases{iPhase}).OFC_pot = restrict(OFC_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).NAc_pot = restrict(NAc_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).CG_pot = restrict(CG_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).PL_pot = restrict(PL_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        
        
        % collect the data
        if isempty(trk_idx)
            data.(PARAMS.Phases{iPhase}).OFC_trk = restrict(OFC_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
            data.(PARAMS.Phases{iPhase}).NAc_trk = restrict(NAc_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
            data.(PARAMS.Phases{iPhase}).CG_trk = restrict(CG_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
            data.(PARAMS.Phases{iPhase}).PL_trk = restrict(PL_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
        else
            data.(PARAMS.Phases{iPhase}).OFC_trk = restrict(OFC_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            data.(PARAMS.Phases{iPhase}).NAc_trk = restrict(NAc_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            data.(PARAMS.Phases{iPhase}).CG_trk = restrict(CG_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            data.(PARAMS.Phases{iPhase}).PL_trk = restrict(PL_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            
        end
        % PL, NAc, CG, and piriform cortex (PC).
    elseif strcmp(ExpKeys.ratID, 'R102')
        cfg_load = [];
        % get the PL
        cfg_load.fc = ExpKeys.Chan_to_use(1);
        PL_csc = LoadCSC(cfg_load);
        PL_csc.data_label = ExpKeys.Chan_to_use_labels{1};
        % get the OFC
        cfg_load.fc = ExpKeys.Chan_to_use(2);
        PC_csc = LoadCSC(cfg_load);
        PC_csc.data_label = ExpKeys.Chan_to_use_labels{2};
        % get the NAc
        cfg_load.fc = ExpKeys.Chan_to_use(3);
        NAc_csc = LoadCSC(cfg_load);
        NAc_csc.data_label = ExpKeys.Chan_to_use_labels{3};
        %
        cfg_load.fc = ExpKeys.Chan_to_use(4);
        CG_csc = LoadCSC(cfg_load);
        CG_csc.data_label = ExpKeys.Chan_to_use_labels{4};
        
        
        % split into pot and track.
        data.(PARAMS.Phases{iPhase}).PC_pot = restrict(PC_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).NAc_pot = restrict(NAc_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).CG_pot = restrict(CG_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        data.(PARAMS.Phases{iPhase}).PL_pot = restrict(PL_csc,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        
        
        % collect the data
        if isempty(trk_idx)
            data.(PARAMS.Phases{iPhase}).PC_trk = restrict(PC_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
            data.(PARAMS.Phases{iPhase}).NAc_trk = restrict(NAc_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
            data.(PARAMS.Phases{iPhase}).CG_trk = restrict(CG_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
            data.(PARAMS.Phases{iPhase}).PL_trk = restrict(PL_csc,evt.t{start_idx}(2),evt.t{stop_idx}(2));
        else
            data.(PARAMS.Phases{iPhase}).PC_trk = restrict(PC_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            data.(PARAMS.Phases{iPhase}).NAc_trk = restrict(NAc_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            data.(PARAMS.Phases{iPhase}).CG_trk = restrict(CG_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            data.(PARAMS.Phases{iPhase}).PL_trk = restrict(PL_csc,evt.t{start_idx}(trk_idx),evt.t{stop_idx}(trk_idx));
            
        end
    end
    % get the position for each phase
    cfg_load = [];
    pos = LoadPos(cfg_load);
    data.(PARAMS.Phases{iPhase}).pos.pot = restrict(pos,evt.t{start_idx}(1),evt.t{stop_idx}(1));
    data.(PARAMS.Phases{iPhase}).pos.trk = restrict(pos,evt.t{start_idx}(2),evt.t{stop_idx}(2));
    
end
