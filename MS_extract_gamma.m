function evts = MS_extract_gamma(cfg_in, data)

%% MS_extractgamma : find and extract all the gamma events in an LFP trace
%
%
%
%          Inputs:
%           - cfg [struct]: cfg.chan_to_use
%           - csc [TSD]:
%          Outputs:
%           - gamma_evts [struct]: contains all the gamma events
%           -
%           -
%
% EC - 2016-05-06
%%


%% Overall
cfg_def = [];
cfg_def.debug = 0;
cfg_def.f_bandpass = {[40 55],[70 100]};
cfg_all = ProcessConfig2(cfg_def, cfg_in);
%% compare OFCs
Phases = {'pre', 'ipsi', 'contra', 'post'};

%% count the gamma event occurance
% Phases = {'control', 'ipsi', 'contra'};
sites = fieldnames(data.pre);
sites = sites(1:end-2); % removes the position data subfield from the list.

for iSite = 1:length(sites)
    Naris.(sites{iSite}).control = UnionTSD([], data.pre.(sites{iSite}), data.post.(sites{iSite}));
    Naris.(sites{iSite}).control.cfg.hdr = data.pre.(sites{iSite}).cfg.hdr;
    Naris.(sites{iSite}).pre = data.pre.(sites{iSite});
    Naris.(sites{iSite}).post = data.post.(sites{iSite});
    Naris.(sites{iSite}).ipsi = data.ipsi.(sites{iSite});
    Naris.(sites{iSite}).contra = data.contra.(sites{iSite});
end
evts = [];
% check union
if cfg_all.debug
    figure(100)
    hold on
    plot(Naris.NAc_pot.control.tvec, Naris.NAc_pot.control.data, '--r')
    plot(data.pre.NAc_pot.tvec, data.pre.NAc_pot.data, '-.b', data.post.NAc_pot.tvec, data.post.NAc_pot.data, '-.k')
end

%% get the gamma events as per the the vStr naris method
exp = {'control','pre', 'ipsi', 'contra', 'post'}; %v 'pre', , 'post'
for iExp = 1:length(exp)
    for iSite = 1:length(sites)
        if strcmp(exp{iExp}, 'control') ==1
            [evts.(sites{iSite}).(exp{iExp}), ~, evt_thr.(sites{iSite}).control] = MS_DetectEvents_thresholds([], Naris.(sites{iSite}).(exp{iExp}), data.pre.ExpKeys);
            cfg= []; cfg.detect_method = 'raw'; cfg.detect_thr = evt_thr.(sites{iSite}).control;
            cfg.f_bandpass = cfg_all.f_bandpass;
            [evts.(sites{iSite}).(exp{iExp}), ~, ~] = MS_DetectEvents(cfg, Naris.(sites{iSite}).(exp{iExp}), data.pre.ExpKeys);
        else
            cfg= []; cfg.detect_method = 'raw'; cfg.detect_thr = evt_thr.(sites{iSite}).control;
            cfg.f_bandpass = cfg_all.f_bandpass;
            [evts.(sites{iSite}).(exp{iExp}), ~, ~] = MS_DetectEvents(cfg, Naris.(sites{iSite}).(exp{iExp}), data.(exp{iExp}).ExpKeys);
        end
        bands = fieldnames(evts.(sites{iSite}).(exp{iExp}));
        for iBand = 1:length(bands)
            if strcmp(exp{iExp}, 'control')
                evts.(sites{iSite}).(exp{iExp}).(bands{iBand}).rate = length(evts.(sites{iSite}).(exp{iExp}).(bands{iBand}).tstart)/...
                    (((Naris.(sites{iSite}).(exp{2}).tvec(end)-Naris.(sites{iSite}).(exp{2}).tvec(1))+(Naris.(sites{iSite}).(exp{5}).tvec(end)-Naris.(sites{iSite}).(exp{5}).tvec(1)))/60);
            else
                evts.(sites{iSite}).(exp{iExp}).(bands{iBand}).rate = length(evts.(sites{iSite}).(exp{iExp}).(bands{iBand}).tstart)/((Naris.(sites{iSite}).(exp{iExp}).tvec(end)-Naris.(sites{iSite}).(exp{iExp}).tvec(1))/60);
            end
        end
    end
end


%% report the outcome
fprintf('\n\n############# EVENTS DETECTED #############')
for iSite = 1:length(sites)
    if isfield(data.pre.ExpKeys, 'ratID')
        fprintf(['\n\n' data.pre.ExpKeys.ratID ' - ' data.pre.ExpKeys.date '\n'])
    else
        fprintf(['\n\n' data.pre.ExpKeys.Subject ' - ' data.pre.ExpKeys.date '\n'])
    end
    disp(['Events:  ' sites{iSite}])
    for iBand = 1:length(bands)
        fprintf(['\n' bands{iBand}(1) '      '])
        for iPhase = 1:length(Phases)
            fprintf([Phases{iPhase}(1:2) ':' num2str(length(evts.(sites{iSite}).(Phases{iPhase}).(bands{iBand}).tstart)) ' ; '])
        end
    end
end
%% get the phase lag for each pair
% pairs = ExpKeys.GoodPairs;
% type = {'low', 'high'};
%
% if isempty(pairs{1})
%     for iPhase = 1:length(Phases)
%         phase_lag_out.(Phases{iPhase}).low = NaN*ones(4,4);
%         phase_lag_out.(Phases{iPhase}).high = NaN*ones(4,4);
%     end
% else
%     % parameters
%     f_bandpass = {[40 55]; [70 85]};
%     cfg = [];
%     cfg.spec_window = 128;
%     cfg.NFFT = 1024;
%     cfg.ref_loc = 'vl'; % which corner of the probe to use as a reference (can be 'vl' (ventrolateral), 'vm' (ventromedial), 'dl' (dorsolateral), 'dm' (dorsal medial)
%     cfg.output = 'rad'; % output in radians ('rad') or degrees ('deg')
%     % compute the phase difference between each channel and the reference using cross-spectral power desnity.
%     % a negative values represent lags, while positives are leads.
%     for iPairs = 1:length(pairs)
%         for itype = 1:length(type)
%             for iPhase = 1:length(Phases)
%                 all_coh_spec.(Phases{iPhase}).(type{itype}) = [];
%                 all_coh_spec_21.(Phases{iPhase}).(type{itype}) = [];
%                 site_12 = {strrep(pairs{iPairs}(1:3), '_', ''), strrep(pairs{iPairs}(end-2:end), '_', '')};
%                 [~, id] = max([length(evts.([site_12{1} '_pot']).(Phases{iPhase}).(type{itype}).tstart), length(evts.([site_12{2} '_pot']).(Phases{iPhase}).(type{itype}).tstart)]); % take the larger of the two gamma events
%                 trials = evts.([site_12{id} '_pot']).(Phases{iPhase}).(type{itype});
%                 tstart_idx = nearest_idx3(trials.tstart,data.(Phases{iPhase}).([site_12{1} '_pot']).tvec);
%                 tend_idx = nearest_idx3(trials.tend,data.(Phases{iPhase}).([site_12{1} '_pot']).tvec);
%
%                 % filter the entire signal
%                 cfg_filter = [];
%                 cfg_filter.f = f_bandpass{itype};
%                 cfg_filter.type = 'cheby1';
%                 cfg_filter.order = 5;
%                 data_in_1 = FilterLFP(cfg_filter, data.(Phases{iPhase}).([site_12{1} '_pot']));
%                 data_in_2 = FilterLFP(cfg_filter, data.(Phases{iPhase}).([site_12{2} '_pot']));
%                 Add_ft()
%
%                 % loop over events and compare the phase of
%                 for itrial = 1:length(trials.tstart)
%                     data_in_cycle = restrict( data_in_1, trials.tstart(itrial), trials.tend(itrial));
%                     cfg_def.data_smooth = 10;
%                     cfg_def.session_type = 'pre';
%                     cfg_def.smooth = 10;
%                     cfg_def.plot = 'off';
%                     cfg_def.MPD = data_in_cycle.cfg.hdr{1}.SamplingFrequency/200;
%                     cfg_def.peak_to_use = 'next';
%                     cycle = AMPX_extract_max_cycle(cfg_def, data_in_cycle);
%                     close all
%                     if isempty(cycle)
%                         continue
%                     end
%                     cycle_idx = tstart_idx(itrial):tend_idx(itrial);
%
%                     [Cxy,F] = cpsd(data_in_1.data(cycle_idx(cycle.prev.idx(1):cycle.next.idx(2))),data_in_2.data(cycle_idx(cycle.prev.idx(1):cycle.next.idx(2))),hamming(cfg.spec_window/2),cfg.spec_window/4,cfg.NFFT,data_in_1.cfg.hdr{1}.SamplingFrequency);
%                     coh_spec_fake = -angle(Cxy); %higher value means leading. outputs radians
%                     all_coh_spec.(Phases{iPhase}).(type{itype})(itrial) =circ_mean(coh_spec_fake(nearest(F,cfg_filter.f(1)):nearest(F,cfg_filter.f(2))));
%
%                     [Cxy,F] = cpsd(data_in_2.data(cycle_idx(cycle.prev.idx(1):cycle.next.idx(2))),data_in_1.data(cycle_idx(cycle.prev.idx(1):cycle.next.idx(2))),hamming(cfg.spec_window/2),cfg.spec_window/4,cfg.NFFT,data_in_1.cfg.hdr{1}.SamplingFrequency);
%                     coh_spec_fake = -angle(Cxy); %higher value means leading. outputs radians
%                     all_coh_spec_21.(Phases{iPhase}).(type{itype})(itrial) = circ_mean(coh_spec_fake(nearest(F,cfg_filter.f(1)):nearest(F,cfg_filter.f(2))));
%
%                 end
%                 all_coh_spec.(Phases{iPhase}).(type{itype})(all_coh_spec.(Phases{iPhase}).(type{itype})==0) =[];
%                 all_coh_spec_21.(Phases{iPhase}).(type{itype})(all_coh_spec_21.(Phases{iPhase}).(type{itype})==0) =[];
%                 Remove_ft()
%                 if strcmp(pairs{iPairs}, 'OFC_NAc')
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(1,2) = rad2deg(circ_mean(all_coh_spec.(Phases{iPhase}).(type{itype})'));
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(2,1) = rad2deg(circ_mean(all_coh_spec_21.(Phases{iPhase}).(type{itype})'));
%                 elseif strcmp(pairs{iPairs}, 'OFC_CG')
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(1,3) = rad2deg(circ_mean(all_coh_spec.(Phases{iPhase}).(type{itype})'));
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(3,1) = rad2deg(circ_mean(all_coh_spec_21.(Phases{iPhase}).(type{itype})'));
%                 elseif strcmp(pairs{iPairs}, 'PL_OFC')
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(1,4) = rad2deg(circ_mean(all_coh_spec.(Phases{iPhase}).(type{itype})'));
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(4,1) = rad2deg(circ_mean(all_coh_spec_21.(Phases{iPhase}).(type{itype})'));
%                 elseif strcmp(pairs{iPairs}, 'NAc_Cg')
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(2,3) = rad2deg(circ_mean(all_coh_spec.(Phases{iPhase}).(type{itype})'));
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(3,2) = rad2deg(circ_mean(all_coh_spec_21.(Phases{iPhase}).(type{itype})'));
%                 elseif strcmp(pairs{iPairs}, 'PL_NAc')
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(2,3) = rad2deg(circ_mean(all_coh_spec.(Phases{iPhase}).(type{itype})'));
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(3,2) = rad2deg(circ_mean(all_coh_spec_21.(Phases{iPhase}).(type{itype})'));
%                 elseif strcmp(pairs{iPairs}, 'PL_CG')
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(3,4) = rad2deg(circ_mean(all_coh_spec.(Phases{iPhase}).(type{itype})'));
%                     phase_lag_out.(Phases{iPhase}).(type{itype})(4,3) = rad2deg(circ_mean(all_coh_spec_21.(Phases{iPhase}).(type{itype})'));
%                 end
%             end
%         end
%     end
% end
%
% %% make a plot of the gamma events over time of recording for each phase/site
% count = [];
% exp = {'control', 'ipsi', 'contra'};
% type = {'low', 'high'};
% for iExp = 1:length(exp)
%     for iSite = 1:length(sites)
%         if strcmp(exp{iExp}, 'control')
%             if isempty(evts.(sites{iSite}).pre.(type{1}).tstart)
%                 pre = 0;
%             else
%                 pre = length(evts.(sites{iSite}).pre.(type{1}).tstart)/((length(Naris.(sites{iSite}).pre.data)/Naris.(sites{iSite}).pre.cfg.hdr{1}.SamplingFrequency)/60);
%             end
%             post = length(evts.(sites{iSite}).post.(type{1}).tstart)/((length(Naris.(sites{iSite}).post.data)/Naris.(sites{iSite}).post.cfg.hdr{1}.SamplingFrequency)/60);
%             count_L(iSite, iExp) = median([pre, post]);
%             if isempty(evts.(sites{iSite}).pre.(type{2}).tstart)
%                 pre = 0;
%             else
%                 pre = length(evts.(sites{iSite}).pre.(type{2}).tstart)/((length(Naris.(sites{iSite}).pre.data)/Naris.(sites{iSite}).pre.cfg.hdr{1}.SamplingFrequency)/60);
%             end
%             post = length(evts.(sites{iSite}).post.(type{2}).tstart)/((length(Naris.(sites{iSite}).post.data)/Naris.(sites{iSite}).post.cfg.hdr{1}.SamplingFrequency)/60);
%             count_H(iSite, iExp) = median([pre, post]);
%         else
%             count_L(iSite, iExp) = length(evts.(sites{iSite}).(exp{iExp}).(type{1}).tstart)/((length(Naris.(sites{iSite}).(exp{iExp}).data)/Naris.(sites{iSite}).(exp{iExp}).cfg.hdr{1}.SamplingFrequency)/60);
%             count_H(iSite, iExp) = length(evts.(sites{iSite}).(exp{iExp}).(type{2}).tstart)/((length(Naris.(sites{iSite}).(exp{iExp}).data)/Naris.(sites{iSite}).(exp{iExp}).cfg.hdr{1}.SamplingFrequency)/60);
%         end
%         labels{iSite,1} = sites{iSite};
%
%     end
% end
% % out = [label'; num2cell(count)]
% disp(fname)
% % disp(count)
% disp(labels)
% load('G:\JK_recordings\Naris\Naris_count_out_new_f2.mat')
% if strcmp(fname(1:3), 'R10')
%     out.(fname(1:4)).(['D_' strrep(fname(end-9:end), '-', '_')]).low = count_L;
%     out.(fname(1:4)).(['D_' strrep(fname(end-9:end), '-', '_')]).high = count_H;
%     out.(fname(1:4)).(['D_' strrep(fname(end-9:end), '-', '_')]).phase_lag = phase_lag_out;
% else
%     out.(fname(1:7)).(['D_' strrep(fname(end-9:end), '-', '_')]).low = count_L;
%     out.(fname(1:7)).(['D_' strrep(fname(end-9:end), '-', '_')]).high = count_H;
%     out.(fname(1:7)).(['D_' strrep(fname(end-9:end), '-', '_')]).phase_lag = phase_lag_out;
% end
% save('G:\JK_recordings\Naris\Naris_count_out_new_f2.mat', 'out')
