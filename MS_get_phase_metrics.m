function [mat_out] = MS_get_phase_metrics(cfg_in, Events, Data);
%% MS_get_phase_metrics: extracts the phase measurements.  This can be the
%       phase coherence, amplitude cross correlation, or the phase slope.
%
%
% INPUTS
% - cfg:   [struct] contains configuration parameters.
% - data1: [N x 1]
%
% EC 2018
cfg_def = [];
cfg_def.resize = [];
cfg_def.resize.d = [-0.1 0.1];


cfg = ProcessConfig2(cfg_def, cfg_in);

global PARAMS
bands = {'low', 'high'};
pairs = Data.pre.ExpKeys.GoodPairs;
%% set up the output matrix
labels = Data.pre.ExpKeys.Chan_to_use_labels;
for ii = 1:length(labels)
    for jj = 1:length(labels)
        mat_out.labels{ii,jj} = [labels{ii}(1:end-1) '_' labels{jj}(1:end-1)];
    end
end
for iPhase = 1:4
    for iBand = 1:2
        %         mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).sess_coh = cell(size(mat_out.labels));
        %         mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_coh = cell(size(mat_out.labels));
        %         mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).sess_amp = cell(size(mat_out.labels));
        %         mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_amp_ac = cell(size(mat_out.labels));
        %         mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_amp_lag = cell(size(mat_out.labels));
        mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).PS = cell(size(mat_out.labels));
    end
end

%% cycle through events
type = 'pot';

for iPair = 1:length(pairs)
    [S_out, ~] = strsplit(pairs{iPair}, '_');
    for Round = 1:2
        if Round ==1
            S1 = S_out{1}; % good site 1
            S2 = S_out{2}; % good site 2
        elseif Round ==2
            S1 = S_out{2}; % good site 1
            S2 = S_out{1}; % good site 2
        end
        % get the index of the current pair of sites
        idx = strfind(mat_out.labels, [S1 '_' S2]);
        [x_idx,y_idx] = find(not(cellfun('isempty', idx)));
        for iPhase = 1:length(PARAMS.Phases)
            for iBand = 1:length(bands)
                evts = Events.([S1 '_' type]).(PARAMS.Phases{iPhase}).(bands{iBand});
                this_data_1 = Data.(PARAMS.Phases{iPhase}).([S1 '_' type]);
                this_data_2 = Data.(PARAMS.Phases{iPhase}).([S2 '_' type]);
                phase_slopes = {};
                % cycle through all events
                for iEvent = length(evts.tstart):-1:1
                    this_event = SelectIV([],evts,iEvent);
                    this_event = ResizeIV(cfg.resize,this_event);
                    
                    d1 = restrict(this_data_1,this_event);
                    d2 = restrict(this_data_2,this_event);
                    %try the function
                    cfg_phase = [];
                    cfg_phase.Fs = d1.cfg.hdr{1}.SamplingFrequency;
                    cfg_phase.debug = 0;
                    phase_slopes{iEvent} = MS_phase_slope(cfg_phase, d1.data, d2.data);
                    close all
                end
                mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).PS{x_idx, y_idx} = phase_slopes;
            end
        end
    end
end

end
