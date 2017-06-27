function [evts_out] = MS_event_pairs(cfg_in, evts_in, data)
%% MS_event_pairs: uses the detected events for each channel and extracts the IVs for the corresponding pair sites
%
%
%
%
%
%
%% default parameters


cfg_def = [];


cfg = ProcessConfig2(cfg_def, cfg_in);

global PARAMS

%% cycle through phases and recording cite pairs (uses only the "good pairs" specified in the ExpKeys)
pairs = data.pre.ExpKeys.GoodPairs;
types = {'_pot', '_trk'};
evts_out = evts_in;
phases = PARAMS.Phases;
bands = {'low', 'high'};

for iType = 1:length(types)
    for iPair = 1:length(pairs)
        sites = strsplit(pairs{iPair}, '_');
        if length(sites) >2
            newsites{1} = [sites{1} sites{2}];
            newsites{2} = sites{3};
            sites = newsites;
        end
        for iPhase = 1:length(PARAMS.Phases)
            for iBand = 1:length(bands)
                evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts = [];
                % collect the events using the current site
                if isfield(evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}), 'evts')
                evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.(sites{1}) = evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts; 
                for iEvt = length(evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).tstart):-1:1
                    evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.(sites{2}){iEvt} = restrict(data.(phases{iPhase}).([sites{1} types{iType}]), evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).tstart(iEvt),evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).tend(iEvt));
                end
                else
                   evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.(sites{1}) = {}; 
                   evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.(sites{2}) = {};
                end
            end
        end
    end
end

