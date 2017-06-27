function [phase_coh_out] = MS_coh_matrix(trials, ExpKeys)
%% MS_coh_matrix: 
%
%
%
%
%         OUTPUT MATRIX
%
%         |OFC | NAC | PL | CG | Piri_OFC | Piri_NAc
%     OFC | % 
%     NAC |       %           [High gamma]
%     PL  |            %
%     CG  |                 %
%Piri_OFC |                         %
%Piri_NAc |      [low gamma]                  %
%
%
%
%
%
pairs = ExpKeys.GoodPairs;
bands = {'low', 'high'};
phases = {'pre', 'ipsi', 'contra', 'post'};
types = {'_pot', '_trk'};
    
if isempty(pairs) || isempty(pairs{1})
    for iPhase = 1:length(phases)
        phase_coh_out.(phases{iPhase}).low = NaN*ones(5,5);
        phase_coh_out.(phases{iPhase}).high = NaN*ones(5,5);
    end
    all_trials = [];
else
    % parameters
    f_bandpass = {[45 65]; [70 90]};
   cfg = [];
   cfg.Fs =2000;
   cfg.wsize = 128;
    % compute the phase difference between each channel and the reference using cross-spectral power desnity.
    % a negative values represent lags, while positives are leads.v
                        Add_ft()

    for iPairs = 1:length(pairs)
        sites = strsplit(pairs{iPairs}, '_');
        for iBand = 1:length(bands)
            for iType = 1:length(types)
            for iPhase = 1:length(phases)
                if ~isfield(trials.([sites{iSite} types{iType}]).(phases{iPhase}), (bands{iBand}))
                    continue
                else
                    % filter the entire signal
%                     keep = any(~cellfun('isempty',trials.([sites{iSite} types{iType}]).(phases{iPhase}).(bands{iBand}).evts), 1);
%                     trials.(phases{iPhase}).(pairs{iPairs}).(bands{iBand}) = trials.(phases{iPhase}).(pairs{iPairs}).(bands{iBand})(:,keep);
                    % loop over events and compare the phase of
                    for itrial = 1:length(trials.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts)
                        
                        [COH,F] = mscohere(trials.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts{itrial}.data,trials.([sites{2} types{iType}]).(phases{iPhase}).(bands{iBand}).evts{itrial}.data,hanning(cfg.wsize/2),cfg.wsize/4,cfg.wsize,cfg.Fs);
                        all_coh_spec.(phases{iPhase}).(bands{iBand})(itrial) =nanmean(COH(nearest(F,f_bandpass{iBand}(1)):nearest(F,f_bandpass{iBand}(2))));
                        
                        [COH,F] = mscohere(trials.(phases{iPhase}).(pairs{iPairs}).(bands{iBand}){itrial}(1,:),trials.(phases{iPhase}).(pairs{iPairs}).(bands{iBand}){itrial}(2,:),hanning(cfg.wsize/2),cfg.wsize/4,cfg.wsize,cfg.Fs);
                        all_coh_spec_21.(phases{iPhase}).(bands{iBand})(itrial) = nanmean(COH(nearest(F,f_bandpass{iBand}(1)):nearest(F,f_bandpass{iBand}(2))));
                        %                     all_trials.(Phases{iPhase}).(pairs{iPairs}).(type{itype}){itrial} = [data_in_1.data(cycle_idx(cycle.prev.idx(1):cycle.next.idx(2))); data_in_2.data(cycle_idx(cycle.prev.idx(1):cycle.next.idx(2)))];
                    end
                    all_coh_spec.(phases{iPhase}).(bands{iBand})(all_coh_spec.(phases{iPhase}).(bands{iBand})==0) =[];
                    all_coh_spec_21.(phases{iPhase}).(bands{iBand})(all_coh_spec_21.(phases{iPhase}).(bands{iBand})==0) =[];
                    if strcmp(pairs{iPairs}, 'OFC_NAc')
                        if isempty(all_coh_spec.(phases{iPhase}).(bands{iBand}))
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(1,2) = NaN;
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(2,1) = NaN;
                        else
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(1,2) = rad2deg(circ_mean(all_coh_spec.(phases{iPhase}).(bands{iBand})'));
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(2,1) = rad2deg(circ_mean(all_coh_spec_21.(phases{iPhase}).(bands{iBand})'));
                        end
                    elseif strcmp(pairs{iPairs}, 'OFC_CG')
                        if isempty(all_coh_spec.(phases{iPhase}).(bands{iBand}))
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(1,3) = NaN;
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(3,1) = NaN;
                        else
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(1,3) = rad2deg(circ_mean(all_coh_spec.(phases{iPhase}).(bands{iBand})'));
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(3,1) = rad2deg(circ_mean(all_coh_spec_21.(phases{iPhase}).(bands{iBand})'));
                        end
                    elseif strcmp(pairs{iPairs}, 'PL_OFC')
                        if isempty(all_coh_spec.(phases{iPhase}).(bands{iBand}))
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(1,4) = NaN;
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(4,1) = NaN;
                        else
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(1,4) = rad2deg(circ_mean(all_coh_spec.(phases{iPhase}).(bands{iBand})'));
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(4,1) = rad2deg(circ_mean(all_coh_spec_21.(phases{iPhase}).(bands{iBand})'));
                        end
                    elseif strcmp(pairs{iPairs}, 'NAc_CG')
                        if isempty(all_coh_spec.(phases{iPhase}).(bands{iBand}))
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(2,3) = NaN;
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(3,2) = NaN;
                        else
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(2,3) = rad2deg(circ_mean(all_coh_spec.(phases{iPhase}).(bands{iBand})'));
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(3,2) = rad2deg(circ_mean(all_coh_spec_21.(phases{iPhase}).(bands{iBand})'));
                        end
                    elseif strcmp(pairs{iPairs}, 'PL_NAc')
                        if isempty(all_coh_spec.(phases{iPhase}).(bands{iBand}))
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(2,4) = NaN;
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(4,2) = NaN;
                        else
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(2,4) = rad2deg(circ_mean(all_coh_spec.(phases{iPhase}).(bands{iBand})'));
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(4,2) = rad2deg(circ_mean(all_coh_spec_21.(phases{iPhase}).(bands{iBand})'));
                        end
                    elseif strcmp(pairs{iPairs}, 'PL_CG')
                        if isempty(all_coh_spec.(phases{iPhase}).(bands{iBand}))
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(3,4) = NaN;
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(4,3) = NaN;
                        else
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(3,4) = rad2deg(circ_mean(all_coh_spec.(phases{iPhase}).(bands{iBand})'));
                            phase_coh_out.(phases{iPhase}).(bands{iBand})(4,3) = rad2deg(circ_mean(all_coh_spec_21.(phases{iPhase}).(bands{iBand})'));
                        end
                    end
                end
            end
        end
    end
                        Remove_ft()

end