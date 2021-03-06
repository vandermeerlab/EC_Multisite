function [mat_out] = MS_get_phase_metrics_serial(cfg_in, Events, Data)
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
cfg_def.resize.d = [-0.1 0.1]; % add data to each end of the event. 
cfg_def.f = [45 65; 70 90];
cfg_def.nShuffle = 10000;
cfg_def.debug = 0; % used to check each event with a plot showing both the filtered data + envelope and the Amp_xcorr
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
        mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_coh = cell(size(mat_out.labels));
%         mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).sess_amp = cell(size(mat_out.labels));
        mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_amp_ac = cell(size(mat_out.labels));
        mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_amp_lag = cell(size(mat_out.labels));
        mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).PS = cell(size(mat_out.labels));
    end
end

%% cycle through events
t_type = 'pot';

for iPair = 1:length(pairs)
    [S_out, ~] = strsplit(pairs{iPair}, '_');
    if strfind(S_out{1}, 'PiriO')  % walk around if the pair name for piriform sites don't match the data names beacuse of an underscore
         S_out{1} = strrep(S_out{1}, 'PiriO', 'Piri_O');
    elseif strfind(S_out{2}, 'PiriO')
        S_out{2} = strrep(S_out{2}, 'PiriO', 'Piri_O');
    end
    if strfind(S_out{1}, 'PiriN')  % walk around if the pair name for piriform sites don't match the data names beacuse of an underscore
         S_out{1} = strrep(S_out{1}, 'PiriN', 'Piri_N');
    elseif strfind(S_out{2}, 'PiriN')
         S_out{2} = strrep(S_out{2}, 'PiriN', 'Piri_N');
    end
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
        for iPhase = 2:length(PARAMS.Phases)
            for iBand = 1:length(bands)
                evts = Events.([S1 '_' t_type]).(PARAMS.Phases{iPhase}).(bands{iBand});
                this_data_1 = Data.(PARAMS.Phases{iPhase}).([S1 '_' t_type]);
                this_data_2 = Data.(PARAMS.Phases{iPhase}).([S2 '_' t_type]);
                
                
                % Filter the data for the amplitude x_corr
                cfg_filter = [];
                cfg_filter.f = [cfg.f(iBand,1) cfg.f(iBand,2)];
                
                % filter the two signals
                d_filter1 = FilterLFP(cfg_filter, this_data_1);
                d_filter2 = FilterLFP(cfg_filter, this_data_2);
                
                % get the evelope
                d_abs_1 = d_filter1;
                d_abs_2 = d_filter2;
                
                d_abs_1.data = abs(hilbert(d_abs_1.data));
                d_abs_2.data = abs(hilbert(d_abs_2.data));
                
                
                % %% phase slope index
                % phase_slopes = {};
                % cycle through all events
                all_ac = []; all_coh = []; all_lag = []; phase_slopes = []; 
                for iEvent = length(evts.tstart):-1:1
                    this_event = SelectIV([],evts,iEvent);
                    this_event = ResizeIV(cfg.resize,this_event);
                    
                    d_f_1 = restrict(d_filter1,this_event);
                    d_f_2 = restrict(d_filter2,this_event);
                    
                    d1 = restrict(this_data_1,this_event);
                    d2 = restrict(this_data_2,this_event);
                    
                    %% get the phase angle within each event. 
                    cfg_coh.spec_window = 256;
                    cfg_coh.NFFT = 1024;
                    
                     [Cxy,F] = cpsd(d1.data,d2.data,hamming(cfg_coh.spec_window),cfg_coh.spec_window/2,cfg_coh.NFFT,d1.cfg.hdr{1}.SamplingFrequency);
                     coh_spec_phase= -angle(Cxy); %higher value means leading. outputs radians
                     all_coh(iEvent) =rad2deg(circ_mean(coh_spec_phase(nearest_idx(cfg_filter.f(1), F):nearest_idx(cfg_filter.f(2), F))));
                        
                    
                    %% get the amp x-corr
                    amp_ev1 = restrict(d_abs_1, this_event);
                    amp_ev1 = amp_ev1.data - nanmean(amp_ev1.data);
                    amp_ev2 = restrict(d_abs_2, this_event);
                    amp_ev2 = amp_ev2.data - nanmean(amp_ev2.data);
                    
                    % get the amplitue cross-correlation "amp_ac" and lag
                    % time
                    [ac, lag] = xcov(amp_ev1,amp_ev2,100,  'coeff');
                    lag = lag * 1/d1.cfg.hdr{1}.SamplingFrequency;
                    %
                    if cfg.debug
                        figure(100)
                        subplot(2,1,1)
                        plot(1:length(amp_ev1), amp_ev1, 'r', 1:length(amp_ev2), amp_ev2, 'b')
                        hold on
                        plot(1:length(d_f_1.tvec), d_f_1.data, 'r', 1:length(d_f_2.tvec), d_f_2.data, 'b')
                        subplot(2,1,2)
                        plot(lag, ac,'b')
                        pause(.5)
                        close all
                    end
                    if cfg.nShuffle > 0
                        clear shuf_max_xcov;
                        
                        for iShuf = cfg.nShuffle:-1:1
                            %temp_whitenoise = rand(size(vStr_lg))-0.5;
                            temp_whitenoise = AUX_shuffle_phases(amp_ev2);
                            temp_corrvalues = xcov(amp_ev1,temp_whitenoise,'coeff');
                            shuf_max_xcov(iShuf) = max(temp_corrvalues);
                        end
                        if max(ac) >= std(shuf_max_xcov)  % if the max ac value is greater than the max than 1sd of the shuffle max values then keep the event
                            all_ac(iEvent,:) = ac; all_lag(iEvent,:) = lag;
                        else
                            all_ac(iEvent,:) = NaN; all_lag(iEvent,:) = NaN;
                            fprintf(['Event failed shuffle: ' num2str(iEvent) ' in ' bands{iBand} '  ' PARAMS.Phases{iPhase} '\n'])
                        end
                    else
                        all_ac(iEvt,:) = ac; all_lag(iEvt,:) = lag;
                    end
                    
                    %try the function
                    if max(ac) >= std(shuf_max_xcov)
                        cfg_phase = [];
                        cfg_phase.Fs = d1.cfg.hdr{1}.SamplingFrequency;
                        cfg_phase.debug = 0;
                        cfg_phase.circ_reg.window_size = 9; 
                        cfg_phase.circ_reg.offset_delta = pi/16; 
                        phase_slopes(iEvent) = MS_phase_slope(cfg_phase, d1.data, d2.data);
                    else
                        phase_slopes(iEvent) = NaN;
                    end
                    close all
                end
                fprintf(['Phase Measures: ' num2str(length(evts.tstart))  ' Events in, ' num2str(length(evts.tstart) - sum(isnan(phase_slopes))) ' Events out '  bands{iBand} '  ' PARAMS.Phases{iPhase} '\n'])
                mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_coh{x_idx, y_idx} = all_coh;
                mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_amp_lag{x_idx, y_idx} = all_lag;
                mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).evt_amp_ac{x_idx, y_idx} = all_ac;
                mat_out.(PARAMS.Phases{iPhase}).(bands{iBand}).PS{x_idx, y_idx} = phase_slopes;
            end
        end
    end
end

end