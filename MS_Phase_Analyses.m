function mat_out = MS_Phase_Analyses(cfg_in, data, Events)
%% MS_Phase_Analyses: completes a series of phase related analyses:
%   - coherence: measures the degree of similarity in the change of two
%   signals in time
%
%   - phase offset: used to measure the different in the phase of two
%   signals within a certain frequency range.
%
%   - amplitude x_cor: measure of the timing and correlation between the
%   power envelope of two signals.  Can reveal limited lead/lag
%   relationships
%
%   - Phase Slope (PS): a more robust measure of the lead lag relationship in the phase of two signals.  Can reveal lead/lag relationships.
%
%
%
% Inputs:
%   - cfg_in: [struct]   contains any configuration paramters that deviated from the defaults.
%   - data:   [struct]   this is for a single subject for a single session.  data structure containing the tvec, data, labels, cfg
%   (with hdr) which results from MS_load_data_fast.
%   - events: [struct]   events IV structure for a single subject for a single
%   session. this will also contain the time intervals for all of the
%   events of interest.  In this case for the "low" and "high" gamma
%   events.
%
% Outputs:
%   phase_vals: [struct]
%
%%  Set up defaults
global PARAMS
cfg_def = [];
cfg_def.type = '_pot'; % can also be '_trk'.  For now only works with pot or track and not both.  to get both just run it again with the other cfg.type
cfg_def.check = 1; % used to create a plot for verification.
cfg_def.check_dir = [PARAMS.inter_dir '/phase_check']; %where to save the check figure.


%cfgs for amp x-corr
cfg_def.cfg_amp = [];
cfg_def.cfg_amp.count = 100;
cfg_def.cfg_amp.nShuffle = 100;

% cfgs for power spectral densities (only just for checking plot
cfg_def.cfg_psd = [];
cfg_def.cfg_psd.hann_win = 256;

%resize the events to allow for a more accurate measure given other
%frequencies
cfg_def.cfg_re = [];
cfg_def.cfg_re.d = [-.1 .1];

% Filter the data for the amplitude x_corr
cfg_def.cfg_filter = [];
cfg_def.cfg_filter1.f = [45 65];
cfg_def.cfg_filter2.f = [70 90];


% phase slope congifurations
cfg_def.cfg_phase = [];
cfg_def.cfg_phase.debug = 0;
cfg_def.cfg_phase.freq = [30 100];
cfg_def.cfg_phase.circ_reg.window_size = 9;
cfg_def.cfg_phase.circ_reg.offset_delta = pi/16;

% phase offset using cpsd
cfg_def.cfg_coh.spec_window = 128;
cfg_def.cfg_coh.NFFT = cfg_def.cfg_coh.spec_window*4;

cfg = ProcessConfig2(cfg_def, cfg_in);

%% rename any "piri_x_..." data labels.
for iPhase = 1:length(PARAMS.Phases)
    f_names = fieldnames(data.(PARAMS.Phases{iPhase}));
    for iF = 1:length(f_names)
        if strcmp(f_names{iF}, 'Piri_O_pot')
            data.(PARAMS.Phases{iPhase}).PiriO_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_O_pot');
        elseif strcmp(f_names{iF}, 'Piri_O_trk')
            data.(PARAMS.Phases{iPhase}).PiriO_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_O_trk');
        elseif strcmp(f_names{iF}, 'Piri_N_pot')
            data.(PARAMS.Phases{iPhase}).PiriN_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_N_pot');
        elseif strcmp(f_names{iF}, 'Piri_N_trk')
            data.(PARAMS.Phases{iPhase}).PiriN_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_N_trk');
        end
    end
end

% same for Events
f_names = fieldnames(Events);
for iF = 1:length(f_names)
    if strcmp(f_names{iF}, 'Piri_O_pot')
        Events.PiriO_pot = Events.(f_names{iF});
        Events = rmfield(Events, 'Piri_O_pot');
    elseif strcmp(f_names{iF}, 'Piri_O_trk')
        Events.PiriO_trk = Events.(f_names{iF});
        Events = rmfield(Events, 'Piri_O_trk');
    elseif strcmp(f_names{iF}, 'Piri_N_pot')
        Events.PiriN_pot = Events.(f_names{iF});
        Events = rmfield(Events, 'Piri_N_pot');
    elseif strcmp(f_names{iF}, 'Piri_N_trk')
        Events.PiriN_trk = Events.(f_names{iF});
        Events = rmfield(Events, 'Piri_N_trk');
    end
end
%% get only the specific type.
site = fieldnames(data.post);
expStr = ['*' cfg.type];
regStr = ['^',strrep(strrep(expStr,'?','.'),'*','.{0,}'),'$'];
starts = regexpi(site, regStr);
iMatch = ~cellfun(@isempty, starts);
idx = find(iMatch);
site_list = site(idx);


bands = {'low', 'high'};
%%  create a filtered, envelope, and resized version of the events for each channel
% set the parameters
cfg_re = cfg.cfg_re;
%
for iPhase = 1:length(PARAMS.Phases)
    for iSite = 1:length(site_list)
        for iBand = 1:length(bands)
            % resize the events for that channel, phase, band
            Re_Events.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand}) = ResizeIV(cfg_re, Events.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand}));
            % filter the signal
            cfg_filter = cfg.(['cfg_filter' num2str(iBand)]);
            d_filter.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand}) = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).(site_list{iSite}));
            % create a struct for the amplitude to go into
            d_amp.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand}) = d_filter.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand});
            % put the amplitude data in to the structure.
            d_amp.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand}).data = abs(hilbert(d_filter.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand}).data));
        end
    end
end
%% create a list of pairs
pairs = {};
loop_num = 0;
numS =1:length(site_list);
for iSite = 1:length(site_list)
    others = numS(numS~=iSite);
    for iOther = others
        loop_num = loop_num+1;
        pairs{loop_num} = [site_list{iSite}(1:end-4) '_' site_list{iOther}(1:end-4)];
    end
end



for ii = 1:length(PARAMS.all_sites)
    for jj = 1:length(PARAMS.all_sites)
        mat_out.labels{ii,jj} = [PARAMS.all_sites{ii} '_' PARAMS.all_sites{jj}];
    end
end

% pre-make the output matrices
for iPhase = 1:length(PARAMS.Phases)
    for iBand = 1:length(bands)
        all_evts = [];
        for iSite = 1:length(site_list)
            all_evts = cat(2,all_evts, length(Re_Events.(site_list{iSite}).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart));
        end
        mat_out.(PARAMS.Phases{iPhase}).COH_cxx.(bands{iBand})        = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));   % for the coherence values
        mat_out.(PARAMS.Phases{iPhase}).COH_fxx.(bands{iBand})        = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));   % for the coherence frequency
        mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand})         = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));  % for the amplitude x-cor
        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand})        = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));  % for the amplitude lag
        mat_out.(PARAMS.Phases{iPhase}).AMP_AC_max.(bands{iBand})         = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));  % for the amplitude x-cor
        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG_max.(bands{iBand})        = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));  % for the amplitude lag
        mat_out.(PARAMS.Phases{iPhase}).Phase_lag_cxy.(bands{iBand})  = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));   % for the phase slope
        mat_out.(PARAMS.Phases{iPhase}).Phase_lag_F.(bands{iBand})    = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));   % for the phase slope
        mat_out.(PARAMS.Phases{iPhase}).Phase_lag_mean.(bands{iBand}) = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));   % for the phase slope
        mat_out.(PARAMS.Phases{iPhase}).PS_slope.(bands{iBand})       = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));   % for the phase slope values
        mat_out.(PARAMS.Phases{iPhase}).PS_F.(bands{iBand})           = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));   % for the phase slope frequency
        mat_out.(PARAMS.Phases{iPhase}).EVT_COUNT.(bands{iBand})         = cell(size(mat_out.labels,1),size(mat_out.labels,2), max(all_evts));  % for the amplitude x-cor
    end
end
%% Loop over site pairs for comparisons.  Loop over each event.
for iPhase = 1:length(PARAMS.Phases)
    for iPair = 1:length(pairs)
        for iBand = 1:length(bands)
            S = strsplit(pairs{iPair}, '_');
            for iEvt = length(Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart):-1:1
                % get the raw trace from that channel for each site in
                % the pair
                RAW.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(data.(PARAMS.Phases{iPhase}).([S{1} cfg.type]), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                RAW.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(data.(PARAMS.Phases{iPhase}).([S{2} cfg.type]), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                
                % get the filtered signal for each
                FILT.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(d_filter.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                FILT.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(d_filter.([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                
                % get the power evelop signal for each
                AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(d_amp.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                AMP.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(d_amp.([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                
                AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data = AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data - nanmean(AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data);
                AMP.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data = AMP.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data - nanmean(AMP.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data);
                
                % get the PSD if running a check.
                if cfg.check
                    [PSD.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.pxx, PSD.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.f]   = pwelch(RAW.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data, hanning(cfg.cfg_psd.hann_win), cfg.cfg_psd.hann_win/2, 2*cfg.cfg_psd.hann_win, RAW.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.cfg.hdr{1}.SamplingFrequency);
                    [PSD.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.pxx, PSD.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.f]   = pwelch(RAW.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data, hanning(cfg.cfg_psd.hann_win), cfg.cfg_psd.hann_win/2, 2*cfg.cfg_psd.hann_win, RAW.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.cfg.hdr{1}.SamplingFrequency);
                end
                
            end
        end
    end
end
clear data d_filter d_amp
%%  get the Amplitude x-corr and use to determine if this event should be used.
for iPhase = 1:length(PARAMS.Phases)
    for iPair = 1:length(pairs)
        for iBand = 1:length(bands)
            S = strsplit(pairs{iPair}, '_');
            idx = strfind(mat_out.labels, [S{1} '_' S{2}]);
            [x_idx,y_idx] = find(not(cellfun('isempty', idx)));
            for iEvt = length(Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart):-1:1
                % rename the current event to a shorthand
                this_RAW1 = RAW.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                this_RAW2 = RAW.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                
                this_AMP1 = AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                this_AMP2 = AMP.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                
                this_FILT1 = FILT.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                this_FILT2 = FILT.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                
                this_PSD1 = PSD.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                this_PSD2 = PSD.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt};
                
                %amp xcorr
                [ac, lag] = xcov(this_AMP1.data,this_AMP2.data,100,  'coeff');
                lag = lag * 1/this_AMP1.cfg.hdr{1}.SamplingFrequency;
                
                % determine if the event
                if cfg.cfg_amp.nShuffle > 0
                    clear shuf_max_xcov;
                    
                    for iShuf = cfg.cfg_amp.nShuffle:-1:1
                        %temp_whitenoise = rand(size(vStr_lg))-0.5;
                        temp_whitenoise = AUX_shuffle_phases(this_AMP2.data);
                        temp_corrvalues = xcov(this_AMP1.data,temp_whitenoise,'coeff');
                        shuf_max_xcov(iShuf) = max(temp_corrvalues);
                    end
                    if max(ac) >= std(shuf_max_xcov)  % if the max ac value is greater than the max than 1sd of the shuffle max values then keep the event
                        mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand}){x_idx, y_idx, iEvt} = ac;
                        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand}){x_idx, y_idx, iEvt} = lag;
                        mat_out.(PARAMS.Phases{iPhase}).EVT_COUNT.(bands{iBand}){x_idx, y_idx, iEvt} = 1;
                        [t_max, idx] = max(ac);
                        mat_out.(PARAMS.Phases{iPhase}).AMP_AC_max.(bands{iBand}){x_idx, y_idx, iEvt} = t_max;
                        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG_max.(bands{iBand}){x_idx, y_idx, iEvt} = lag(idx);

                    else
                        mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(ac));
                        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(lag));
                        mat_out.(PARAMS.Phases{iPhase}).AMP_AC_max.(bands{iBand}){x_idx, y_idx, iEvt} = NaN;
                        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG_max.(bands{iBand}){x_idx, y_idx, iEvt} = NaN;
                        mat_out.(PARAMS.Phases{iPhase}).EVT_COUNT.(bands{iBand}){x_idx, y_idx, iEvt} = 0;

                        fprintf(['Event failed shuffle: ' num2str(iEvt) ' in ' bands{iBand} '  ' PARAMS.Phases{iPhase} '\n'])
                    end
                else
                    mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand}){x_idx, y_idx, iEvt} = ac;
                    mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand}){x_idx, y_idx, iEvt} = lag;
                    [t_max, idx] = max(ac);
                    mat_out.(PARAMS.Phases{iPhase}).AMP_AC_max.(bands{iBand}){x_idx, y_idx, iEvt} = t_max;
                    mat_out.(PARAMS.Phases{iPhase}).AMP_LAG_max.(bands{iBand}){x_idx, y_idx, iEvt} = lag(idx);
                end
                
                %% get the coherence the same way.
                if max(ac) >= std(shuf_max_xcov)
                    
                    [cxx, fxx] = mscohere(this_RAW1.data,this_RAW2.data,hanning(cfg.cfg_coh.spec_window),cfg.cfg_coh.spec_window/2, cfg.cfg_coh.NFFT, this_RAW1.cfg.hdr{1}.SamplingFrequency);
                    
                    mat_out.(PARAMS.Phases{iPhase}).COH_cxx.(bands{iBand}){x_idx, y_idx, iEvt} = cxx';
                    mat_out.(PARAMS.Phases{iPhase}).COH_fxx.(bands{iBand}){x_idx, y_idx, iEvt} = fxx';
                    
                else
                    [cxx, fxx] = mscohere(this_RAW1.data,this_RAW2.data,hanning(cfg.cfg_coh.spec_window),cfg.cfg_coh.spec_window/2, cfg.cfg_coh.NFFT, this_RAW1.cfg.hdr{1}.SamplingFrequency);
                    mat_out.(PARAMS.Phases{iPhase}).COH_cxx.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(cxx'));
                    mat_out.(PARAMS.Phases{iPhase}).COH_fxx.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(fxx'));
                end
                %% get the phase lag
                if max(ac) >= std(shuf_max_xcov)
                    
                    [Cxy,F] = cpsd(this_RAW1.data, this_RAW2.data,length(this_RAW1.data),0,cfg.cfg_coh.NFFT,this_RAW1.cfg.hdr{1}.SamplingFrequency);
                    
                    coh_spec_phase= -angle(Cxy); %higher value means leading. outputs radians
                    phase_off =circ_mean(coh_spec_phase(nearest_idx(cfg.(['cfg_filter' num2str(iBand)]).f(1), F):nearest_idx(cfg.(['cfg_filter' num2str(iBand)]).f(2), F)));
                    
                    mat_out.(PARAMS.Phases{iPhase}).Phase_lag_cxy.(bands{iBand}){x_idx, y_idx, iEvt} = coh_spec_phase';
                    mat_out.(PARAMS.Phases{iPhase}).Phase_lag_F.(bands{iBand}){x_idx, y_idx, iEvt} = F';
                    mat_out.(PARAMS.Phases{iPhase}).Phase_lag_mean.(bands{iBand}){x_idx, y_idx, iEvt} =phase_off;
                    
                else
                    [Cxy,F] = cpsd(this_RAW1.data, this_RAW2.data,length(this_RAW1.data),0,cfg.cfg_coh.NFFT,this_RAW1.cfg.hdr{1}.SamplingFrequency);
                    mat_out.(PARAMS.Phases{iPhase}).Phase_lag_cxy.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(Cxy'));
                    mat_out.(PARAMS.Phases{iPhase}).Phase_lag_F.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(F'));
                    mat_out.(PARAMS.Phases{iPhase}).Phase_lag_mean.(bands{iBand}){x_idx, y_idx, iEvt} = NaN;
                end
                %% get the phase slope index.
                    cfg.cfg_phase.Fs = this_RAW1.cfg.hdr{1}.SamplingFrequency;
                    [phase_slopes, F_PS] = MS_phase_slope(cfg.cfg_phase,this_RAW1.data,this_RAW2.data);
               if max(ac) >= std(shuf_max_xcov)
                    mat_out.(PARAMS.Phases{iPhase}).PS_slope.(bands{iBand}){x_idx, y_idx, iEvt} = phase_slopes;
                    mat_out.(PARAMS.Phases{iPhase}).PS_F.(bands{iBand}){x_idx, y_idx, iEvt} = F_PS';
                    
               else               
                    mat_out.(PARAMS.Phases{iPhase}).PS_slope.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(phase_slopes));
                    mat_out.(PARAMS.Phases{iPhase}).PS_F.(bands{iBand}){x_idx, y_idx, iEvt} = NaN(size(F_PS'));
                end
            end
            if cfg.check
                close all
                figure(111)
                
                % plot the raw event
                subplot(2,4,1:2)
                hold on
                plot(this_RAW1.tvec, this_RAW1.data, 'b')
                plot(this_RAW2.tvec, this_RAW2.data, 'r')
                xlim([this_RAW1.tvec(1) this_RAW1.tvec(end)])
                legend({S{1} S{2}}, 'location', 'northwest')
                xlabel('Time (s)')
                
                % plot the PSD
                subplot(2,4,3)
                hold on
                plot(this_PSD1.f, 10*log10(this_PSD1.pxx), 'b');
                plot(this_PSD2.f, 10*log10(this_PSD2.pxx), 'r');
                %                 legend({S{1} S{2}}, 'location', 'northwest')
                xlim([0 100])
                text(5, max(10*log10(this_PSD1.pxx)) +10, [S{1} '-' S{2} ' Event: ' num2str(iEvt) ' ' bands{iBand}])
                ylabel('Power')
                legend({S{1} S{2}}, 'location', 'northwest')
                xlabel('Frequency')
                
                % plot the phase and amp
                subplot(2,4,5:6)
                hold on
                plot(this_FILT1.tvec, this_FILT1.data, 'b')
                plot(this_AMP1.tvec, this_AMP1.data, 'b')
                
                plot(this_FILT2.tvec, this_FILT2.data, 'r')
                plot(this_AMP2.tvec, this_AMP2.data, 'r')
                
                xlim([this_AMP1.tvec(1) this_AMP1.tvec(end)])
                y_l = get(gca, 'ylim');
                [ac_max, idx] =  max(ac);
                s_amp = sprintf('Amp xcorr = %.2f\nLag = %.2fms', ac_max, lag(idx)*1000);
                text(this_AMP1.tvec(5), y_l(2)-(y_l(2)*.1), s_amp, 'fontsize', 12);
                xlabel('Frequency')
                
                % coherence between
                subplot(2,4,7)
                plot(fxx, cxx)
                xlim([0 100])
                vline([cfg.cfg_filter1.f,cfg.cfg_filter2.f], {'--b', '--b', '--g', '--g'})
                ylabel('Coherence')
                xlabel('Frequency')
                % get the phase offset value and plot in the corner
                subplot(2,4,4)
                hold on
                plot(F, rad2deg(coh_spec_phase), 'b');
                vline([cfg.cfg_filter1.f,cfg.cfg_filter2.f], {'--b', '--b', '--g', '--g'})
                xlim([0 100])
                plot([0:100], zeros(101), 'k--')
                text(5, 20, ['Offset = ' num2str(rad2deg(phase_off),3) ' Degrees'], 'fontsize', 12);
                xlabel('Frequency')
                ylabel('Phase offset (deg)')
                
                % plot the phase slope
                subplot(2,4,8)
                plot(F_PS, phase_slopes)
                xlim(cfg.cfg_phase.freq)
                vline([cfg.cfg_filter1.f,cfg.cfg_filter2.f], {'--b', '--b', '--g', '--g'})
                hline(0, '--k')
                ylabel('phase/deg')
                xlabel('Frequency')
                cfg_fig.pos = [145 45 1200 750];
                % SetFigure(cfg_fig, gcf)
                saveas(gcf, [PARAMS.inter_dir 'Phase_checks/' cfg.Subject_id '_' S{1} '_' S{2} '_Event_' num2str(iEvt) '_' bands{iBand} '.png'])
            end
        end
    end
end




