function MS_Phase_Analyses(cfg_in, data, Events)
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



for ii = 1:length(site_list)
    for jj = 1:length(site_list)
        mat_out.labels{ii,jj} = [site_list{ii}(1:end-4) '_' site_list{jj}(1:end-4)];
    end
end

% pre-make the output matrices
for iPhase = 1:length(PARAMS.Phases)
    for iBand = 1:length(bands)
        mat_out.(PARAMS.Phases{iPhase}).COH.(bands{iBand}) = cell(size(mat_out.labels));   % for the coherence
        mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand}) = cell(size(mat_out.labels));  % for the amplitude x-cor
        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand}) = cell(size(mat_out.labels));  % for the amplitude lag
        mat_out.(PARAMS.Phases{iPhase}).PS.(bands{iBand}) = cell(size(mat_out.labels));   % for the phase slope
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
                AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(d_filter.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                AMP.(pairs{iPair}).([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt} = restrict(d_filter.([S{2} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart(iEvt), Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tend(iEvt));
                
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
%%  get the Amplitude x-corr and use to determine if this event should be used.
for iPhase = 1:length(PARAMS.Phases)
    for iPair = 1:length(pairs)
        for iBand = 1:length(bands)
            S = strsplit(pairs{iPair}, '_');
            idx = strfind(mat_out.labels, [S{1} '_' S{2}]);
            [x_idx,y_idx] = find(not(cellfun('isempty', idx)));
            for iEvt = length(Re_Events.([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}).tstart):-1:1
                
                %amp xcorr
                [ac, lag] = xcov(AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data,AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data,100,  'coeff');
                lag = lag * 1/AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.cfg.hdr{1}.SamplingFrequency;
                
                % determine if the event 
                if cfg.cfg_amp.nShuffle > 0
                    clear shuf_max_xcov;
                    
                    for iShuf = cfg.cfg_amp.nShuffle:-1:1
                        %temp_whitenoise = rand(size(vStr_lg))-0.5;
                        temp_whitenoise = AUX_shuffle_phases(AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data);
                        temp_corrvalues = xcov(AMP.(pairs{iPair}).([S{1} cfg.type]).(PARAMS.Phases{iPhase}).(bands{iBand}){iEvt}.data,temp_whitenoise,'coeff');
                        shuf_max_xcov(iShuf) = max(temp_corrvalues);
                    end
                    if max(ac) >= std(shuf_max_xcov)  % if the max ac value is greater than the max than 1sd of the shuffle max values then keep the event
                        mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand}){x_idx, y_idx, iEvt} = ac;
                        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand}){x_idx, y_idx, iEvt} = lag;
                    else
                        mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand}){x_idx, y_idx, iEvt} = NaN;
                        mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand}){x_idx, y_idx, iEvt} = NaN;
                        fprintf(['Event failed shuffle: ' num2str(iEvent) ' in ' bands{iBand} '  ' PARAMS.Phases{iPhase} '\n'])
                    end
                else
                    mat_out.(PARAMS.Phases{iPhase}).AMP_AC.(bands{iBand}){x_idx, y_idx, iEvt} = ac;
                    mat_out.(PARAMS.Phases{iPhase}).AMP_LAG.(bands{iBand}){x_idx, y_idx, iEvt} = lag;
                end
            end
        end
    end
end




