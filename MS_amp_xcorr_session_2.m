function Amp = MS_amp_xcorr_session_2(cfg_in, data)

%% temp
% cfg_in = [];
% global PARAMS
%
% load([PARAMS.inter_dir 'R102_Data.mat'])
% data = data.R102_2016_09_24;
%
% load([PARAMS.inter_dir 'R102_Events.mat'])

%% MS_amp_xcorr: gets the amplitude cross correlation for eahc event across a range of frequencies specified in cfg.freq in steps of cfg.freq_step (default is 1Hz) for a given pair of channels.
%
%
%   - amplitude x_cor: measure of the timing and correlation between the
%   power envelope of two signals.  Can reveal limited lead/lag
%   relationships
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
cfg_def.check_dir = [PARAMS.inter_dir 'phase_check']; %where to save the check figure.

cfg_def.cfg_filter = [];
cfg_def.cfg_filter.freq = [3 100];
cfg_def.cfg_filter.freq_step = 1;

%cfgs for amp x-corr
cfg_def.cfg_amp = [];
cfg_def.cfg_amp.count = 100;
cfg_def.cfg_amp.nShuffle = 100;

% cfgs for amp-corr-ogram
cfg_def.cfg_amp_cor_ogram = [];
cfg_def.cfg_amp_cor_ogram.dT = 0.5;

cfg = ProcessConfig2(cfg_def, cfg_in);
%% rename any "piri_x_..." data labels.
for iPhase = 1:length(PARAMS.Phases)
    f_names = fieldnames(data.(PARAMS.Phases{iPhase}));
    for iF = 1:length(f_names)
         if strcmp(f_names{iF}, 'Piri_O_pot')
            data.(PARAMS.Phases{iPhase}).PiriO_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_O_pot');
        elseif   strcmp(f_names{iF}, 'Piri_OFC_pot')
            data.(PARAMS.Phases{iPhase}).PiriO_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_OFC_pot');
        elseif strcmp(f_names{iF}, 'Piri_O_trk')
            data.(PARAMS.Phases{iPhase}).PiriO_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_O_trk');
        elseif   strcmp(f_names{iF}, 'Piri_OFC_trk')
            data.(PARAMS.Phases{iPhase}).PiriO_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_OFC_trk');
        elseif strcmp(f_names{iF}, 'Piri_N_pot')
            data.(PARAMS.Phases{iPhase}).PiriN_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_N_pot');
        elseif strcmp(f_names{iF}, 'Piri_NAc_pot')
            data.(PARAMS.Phases{iPhase}).PiriN_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_NAc_pot');
        elseif strcmp(f_names{iF}, 'Piri_N_trk')
            data.(PARAMS.Phases{iPhase}).PiriN_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_N_trk');
        elseif strcmp(f_names{iF}, 'Piri_NAc_trk')
            data.(PARAMS.Phases{iPhase}).PiriN_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_NAc_trk');
        end
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


% temporary
% pairs = {'PL_NAc', 'PL_CG', 'PL_PiriO','NAc_CG', 'NAc_PiriO', 'CG_PiriO'};


%% loop over pairs for each frequency band
for iPhase = 1:length(PARAMS.Phases)
    completed_pairs = {};
    for  iPair = 1:length(pairs)
        Same_flag = 0;
        %     iPair = 5;
        S = strsplit(pairs{iPair}, '_');
        % make sure oyu only have one copy of each pair.
        for iC = 1:length(completed_pairs)
            if strcmp(completed_pairs{iC}, [S{2} '_' S{1}])
                Same_flag = 1;
            end
        end
        
       if Same_flag ~=1
            % set up an NaN array
            Freq_list = cfg.cfg_filter.freq(1): cfg.cfg_filter.freq_step:cfg.cfg_filter.freq(end);
            
            times = data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).tvec(1):cfg.cfg_amp_cor_ogram.dT:data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).tvec(end);
            
            amp_ac.(pairs{iPair}).(PARAMS.Phases{iPhase}) =NaN(length(Freq_list), length(times));
            
            % S = strsplit(pairs{iPair}, '_');
            fprintf(['Processing: Frequency (Hz)  '])
            Amp.f.(pairs{iPair}).(PARAMS.Phases{iPhase}) = NaN(1,length(Freq_list));
            Amp.ac.(pairs{iPair}).(PARAMS.Phases{iPhase}) = NaN(1,length(Freq_list));
            Amp.lag.(pairs{iPair}).(PARAMS.Phases{iPhase}) = NaN(1,length(Freq_list));
            
            for iF = 1:length(Freq_list)
                
                this_F = Freq_list(iF);
                
                
                if iF>1
                    for j=0:log10(this_F-1)
                        fprintf('\b'); % delete previous counter display
                    end
                end
                fprintf('%d', this_F);
                %%
                cfg_filter= [];
                cfg_filter.display_filter = 0;
                cfg_filter.f = [this_F-1 this_F+1];
%                 cfg_filter.wp = cfg_filter.f*(2/data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).cfg.hdr{1}.SamplingFrequency);
%                 cfg_filter.ws = [cfg_filter.f(1)-2.5 cfg_filter.f(2)+2.5]*(2/data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).cfg.hdr{1}.SamplingFrequency);
%                 cfg_filter.R = 0.5;
                cfg_filter.order = 4;%cheb2ord(cfg_filter.wp, cfg_filter.ws, cfg_filter.R,20);
                cfg_filter.verbose = 0;
                cfg_filter.type = 'fdesign';%'cheby1';

                FILT_1 = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).([S{1} '_pot']));
                FILT_2 = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).([S{2} '_pot']));
                
                AMP_1 = FILT_1;
                AMP_2 = FILT_2;
                
                AMP_1.data = abs(hilbert(FILT_1.data));
                AMP_2.data = abs(hilbert(FILT_2.data));
                
                [ac, lag] = xcov(AMP_1.data,AMP_2.data,200,  'coeff');
                lag = lag * 1/AMP_1.cfg.hdr{1}.SamplingFrequency;
                [ac_max,idx] = max(ac);
                lag_max = lag(idx);
                
                Amp.f.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF) = this_F;
                Amp.ac.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF) = ac_max;
                Amp.lag.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF) = lag_max;
                
                
                %% create a Amp-xcorr-ogram
                
                %         for iT = 1:length(times)-1
                %
                %             D_1_temp = restrict(AMP_1, times(iT), times(iT+1));
                %             D_2_temp = restrict(AMP_2, times(iT), times(iT+1));
                %             %                 Dxx1 = restrict(FILT_1, times(iT), times(iT+1)); checks for debugging
                %             %                 Dxx2 = restrict(FILT_2, times(iT), times(iT+1));
                %
                %
                %             [ac, lag] = xcov(D_1_temp.data,D_2_temp.data,200,  'coeff');
                %             max_ac = ac(lag==0);
                %             amp_ac.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF, iT) = max_ac;
                %             amp_ac_t.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF, iT) = times(iT);
                %
                %         end
                
                
            end;
            fprintf(['\n' S{1} '-' S{2} '...complete\n'])
            completed_pairs{end+1} = pairs{iPair};
        else
            fprintf(['\n' 'Skipping ' S{1} '-' S{2} '  because that pair has already been processed\n'])
        end
    end
end

% amp_out.amp_ac = amp_ac;
% amp_out.amp_ac_t = amp_ac_t;
% amp_out.AMP = Amp;
% amp_out.cfg = cfg;

% save([PARAMS.inter_dir 'R102_amp.mat'], 'amp_out','-v7.3')
