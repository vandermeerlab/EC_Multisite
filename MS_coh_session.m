function Coh = MS_coh_session(cfg_in, data)

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

cfg_def.cfg_coh.spec_window = 512;
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

            fprintf(['Processing: Frequency (Hz)  '])

            % compute the coherence. 
            [cxx, fxx] = mscohere(data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).data , data.(PARAMS.Phases{iPhase}).([S{2} '_pot']).data,...
                hanning(cfg.cfg_coh.spec_window),cfg.cfg_coh.spec_window/2, cfg.cfg_coh.NFFT,...
                data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).cfg.hdr{1}.SamplingFrequency);

                 Coh.fxx.(pairs{iPair}).(PARAMS.Phases{iPhase}) = fxx; % keep the coherence values
                 Coh.cxx.(pairs{iPair}).(PARAMS.Phases{iPhase}) = cxx; % keep the 
                
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
