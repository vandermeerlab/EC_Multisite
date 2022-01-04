function [Naris_out, cfg] = MS_get_bandpower(cfg_in, data_in, Naris_in)
%% MS_get_bandpower : get the area under the psd curve for each experimental condition. This is an update of MS_get_power_ratio
%

%          Inputs:
%           - cfg_in [struct] config file with parameters
%           - data_in [struct] contains the csc data 
%
%          Outputs:
%           - Naris_out [struct] same as the input but with the power ratio
%           appended
%
%
% EC - 2021-05-23 
global PARAMS
cfg_def = [];
cfg_def.gamma_freq = [45, 65; 70, 90];
cfg_def.contrast = [25 45; 90, 110];
cfg = ProcessConfig2(cfg_def, cfg_in);

%% append to Naris_in if given

if nargin < 3
    Naris_in = [];
end

Naris_out = Naris_in; % copy the Naris in struct to append the bandpower

%% compute the bandpower for each session/site/condition
phases = PARAMS.Phases;

sites = fieldnames(data_in.pre);
for iSite = 1:length(sites)
    if strcmpi(sites{iSite}, 'pos') || strcmpi(sites(iSite), 'ExpKeys')
        continue
    else

        for iPhase = 1:length(phases)
            fprintf('\n<strong>%s</strong> Processing Phase: %s Site: %s',mfilename, phases{iPhase}, sites{iSite})
            % get the bandpower for low/high gamma
            Naris_out.(phases{iPhase}).(sites{iSite}).bandpow.low = bandpower(data_in.(phases{iPhase}).(sites{iSite}).data, data_in.(phases{iPhase}).(sites{iSite}).cfg.hdr{1}.SamplingFrequency, cfg.gamma_freq(1,:));
            Naris_out.(phases{iPhase}).(sites{iSite}).bandpow.high = bandpower(data_in.(phases{iPhase}).(sites{iSite}).data, data_in.(phases{iPhase}).(sites{iSite}).cfg.hdr{1}.SamplingFrequency, cfg.gamma_freq(2,:));
            % get the bandpower for low/high contrast bands
            Naris_out.(phases{iPhase}).(sites{iSite}).bandpow.cont_low = bandpower(data_in.(phases{iPhase}).(sites{iSite}).data, data_in.(phases{iPhase}).(sites{iSite}).cfg.hdr{1}.SamplingFrequency, cfg.contrast(1,:));
            Naris_out.(phases{iPhase}).(sites{iSite}).bandpow.cont_high = bandpower(data_in.(phases{iPhase}).(sites{iSite}).data, data_in.(phases{iPhase}).(sites{iSite}).cfg.hdr{1}.SamplingFrequency, cfg.contrast(2,:));
        end
        % make the 'control' coniditon using the mean of the pre and post. 
        Naris_out.control.(sites{iSite}).bandpow.low = mean([Naris_out.pre.(sites{iSite}).bandpow.low; Naris_out.post.(sites{iSite}).bandpow.low]); 
        Naris_out.control.(sites{iSite}).bandpow.high = mean([Naris_out.pre.(sites{iSite}).bandpow.high; Naris_out.post.(sites{iSite}).bandpow.high]); 
               
        Naris_out.control.(sites{iSite}).bandpow.cont_low = mean([Naris_out.pre.(sites{iSite}).bandpow.cont_low; Naris_out.post.(sites{iSite}).bandpow.cont_low]); 
        Naris_out.control.(sites{iSite}).bandpow.cont_high = mean([Naris_out.pre.(sites{iSite}).bandpow.cont_high; Naris_out.post.(sites{iSite}).bandpow.cont_high]); 
 
    end
end
fprintf('\n')
