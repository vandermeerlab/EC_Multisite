function Naris_out = MS_get_power_ratio(cfg_in, Naris_in)
%%        :
%
%
%
%          Inputs:
%           - cfg_in [struct] config file with parameters
%           - Naris [struct] contains the output from MS_collect_psd
%           -
%          Outputs:
%           - Naris [struct] same as the input but with the power ratio
%           appended
%           -
%           -
%
% EC - 2017-01-08

cfg_def = [];
cfg_def.gamma_freq = [45, 65; 70, 90];
cfg_def.contrast = [25 45; 90, 110];
cfg_def.method = 'raw';   % can be "ratio" which uses the best ration between gamma and between 20-30 hz
cfg_def.plot =1; 
cfg = ProcessConfig2(cfg_def, cfg_in);

%% get some psds for the pole of interest
phases = fieldnames(Naris_in);
bands = {'low', 'high'};
for iband = 1:2
    for iPhase = 1:length(phases)
        sites = fieldnames(Naris_in.(phases{iPhase}));
        for iSite = 1:length(sites)
            w_psd = Naris_in.(phases{iPhase}).(sites{iSite}).psd.White_Pxx;
            w_F = Naris_in.(phases{iPhase}).(sites{iSite}).psd.White_F;
            psd = Naris_in.(phases{iPhase}).(sites{iSite}).psd.Pxx;
            F = Naris_in.(phases{iPhase}).(sites{iSite}).psd.F;
            
            f_idx = find(F > cfg.gamma_freq(iband,1) & F <= cfg.gamma_freq(iband,2));
            f_idx_contrast = find(F > cfg.contrast(iband,1) & F <= cfg.contrast (iband,2));
            
            w_f_idx = find(w_F > cfg.gamma_freq(iband,1) & w_F <= cfg.gamma_freq(iband,2));
            w_f_idx_contrast = find(w_F > cfg.contrast(iband,1) & w_F <= cfg.contrast (iband,2));
            
            %% calculate the raw and the ratio scores
            Naris_out.(phases{iPhase}).(sites{iSite}).psd.raw.(bands{iband}) = nanmean(10*log10(psd(f_idx)));
            Naris_out.(phases{iPhase}).(sites{iSite}).psd.ratio.(bands{iband}) = mean(trapz(10*log10(psd(f_idx)))) / mean(trapz(10*log10(psd(f_idx_contrast)))) ;
            Naris_out.(phases{iPhase}).(sites{iSite}).psd.white_raw.(bands{iband}) = nanmean(10*log10(w_psd(w_f_idx)));
            Naris_out.(phases{iPhase}).(sites{iSite}).psd.white_ratio.(bands{iband})= mean(trapz(10*log10(w_psd(w_f_idx))))/ mean(trapz(10*log10(w_psd(w_f_idx_contrast))));
            
            fields = fieldnames(Naris_out.(phases{iPhase}).(sites{iSite}).psd);
            for ifield = 1:length(fields)
                plot_mat.(fields{ifield}).(bands{iband})(iSite, iPhase) = Naris_out.(phases{iPhase}).(sites{iSite}).psd.(fields{ifield}).(bands{iband});
            end
        end
    end
end
Naris_in.(phases{iPhase}).(sites{iSite}).psd.plot_mat = plot_mat; % keep the phase by site matrix of each output for simple plotting. 

%% quick check
if cfg.plot == 1
    for iband = 1:2
        figure(iband)
        for ifield = 1:length(fields)
            subplot(2,2,ifield)
            imagesc(plot_mat.(fields{ifield}).(bands{iband}))
            set(gca, 'ytick', [1:4], 'ytickLabel', sites, 'xticklabel', phases)
            xlabel(fields{ifield})
        end
    end
end



end