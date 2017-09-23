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
%% create a line fit for the power ratio compariso using the average psd
% from the pre and post conditions and using an exponential fit on the data
% from 1 -120Hz.
ii = 200;
type = {'White_Pxx','White_F'};
sites = fieldnames(Naris_in.pre);
for iSite = 1:length(sites)
    Naris_in.control.(sites{iSite}).psd.(type{1}) = mean([Naris_in.pre.(sites{iSite}).psd.(type{1}), Naris_in.post.(sites{iSite}).psd.(type{1})],2);
    Naris_in.control.(sites{iSite}).psd.White_F = Naris_in.pre.(sites{iSite}).psd.(type{2});
    idx_low = nearest_idx3(1, Naris_in.control.(sites{iSite}).psd.(type{2}));
    idx_high = nearest_idx3(ii, Naris_in.control.(sites{iSite}).psd.(type{2}));
    C = fit(Naris_in.control.(sites{iSite}).psd.(type{2})(idx_low:idx_high),10*log10(Naris_in.control.(sites{iSite}).psd.(type{1})(idx_low:idx_high)), 'exp2');
    Control.(sites{iSite}) = C(Naris_in.control.(sites{iSite}).psd.(type{2}));
end

%% temp plot
figure(222); subplot(1,2,1)
phases = fieldnames(Naris_in);
c_ord = linspecer(7); loop = 0;
for iSite = 3%:length(sites);
%         subplot(2,4,iSite)
    title(sites{iSite})
    for iPhase = 1:length(phases)
        hold on
        loop = loop+1;   
        plot(Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{2}), 10*log10(Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{1})),'color', c_ord(loop,:) );
    end
    
    % for ii = [80:20:200];
    plot(Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{2}), Control.(sites{iSite}),'-k')
    %     phases{leg_id} = num2str(ii);
    % end
    legend(phases)
    xlim([0 100])
end
%% try to get the AOC for the phases relative to the exp_curve fit
% close all
bands = {'low', 'high'};
phases = fieldnames(Naris_in);
AOC_out = [];
for iPhase = 1:length(phases)
    for iSite = 1:length(sites)
        for iband = 1:2
            psd = Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{1});
            F = Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{2});
            f_idx = find(F > cfg.gamma_freq(iband,1) & F <= cfg.gamma_freq(iband,2));
            
            f_idx_contrast = find(F > cfg.contrast(iband,1) & F <= cfg.contrast (iband,2));
            
            AOC_val_con(iPhase, iSite) =  trapz(F(f_idx_contrast) ,10*log10(psd(f_idx_contrast)) -Control.(sites{iSite})(f_idx_contrast)); 
            AOC_val.(bands{iband})(iPhase, iSite) =  trapz(F(f_idx) ,10*log10(psd(f_idx)) -Control.(sites{iSite})(f_idx)) ;
        end
    end
end
figure(223)
subplot(3,1,1)
bar(AOC_val.low)
set(gca, 'xticklabel', phases)
title('Low Gamma')
hline(1)
legend(sites)

subplot(3,1,2)
bar(AOC_val.high)
set(gca, 'xticklabel', phases)
title('High Gamma')
hline(1)
legend(sites)

subplot(3,1,3)
bar(AOC_val_con)
set(gca, 'xticklabel', phases)
title('Contrast')
% ylim([0 2]);


%% test space


%% get some psds values for the band of interest
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
            set(gca, 'ytick', [1:length(sites)], 'ytickLabel', sites, 'xticklabel', phases)
            xlabel(fields{ifield})
        end
    end
end



end