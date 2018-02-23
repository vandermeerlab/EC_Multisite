function [Naris_out, cfg] = MS_get_power_ratio(cfg_in, Naris_in)
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
global PARAMS
cfg_def = [];
cfg_def.gamma_freq = [45, 65; 70, 90];
cfg_def.contrast = [25 45; 90, 110];
cfg_def.method = 'raw';   % can be "ratio" which uses the best ration between gamma and between 20-30 hz
cfg_def.plot =1;
cfg_def.id = 'no_id';
cfg = ProcessConfig2(cfg_def, cfg_in);

%% set up I/O
Naris_out = Naris_in;

%% create a line fit for the power ratio compariso using the average psd
% from the pre and post conditions and using an exponential fit on the data
% from 1 -120Hz.
ii = 200;
type = {'Pxx', 'F'; 'White_Pxx','White_F'};
for iType = 1:length(type);
sites = fieldnames(Naris_in.pre);
for iSite = 1:length(sites)
    Naris_in.control.(sites{iSite}).psd.(type{iType,1}) = mean([Naris_in.pre.(sites{iSite}).psd.(type{iType,1}), Naris_in.post.(sites{iSite}).psd.(type{iType,1})],2);
    Naris_in.control.(sites{iSite}).psd.(type{iType,2}) = Naris_in.pre.(sites{iSite}).psd.(type{iType,2});
    idx_low = nearest_idx3(1, Naris_in.control.(sites{iSite}).psd.(type{iType,2}));
    idx_high = nearest_idx3(ii, Naris_in.control.(sites{iSite}).psd.(type{iType,2}));
    C = fit(Naris_in.control.(sites{iSite}).psd.(type{iType,2})(idx_low:idx_high),10*log10(Naris_in.control.(sites{iSite}).psd.(type{iType,1})(idx_low:idx_high)), 'exp2');
    Control.(sites{iSite}) = C(Naris_in.control.(sites{iSite}).psd.(type{iType,2}));
end

%%
wsize = 1024*2;
p_noise = MS_pink_noise(10000000);
% Fs_pink = (10000000 - 1)/
[Pink.Pxx,Pink.F] =  pwelch((p_noise), hanning(wsize), wsize/2, 2*wsize, 2000);
[Pink.White_Pxx,Pink.White_F] =  pwelch(diff(p_noise), hanning(wsize), wsize/2, 2*wsize, 2000);
% 
% 
% hold on
plot(Pink.F, 10*log10(Pink.Pxx), 'm')
% plot(F_white, 10*log10(PX_white), 'b')
% xlim([1 120])

%% temp plot
figure(222); 
phases = PARAMS.Phases;
c_ord = linspecer(4); 
c_ord(5,:) = [.8,.8,.8];
for iSite = 1:length(sites);
    subplot(1,length(sites),iSite)
    title(sites{iSite})
    for iPhase = 1:length(phases)
        hold on
        plot(Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{iType,2}), 10*log10(Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{iType,1})),'color', c_ord(iPhase,:) );
    end
    
    % for ii = [80:20:200];
    % get the intercept 
    l_idx = nearest_idx3(20, Naris_in.control.(sites{iSite}).psd.(type{iType,2})); 
    h_idx = nearest_idx3(40, Naris_in.control.(sites{iSite}).psd.(type{iType,2})); 
    
    y_inter_mean = mean([mean(10*log10(Naris_in.(phases{1}).(sites{iSite}).psd.(type{iType,1})(l_idx:h_idx))), mean(10*log10(Naris_in.(phases{2}).(sites{iSite}).psd.(type{iType,1})(l_idx:h_idx))), ...
        mean(10*log10(Naris_in.(phases{3}).(sites{iSite}).psd.(type{iType,1})(l_idx:h_idx))), mean(10*log10(Naris_in.(phases{4}).(sites{iSite}).psd.(type{iType,1})(l_idx:h_idx)))]);
    
    l_idx = nearest_idx3(20, Pink.(type{iType,2})); 
    h_idx = nearest_idx3(40, Pink.(type{iType,2}));
    
    offset = mean(10*log10(Pink.(type{iType,1})(l_idx:h_idx))) - y_inter_mean;
    
    plot(Pink.(type{iType,2}), 10*log10(Pink.(type{iType,1}))-offset,'-k')
%     plot(Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{iType,2}), Control.(sites{iSite}), '--k')
    
    %     phases{leg_id} = num2str(ii);
    % end
    legend(phases, 'location', 'SouthEast')
    xlim([0 100])
    close all
end
% cfg_fig = []; cfg_fig.pos = [90 386 1256 401];
% SetFigure(cfg_fig, gcf)
% mkdir(PARAMS.inter_dir, 'AOC_fit')
% if isunix
% saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_'  type{iType,1} '_' cfg.id])
% else
%     saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_'  type{iType,1}  '_' cfg.id])
% end
% close all
%% try to get the AOC for the phases relative to the exp_curve fit
% close all
bands = {'low', 'high'};
phases = PARAMS.Phases;
AOC_val = []; AOC_val_con = []; 
for iPhase = 1:length(phases)
    for iSite = 1:length(sites)
        for iband = 1:2
            psd = Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{iType,1});
            F = Naris_in.(phases{iPhase}).(sites{iSite}).psd.(type{iType,2});
            f_idx = find(F > cfg.gamma_freq(iband,1) & F <= cfg.gamma_freq(iband,2));
            
            f_idx_contrast = find(F > cfg.contrast(iband,1) & F <= cfg.contrast (iband,2));
            
            AOC_val_con.(type{iType,1}).([bands{iband} '_con_' num2str(cfg.contrast (iband,1)) '_' num2str(cfg.contrast (iband,2))])(iPhase, iSite) =  trapz(F(f_idx_contrast) ,10*log10(psd(f_idx_contrast)) -Control.(sites{iSite})(f_idx_contrast)); 
            AOC_val.(type{iType,1}).(bands{iband})(iPhase, iSite) =  trapz(F(f_idx) ,10*log10(psd(f_idx)) -Control.(sites{iSite})(f_idx)) ;
        end
    end
end
figure(223)
subplot(4,1,1)
bar(AOC_val.(type{iType,1}).low)
set(gca, 'xticklabel', phases)
title('Low Gamma')
hline(1)
legend(sites)

subplot(4,1,2)
bar(AOC_val.(type{iType,1}).high)
set(gca, 'xticklabel', phases)
title('High Gamma')
hline(1)
legend(sites)

subplot(4,1,3)
bar(AOC_val_con.(type{iType,1}).([bands{1} '_con_' num2str(cfg.contrast (1,1)) '_' num2str(cfg.contrast (1,2))]))
set(gca, 'xticklabel', phases)
title(strrep([bands{1} '_con_' num2str(cfg.contrast (1,1)) '_' num2str(cfg.contrast (1,2))], '_', '-'))


subplot(4,1,4)
bar(AOC_val_con.(type{iType,1}).([bands{2} '_con_' num2str(cfg.contrast (2,1)) '_' num2str(cfg.contrast (2,2))]))
set(gca, 'xticklabel', phases)
title(strrep([bands{2} '_con_' num2str(cfg.contrast (2,1)) '_' num2str(cfg.contrast (2,2))], '_', '-'))
SetFigure([], gcf)

    mkdir(PARAMS.inter_dir, 'AOC_fit')
if isunix
    saveas(gcf, [PARAMS.inter_dir 'AOC_fit/AOC_val_all' type{iType,1}  '_' cfg.id])
    saveas(gcf, [PARAMS.inter_dir 'AOC_fit/AOC_val_all' type{iType,1}  '_' cfg.id], 'png')

else
    saveas(gcf, [PARAMS.inter_dir 'AOC_fit\AOC_val_all'  type{iType,1} '_' cfg.id])
    saveas(gcf, [PARAMS.inter_dir 'AOC_fit\AOC_val_all'  type{iType,1} '_' cfg.id], 'png')
end

Naris_out.ratio.(type{iType,1}) = AOC_val.(type{iType,1});

Naris_out.ratio_con.(type{iType,1}) = AOC_val_con.(type{iType,1});

Naris_out.ratio_labels = sites; 
close all
end
% %% test space
% 
% 
% %% get some psds values for the band of interest
% phases = fieldnames(Naris_in);
% bands = {'low', 'high'};
% for iband = 1:2
%     for iPhase = 1:length(phases)
%         sites = fieldnames(Naris_in.(phases{iPhase}));
%         for iSite = 1:length(sites)
%             w_psd = Naris_in.(phases{iPhase}).(sites{iSite}).psd.White_Pxx;
%             w_F = Naris_in.(phases{iPhase}).(sites{iSite}).psd.White_F;
%             psd = Naris_in.(phases{iPhase}).(sites{iSite}).psd.Pxx;
%             F = Naris_in.(phases{iPhase}).(sites{iSite}).psd.F;
%             
%             f_idx = find(F > cfg.gamma_freq(iband,1) & F <= cfg.gamma_freq(iband,2));
%             f_idx_contrast = find(F > cfg.contrast(iband,1) & F <= cfg.contrast (iband,2));
%             
%             w_f_idx = find(w_F > cfg.gamma_freq(iband,1) & w_F <= cfg.gamma_freq(iband,2));
%             w_f_idx_contrast = find(w_F > cfg.contrast(iband,1) & w_F <= cfg.contrast (iband,2));
%             
%             %% calculate the raw and the ratio scores
%             Naris_out.(phases{iPhase}).(sites{iSite}).psd.raw.(bands{iband}) = nanmean(10*log10(psd(f_idx)));
%             Naris_out.(phases{iPhase}).(sites{iSite}).psd.ratio.(bands{iband}) = mean(trapz(10*log10(psd(f_idx)))) / mean(trapz(10*log10(psd(f_idx_contrast)))) ;
%             Naris_out.(phases{iPhase}).(sites{iSite}).psd.white_raw.(bands{iband}) = nanmean(10*log10(w_psd(w_f_idx)));
%             Naris_out.(phases{iPhase}).(sites{iSite}).psd.white_ratio.(bands{iband})= mean(trapz(10*log10(w_psd(w_f_idx))))/ mean(trapz(10*log10(w_psd(w_f_idx_contrast))));
%             
%             fields = fieldnames(Naris_out.(phases{iPhase}).(sites{iSite}).psd);
%             for ifield = 1:length(fields)
%                 plot_mat.(fields{ifield}).(bands{iband})(iSite, iPhase) = Naris_out.(phases{iPhase}).(sites{iSite}).psd.(fields{ifield}).(bands{iband});
%             end
%         end
%     end
% end
% Naris_in.(phases{iPhase}).(sites{iSite}).psd.plot_mat = plot_mat; % keep the phase by site matrix of each output for simple plotting.
% 
% %% quick check
% if cfg.plot == 1
%     for iband = 1:2
%         figure(iband)
%         for ifield = 1:length(fields)
%             subplot(2,2,ifield)
%             imagesc(plot_mat.(fields{ifield}).(bands{iband}))
%             set(gca, 'ytick', [1:length(sites)], 'ytickLabel', sites, 'xticklabel', phases)
%             xlabel(fields{ifield})
%         end
%     end
% end
% 
% 
% 
% end