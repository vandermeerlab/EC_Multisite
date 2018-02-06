function [phase_slope] = MS_phase_slope(cfg_in, data1, data2)
%% MS_phase_slope: calculates tyhe phase slope between data1 and data two
%       for a specified frequency range in cfg.freq field.  Phase
%       differences are caluculated across the full frequency range using
%       cpsd.  The data is then restricted to the frequencies of interest
%       to speed up the circular regression calculation.  The resulting
%       slopes are then averaged across the cfg.freq range and oupt as
%       "phase_slope".
%
% INPUTS
% - cfg:   [struct] contains configuration parameters.
% - data1: [N x 1]
% - data2: [N x 1]
% OUTPUTS
%
% - phase_slopes: linear regression fit slopes within frequency range
% specified in cfg.freq
%
% CONFIGS
%
% Makes use of circ_regress by MvdM 2018
%
%
% EC 2018
%
%% set up config
cfg_def = [];
cfg_def.window = 2048; % keep in base 2 for speed
cfg_def.nOverlap = 1024;
cfg_def.nFFT = 2^16;
cfg_def.freq = [45 65]; %45-65 for low gamma, 70-90 for high gamma
cfg_def.freq_fraction = .5;
cfg_def.xlim = [1 100];
% cfg for circ_regress
cfg_def.circ_reg.window_size = 11;
cfg_def.debug = 1;


cfg = ProcessConfig2(cfg_def, cfg_in);
if ~isfield(cfg, 'Fs')
    error('Requires a sampleing frequency (cfg.Fs)')
end

%% get teh phase offset from cpsd
[Pxy,F] = cpsd(data1,data2,length(data1),0,cfg.nFFT,cfg.Fs);
phase_diff = angle(Pxy);
phase_diff_plus_pi = phase_diff+pi;

if cfg.debug
    figure(99); subplot(3,3,1); cla;
    h1 = plot(data1); hold on; h2 = plot(data2,'r');
    axis tight;
    
    subplot(3,3,2:3);
    hold on
%     plot(F,phase_diff); set(gca,'XLim',cfg.xlim);
    plot(F,phase_diff_plus_pi); set(gca,'XLim',cfg.xlim);
    rectangle('position', [cfg.freq(1), -pi, (cfg.freq(2)-cfg.freq(1)), .02])
    ylabel('phase difference (rad)');
    
    subplot(3,3,4);
    plot(F,abs(Pxy)); set(gca,'XLim',cfg.xlim);
    ylabel('power');
    drawnow; pause(1);
end
%% extract the frequncies and data to 1Hz increments (for speed)
F_to_get = cfg.freq(1):cfg.freq_fraction:cfg.freq(2);

keep_idx = nearest_idx2(F_to_get,F);
Freq_out = F(keep_idx);
Phase_diff_pi = phase_diff_plus_pi(keep_idx);

% [~,f1_idx] = min(abs(F-(cfg.freq(1))));
% [~,f2_idx] = min(abs(F-(cfg.freq(2))));
% % cut out the F of interest
% % F_interest = F(f1_idx:f2_idx);
% Plus_pi = phase_diff_plus_pi(f1_idx:f2_idx);


%% compute the slope using MVDM's circ_regress function.  This uses a sliding window

% calculate the circ regression for that section.

[slopes,intercepts] = circ_regress(cfg.circ_reg,Phase_diff_pi);
%% check if needed
if cfg.debug
    subplot(3,3,5:6)
    plot(Freq_out, slopes, '-r')
    ylabel('Slopes')
    xlim(cfg.freq)
    win_size = floor(cfg.circ_reg.window_size/2);
    for iS = 1+win_size:length(slopes)-win_size
        % plot
        n1 = iS-win_size;
        n2 = iS+win_size;
        subplot(3,3,8:9)
        plot(Freq_out(iS-win_size:iS+win_size), Phase_diff_pi(iS-win_size:iS+win_size),'.k');
        plot([Freq_out(n1), Freq_out(n2)],[slopes(iS)*-win_size+Phase_diff_pi(iS), slopes(iS)*win_size+Phase_diff_pi(iS)],'-r');
        xlabel('Frequency (Hz)'); ylabel('Phase diff')
        xlim(cfg.freq)
        hold on
        
    end
end
%% get the mean slope values

phase_slope = mean(slopes);

end