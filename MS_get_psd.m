function out = MS_get_psd(cfg_in, data_in)
%MS_get_psd  simply computes the PSD using Welch's average for input data.
%
%Inputs:
%   cfg [struct]: can be blank if you want to use the defaults.
%   data_in [struct]: should be in the form of a csc from LoadCSC

%Outputs:
%   out [struct]: contains the Pxx and F from the PSD 

% EC 2017-05-24
%% default parameters
cfg_def.hann_win_fac = 2; % hanning window multiplier
cfg_def.hann_win = 1024*cfg_def.hann_win_fac; % hanning window size, best to keep these in base 2 for speed. 
cfg_def.whitefilter = 'on';
cfg  = ProcessConfig2(cfg_def, cfg_in);

%% calculate the PSD

[out.Pxx, out.F] = pwelch(data_in.data, hanning(cfg.hann_win), cfg.hann_win/2, 2*cfg.hann_win, data_in.cfg.hdr{1}.SamplingFrequency);
out.cfg_psd = cfg;

if strcmp(cfg.whitefilter, 'on')
    %do the same for the "white filtered" data to remove the 1/f trend
    [out.White_Pxx, out.White_F] = pwelch(diff(data_in.data), hanning(cfg.hann_win), cfg.hann_win/2, 2*cfg.hann_win, data_in.cfg.hdr{1}.SamplingFrequency);
end

%% if needed for checking the PSD
% figure(100)
% plot(out.F, 10*log10(out.Pxx));
% xlim([0 120])
% hold on
% plot(White_F, 10*log10(White_Pxx), 'k');
% xlim([0 120])


