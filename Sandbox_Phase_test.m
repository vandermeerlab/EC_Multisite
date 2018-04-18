
% load('R104_2016_09_24_Data.mat')
% load('R104_2016_09_24_Events.mat')
% %%
% clear all
% close all
% load('R107_2017_08_03_Data.mat')
% load('R107_2017_08_03_Events.mat')

%%  x-piri
clear all
close all
MS_initialize
load('R112_2017_08_03_Events.mat')
load('R112_2017_08_03_Data.mat')
cfg_in = [];
% data.post.NAc_pot = data.post.Piri_O_pot;
%% get only the pot sessions
site = fieldnames(data.post);
expStr = '*_pot';
regStr = ['^',strrep(strrep(expStr,'?','.'),'*','.{0,}'),'$'];
starts = regexpi(site, regStr);
iMatch = ~cellfun(@isempty, starts);
idx = find(iMatch);
site = site(idx);

% cfgs
cfg_psd = [];
cfg_psd.hann_win = 256;

%resize
cfg_re = [];
cfg_re.d = [-.1 .1];

% Filter the data for the amplitude x_corr
cfg_filter = [];
cfg_filter.f = [45 65];
cfg_filter2.f = [70 90];


clear Evts PSD FILT AMP NAC d_amp d_filter1
%% get some events
for iSite = 1:length(site)
    % add some padding
    Re_Events.(site{iSite}).post.low = ResizeIV(cfg_re, Events.(site{iSite}).post.low);
    d_filter1.(site{iSite}) = FilterLFP(cfg_filter, data.post.(site{iSite}));
    d_amp.(site{iSite}) = d_filter1.(site{iSite});
    d_amp.(site{iSite}).data = abs(hilbert(d_filter1.(site{iSite}).data));
    
end
for iSite = 1:length(site)
    for iEvt = length(Re_Events.(site{iSite}).post.low.tstart):-1:1
        Evts.(site{iSite}){iEvt} = restrict(data.post.(site{iSite}), Re_Events.(site{iSite}).post.low.tstart(iEvt), Re_Events.(site{iSite}).post.low.tend(iEvt));
        FILT.(site{iSite}){iEvt} = restrict(d_filter1.(site{iSite}), Re_Events.(site{iSite}).post.low.tstart(iEvt), Re_Events.(site{iSite}).post.low.tend(iEvt));
        AMP.(site{iSite}){iEvt} = restrict(d_amp.(site{iSite}), Re_Events.(site{iSite}).post.low.tstart(iEvt), Re_Events.(site{iSite}).post.low.tend(iEvt));
        [PSD.(site{iSite}){iEvt}.pxx,PSD.(site{iSite}){iEvt}.f]   = pwelch(Evts.(site{iSite}){iEvt}.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, 2*cfg_psd.hann_win, Evts.(site{iSite}){iEvt}.cfg.hdr{1}.SamplingFrequency);
        if strcmp(site{iSite}, 'OFC_pot')
%             NAC_raw{iEvt} = restrict(data.post.NAc_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
%             NAC_amp{iEvt} = restrict(d_amp.NAc_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
%             NAC_filt{iEvt} = restrict(d_filter1.NAc_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
%             [NAC_PSD.(site{iSite}){iEvt}.pxx,NAC_PSD.(site{iSite}){iEvt}.f]   = pwelch(NAC_raw{iEvt}.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, 2*cfg_psd.hann_win, NAC_raw{iEvt}.cfg.hdr{1}.SamplingFrequency);
%             
%             
            NAC_raw{iEvt} = restrict(data.post.Piri_O_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
            NAC_amp{iEvt} = restrict(d_amp.Piri_O_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
            NAC_filt{iEvt} = restrict(d_filter1.Piri_O_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
            [NAC_PSD.(site{iSite}){iEvt}.pxx,NAC_PSD.(site{iSite}){iEvt}.f]   = pwelch(NAC_raw{iEvt}.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, 2*cfg_psd.hann_win, NAC_raw{iEvt}.cfg.hdr{1}.SamplingFrequency);

            %             % add in circshifts data
            %             shift_range = -(cfg.shuff_range*NAC_raw{iEvt}.cfg.hdr{1}.SamplingFrequency):cfg.shuff_range*NAC_raw{iEvt}.cfg.hdr{1}.SamplingFrequency;
            %             shifts = datasample(shift_range, cfg.nShuf, 'Replace', false);
            %
            %             for iShuf = length(shifts):-1:1
            %                 NAC_shift{iEvt}{iShuf} = restrict(data.post.NAc_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
            %             end
            
        end
    end
end
%% smae for high gamma amp corr
for iSite = 1:length(site)
    % add some padding
    Re_Events2.(site{iSite}).post.high= ResizeIV(cfg_re, Events.(site{iSite}).post.high);
    d_filter2.(site{iSite}) = FilterLFP(cfg_filter, data.post.(site{iSite}));
    d_amp2.(site{iSite}) = d_filter2.(site{iSite});
    d_amp2.(site{iSite}).data = abs(hilbert(d_filter2.(site{iSite}).data));
    
end
for iSite = 1:length(site)
    for iEvt = length(Re_Events2.(site{iSite}).post.high.tstart):-1:1
        Evts2.(site{iSite}){iEvt} = restrict(data.post.(site{iSite}), Re_Events2.(site{iSite}).post.high.tstart(iEvt), Re_Events2.(site{iSite}).post.high.tend(iEvt));
        FILT2.(site{iSite}){iEvt} = restrict(d_filter2.(site{iSite}), Re_Events2.(site{iSite}).post.high.tstart(iEvt), Re_Events2.(site{iSite}).post.high.tend(iEvt));
        AMP2.(site{iSite}){iEvt} = restrict(d_amp2.(site{iSite}), Re_Events2.(site{iSite}).post.high.tstart(iEvt), Re_Events2.(site{iSite}).post.high.tend(iEvt));
        [PSD2.(site{iSite}){iEvt}.pxx,PSD2.(site{iSite}){iEvt}.f]   = pwelch(Evts2.(site{iSite}){iEvt}.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, 2*cfg_psd.hann_win, Evts2.(site{iSite}){iEvt}.cfg.hdr{1}.SamplingFrequency);
        if strcmp(site{iSite}, 'OFC_pot')
%             NAC_raw2{iEvt} = restrict(data.post.NAc_pot, Re_Events2.OFC_pot.post.high.tstart(iEvt), Re_Events2.OFC_pot.post.high.tend(iEvt));
%             NAC_amp2{iEvt} = restrict(d_amp2.NAc_pot, Re_Events2.OFC_pot.post.high.tstart(iEvt), Re_Events2.OFC_pot.post.high.tend(iEvt));
%             NAC_filt2{iEvt} = restrict(d_filter2.NAc_pot, Re_Events2.OFC_pot.post.high.tstart(iEvt), Re_Events2.OFC_pot.post.high.tend(iEvt));
%             [NAC_PSD2.(site{iSite}){iEvt}.pxx,NAC_PSD2.(site{iSite}){iEvt}.f]   = pwelch(NAC_raw{iEvt}.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, 2*cfg_psd.hann_win, NAC_raw{iEvt}.cfg.hdr{1}.SamplingFrequency);
%            
            NAC_raw2{iEvt} = restrict(data.post.Piri_O_pot, Re_Events2.OFC_pot.post.high.tstart(iEvt), Re_Events2.OFC_pot.post.high.tend(iEvt));
            NAC_amp2{iEvt} = restrict(d_amp.Piri_O_pot, Re_Events2.OFC_pot.post.high.tstart(iEvt), Re_Events2.OFC_pot.post.high.tend(iEvt));
            NAC_filt2{iEvt} = restrict(d_filter2.Piri_O_pot, Re_Events2.OFC_pot.post.high.tstart(iEvt), Re_Events2.OFC_pot.post.high.tend(iEvt));
            [NAC_PSD2.(site{iSite}){iEvt}.pxx,NAC_PSD2.(site{iSite}){iEvt}.f]   = pwelch(NAC_raw2{iEvt}.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, 2*cfg_psd.hann_win, NAC_raw2{iEvt}.cfg.hdr{1}.SamplingFrequency);

%    
        end
    end
end

%% try the PS
tic
for iEvt = length(Evts.OFC_pot):-1:1
    cfg_phase = [];
    cfg_phase.Fs = Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;
    cfg_phase.debug = 0;
    cfg_phase.freq = [30 90];
    cfg_phase.circ_reg.window_size = 9;
    cfg_phase.circ_reg.offset_delta = pi/16;
    [phase_slopes(iEvt,:), F_PS] = MS_phase_slope(cfg_phase, Evts.OFC_pot{iEvt}.data, NAC_raw{iEvt}.data);
    % % %     SetFigure([], gcf)
    % % %     pause(1)
    % % % same for phase lag
    %
    % % coherence
    cfg_coh.spec_window = 128;
    cfg_coh.NFFT = cfg_coh.spec_window*4;
    [cxx(iEvt,:), fxx] = mscohere(Evts.OFC_pot{iEvt}.data, NAC_raw{iEvt}.data, hanning(cfg_coh.spec_window),cfg_coh.spec_window/2, cfg_coh.NFFT, Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency);
    
    % % create a pool of circ shift values between -10s/Fs:1/Fs:10s
    % shift_range = -(cfg.shuff_range*Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency):cfg.shuff_range*Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;
    % shifts = datasample(shift_range, cfg.nShuf, 'Replace', false);
    % for iShuf = length(shifts):-1:1
    %     [shuf_coh(iShuf,:), shift_fxx] = mscohere(Evts.OFC_pot{iEvt}.data, circshift(NAC_raw{iEvt}.data, [shifts(iShuf), 0]), hanning(cfg_coh.spec_window),cfg_coh.spec_window/2, cfg_coh.NFFT, Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency);
    % end
    
    %amp xcorr
    [ac_all(iEvt,:), lag] = xcov(AMP.OFC_pot{iEvt}.data,NAC_amp{iEvt}.data,100,  'coeff');
    lag = lag * 1/AMP.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;
    
    
    % phase lag
    [Cxy,F] = cpsd(Evts.OFC_pot{iEvt}.data, NAC_raw{iEvt}.data,length(Evts.OFC_pot{iEvt}.data),0,cfg_coh.NFFT,Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency);
    coh_spec_phase(iEvt,:)= -angle(Cxy); %higher value means leading. outputs radians
    OFC_NAC_phase_off(iEvt,:) =circ_mean(coh_spec_phase(iEvt, nearest_idx(cfg_filter.f(1), F):nearest_idx(cfg_filter.f(2), F))');
    
end

%% for hg
for iEvt = length(Evts2.OFC_pot):-1:1
    cfg_phase = [];
    cfg_phase.Fs = Evts2.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;
    cfg_phase.debug = 0;
    cfg_phase.freq = [30 90];
    cfg_phase.circ_reg.window_size = 9;
    cfg_phase.circ_reg.offset_delta = pi/16;
    [phase_slopes2(iEvt,:), F_PS2] = MS_phase_slope(cfg_phase, Evts2.OFC_pot{iEvt}.data, NAC_raw2{iEvt}.data);
    % % %     SetFigure([], gcf)
    % % %     pause(1)
    % % % same for phase lag
    %
    % % coherence
    cfg_coh.spec_window = 128;
    cfg_coh.NFFT = cfg_coh.spec_window*4;
    [cxx2(iEvt,:), fxx2] = mscohere(Evts2.OFC_pot{iEvt}.data, NAC_raw2{iEvt}.data, hanning(cfg_coh.spec_window),cfg_coh.spec_window/2, cfg_coh.NFFT, Evts2.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency);
    
    % % create a pool of circ shift values between -10s/Fs:1/Fs:10s
    % shift_range = -(cfg.shuff_range*Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency):cfg.shuff_range*Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;
    % shifts = datasample(shift_range, cfg.nShuf, 'Replace', false);
    % for iShuf = length(shifts):-1:1
    %     [shuf_coh(iShuf,:), shift_fxx] = mscohere(Evts.OFC_pot{iEvt}.data, circshift(NAC_raw{iEvt}.data, [shifts(iShuf), 0]), hanning(cfg_coh.spec_window),cfg_coh.spec_window/2, cfg_coh.NFFT, Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency);
    % end
    
    %amp xcorr
    [ac_all2(iEvt,:), lag2] = xcov(AMP2.OFC_pot{iEvt}.data,NAC_amp2{iEvt}.data,100,  'coeff');
    lag2 = lag2 * 1/AMP2.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;
    
    % phase lag
    [Cxy2,F2] = cpsd(Evts2.OFC_pot{iEvt}.data, NAC_raw2{iEvt}.data,length(Evts2.OFC_pot{iEvt}.data),0,cfg_coh.NFFT,Evts2.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency);
    coh_spec_phase2(iEvt,:)= -angle(Cxy2); %higher value means leading. outputs radians
    OFC_NAC_phase_off2(iEvt,:) =(circ_mean(coh_spec_phase2(iEvt, nearest_idx(cfg_filter2.f(1), F2):nearest_idx(cfg_filter2.f(2), F2))'));
    
end

%%
close all
c_ord = linspecer(4);
figure(10)
subplot(2,2,1)
% plot(fxx, mean(cxx))
hold on
shadedErrorBar(fxx, median(cxx), std(cxx)/sqrt(length(cxx)),'b')
shadedErrorBar(fxx2, median(cxx2), std(cxx2)/sqrt(length(cxx2)),'g')
xlim([0 100])
% vline(cfg_filter.f)
vline([cfg_filter.f,cfg_filter2.f], {'--b', '--b', '--g', '--g'})
% p_95 = prctile(shuf_coh,95, 1);
% plot(shift_fxx, p_95, '--r');
% plot(shift_fxx, p_95, '--r')
xlabel('Frequency (Hz)');
ylabel('Coherence');

subplot(2,2,2)
hold on
% plot(F, rad2deg(circ_mean(coh_spec_phase)), 'b');
shadedErrorBar(F, rad2deg(circ_mean(coh_spec_phase)), rad2deg(circ_std(coh_spec_phase)/sqrt(length(coh_spec_phase))),'b')
shadedErrorBar(F2, rad2deg(circ_mean(coh_spec_phase2)), rad2deg(circ_std(coh_spec_phase2)/sqrt(length(coh_spec_phase2))),'g')
vline([cfg_filter.f,cfg_filter2.f], {'--b', '--b', '--g', '--g'})
% vline(cfg_filter.f)
ylim([-200 200])
xlim([0 100])
plot([0:100], zeros(101), 'k--')
text(50, 195, ['lg:offset = ' num2str(rad2deg(circ_mean(OFC_NAC_phase_off)),3) ' Degrees'], 'fontsize', 12)
text(50,175, ['hg:offset = ' num2str(rad2deg(circ_mean(OFC_NAC_phase_off2)),3) ' Degrees'], 'fontsize', 12)
ylabel('Phase offset (deg)')
xlabel('Frequency (Hz)');

subplot(2,2,3)
% plot(lag, mean(ac_all)) 
hold on
h1= shadedErrorBar(lag, mean((ac_all)), std(ac_all)/sqrt(length(ac_all)),'b');
h2 =shadedErrorBar(lag2, mean(ac_all2), std(ac_all2)/sqrt(length(ac_all2)),'g');
h2.mainLine.Color = c_ord(3,:);
ylim([0 1])
xlabel('lag (ms)')
ylabel('Amplitude x-corr')
y_l = get(gca, 'ylim');
[ac_max, idx] =  max(mean(ac_all));
s_amp = sprintf('lg xcorr = %.2f Lag = %.2fms', ac_max, lag(idx));
text(-0.048, y_l(2)-(y_l(2)*.05), s_amp, 'fontsize', 10);
% high gamma
[ac_max2, idx2] =  max(mean(ac_all2));
s_amp2 = sprintf('hg xcorr = %.2f Lag = %.2fms', ac_max2, lag2(idx2));
text(-0.048, y_l(2)-(y_l(2)*.1), s_amp2, 'fontsize', 10);
% legend('low gamma', 'high gamma')
vline([0,lag2(idx),lag2(idx2)],{'-k','--b', '--g'})


subplot(2,2,4)
hold on
shadedErrorBar(F_PS, median(phase_slopes), std(phase_slopes)/sqrt(length(phase_slopes)),'b')
shadedErrorBar(F_PS2, median(phase_slopes2), std(phase_slopes2)/sqrt(length(phase_slopes2)),'g')
plot(F_PS, zeros(length(F_PS)), 'k--')
xlim([F_PS(1) F_PS(end)])
ylabel('Phase Slope/Hz')
xlabel('Frequency (Hz)');
vline([cfg_filter.f,cfg_filter2.f], {'--b', '--b', '--g', '--g'})
% vline(cfg_filter.f)


SetFigure([], gcf);
toc
%% plot one to be sure
iEvt = 86;
figure(111)
subplot(2,3,1)
hold on
plot(Evts.OFC_pot{iEvt}.tvec, Evts.OFC_pot{iEvt}.data, 'b')
plot(NAC_raw{iEvt}.tvec, NAC_raw{iEvt}.data, 'r')

xlim([Evts.OFC_pot{iEvt}.tvec(1) Evts.OFC_pot{iEvt}.tvec(end)])

subplot(2,3,2)
hold on
plot(PSD.OFC_pot{iEvt}.f, 10*log10(PSD.OFC_pot{iEvt}.pxx));
plot(NAC_PSD.OFC_pot{iEvt}.f, 10*log10(NAC_PSD.OFC_pot{iEvt}.pxx), 'r')
legend({'OFC', 'NAc'})
xlim([0 100])
ylabel('Power')

% get the amp x_corr

[ac, lag] = xcov(AMP.OFC_pot{iEvt}.data,NAC_amp{iEvt}.data,100,  'coeff');
lag = lag * 1/AMP.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;

subplot(2,3,3)
hold on
plot(FILT.OFC_pot{iEvt}.tvec, FILT.OFC_pot{iEvt}.data, 'b')
plot(AMP.OFC_pot{iEvt}.tvec, AMP.OFC_pot{iEvt}.data, 'b--')

plot(NAC_filt{iEvt}.tvec, NAC_filt{iEvt}.data, 'r')
plot(NAC_amp{iEvt}.tvec, NAC_amp{iEvt}.data, 'r--')
xlim([AMP.OFC_pot{iEvt}.tvec(1) AMP.OFC_pot{iEvt}.tvec(end)])
y_l = get(gca, 'ylim');
[ac_max, idx] =  max(ac);
s_amp = sprintf('Amp xcorr = %.2f\nLag = %.2fms', ac_max, lag(idx)*1000);
text(AMP.OFC_pot{iEvt}.tvec(5), y_l(2)-(y_l(2)*.1), s_amp, 'fontsize', 12);



% coherence between

subplot(2,3,4)
plot(fxx, cxx(iEvt,:))
xlim([0 100])
ylim([.5 1])
vline(cfg_filter.f)
ylabel('Coherence')

% get the phase offset value and plot in the corner


subplot(2,3,5)
hold on
plot(F, rad2deg(coh_spec_phase(iEvt,:)), 'b');
vline(cfg_filter.f)
xlim([0 100])
plot([0:100], zeros(101), 'k--')
text(5, 20, ['Offset = ' num2str(rad2deg(OFC_NAC_phase_off(iEvt,:)),3) ' Degrees'], 'fontsize', 12)


%
subplot(2,3,6)
plot(F_PS, phase_slopes(iEvt,:))
xlim([cfg_filter.f(1) cfg_filter.f(2)])
ylabel('phase/deg')

% %% 
% degrees = 0:1:360;
% vales = zeros(1,length(degrees));
% figure(11)
% compass(degrees, vales)
%%
load('phase_off_R107_112.mat')
n = 100;
figure(2222)
subplot(2,2,1)
[t,r] = rose(coh.R107_2017_08_03_OFC_NAc.lg,n);
mean(
h1 = polar(t,r);
set(h1, 'color',c_ord(1,:), 'linewidth', 2)
view([90 -90])


subplot(2,2,3)
[t2,r2] = rose(coh.R112_2017_08_03_xpiri.lg,n);
h3 = polar(t2,r2);
set(h3, 'color',c_ord(4,:), 'linewidth', 2)
view([90 -90])


subplot(222)
[t,r] = rose(coh.R107_2017_08_03_OFC_NAc.hg,n);
h2 = polar(t,r);
set(h2, 'color',c_ord(3,:), 'linewidth', 2)
view([90 -90])
% hold on
subplot(2,2,4)
[t2,r2] = rose(coh.R112_2017_08_03_xpiri.hg,n);
h4 = polar(t2,r2);
set(h4, 'color',c_ord(4,:), 'linewidth', 2)

view([90 -90])

SetFigure([], gcf)
% tightfig