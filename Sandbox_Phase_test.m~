
load('R104_2016_09_24_Data.mat')
load('R104_2016_09_24_Events.mat')

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


cfg_re = [];
cfg_re.d = [-.1 .1];

% Filter the data for the amplitude x_corr
cfg_filter = [];
cfg_filter.f = [45 65];


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
            NAC_raw{iEvt} = restrict(data.post.NAc_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
            NAC_amp{iEvt} = restrict(d_amp.NAc_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
            NAC_filt{iEvt} = restrict(d_filter1.NAc_pot, Re_Events.OFC_pot.post.low.tstart(iEvt), Re_Events.OFC_pot.post.low.tend(iEvt));
            [NAC_PSD.(site{iSite}){iEvt}.pxx,NAC_PSD.(site{iSite}){iEvt}.f]   = pwelch(NAC_raw{iEvt}.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, 2*cfg_psd.hann_win, NAC_raw{iEvt}.cfg.hdr{1}.SamplingFrequency);

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

% create a pool of circ shift values between -10s/Fs:1/Fs:10s

% for iShuf = length(nShuf) :-1:1
%     
%    
%     shuf_coh(iEvt)
% end

%amp xcorr
[ac_all(iEvt,:), lag] = xcov(AMP.OFC_pot{iEvt}.data,NAC_amp{iEvt}.data,100,  'coeff');
lag = lag * 1/AMP.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency;

% phase lag
[Cxy,F] = cpsd(Evts.OFC_pot{iEvt}.data, NAC_raw{iEvt}.data,length(Evts.OFC_pot{iEvt}.data),0,cfg_coh.NFFT,Evts.OFC_pot{iEvt}.cfg.hdr{1}.SamplingFrequency);
coh_spec_phase(iEvt,:)= -angle(Cxy); %higher value means leading. outputs radians
OFC_NAC_phase_off(iEvt,:) =rad2deg(circ_mean(coh_spec_phase(iEvt, nearest_idx(cfg_filter.f(1), F):nearest_idx(cfg_filter.f(2), F))'));

end


figure(10)
subplot(2,2,1)
% plot(fxx, mean(cxx))
shadedErrorBar(fxx, median(cxx), std(cxx)/sqrt(length(cxx)))
xlim([0 100])
vline(cfg_filter.f)
ylabel('Coherence')

subplot(2,2,2)
hold on
plot(F, rad2deg(circ_mean(coh_spec_phase)), 'b');
vline(cfg_filter.f)
xlim([0 100])
plot([0:100], zeros(101), 'k--')
text(5, -5, ['Offset = ' num2str(circ_mean(OFC_NAC_phase_off),2) ' Degrees'], 'fontsize', 12)
ylabel('Phase offset (deg)')

subplot(2,2,3)
plot(lag, mean(ac_all))
vline(0)
% xlabel('lag'

subplot(2,2,4)
hold on
shadedErrorBar(F_PS, median(phase_slopes), std(phase_slopes)/sqrt(length(phase_slopes)))
plot(F_PS, zeros(length(F_PS)), 'k--')
xlim([F_PS(1) F_PS(end)])
ylabel('Phase Slope/Hz')
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
text(5, 20, ['Offset = ' num2str(OFC_NAC_phase_off(iEvt,:),2) ' Degrees'], 'fontsize', 12)


%
subplot(2,3,6)
plot(F_PS, phase_slopes(iEvt,:))
xlim([cfg_filter.f(1) cfg_filter.f(2)])
ylabel('phase/deg')


