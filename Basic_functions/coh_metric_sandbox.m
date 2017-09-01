%% sandbox for coherence

% generate two sine waves
clear d d_int d_comp temp temp_2
g_freq = [45 65];
% wsize = 2^8;
freq=60;
amp=0;
ts=2000;
t_l=.25;
p_off = 10;
sd = 6;
noise = 1;
d.tvec=0:1/ts:t_l;
wsize = length(d.tvec);
d.data(1,:) = sin(2*pi*freq*d.tvec);
d.data(2,:) = sin(2*pi*freq*d.tvec+degtorad(p_off));
d.data(3,:) = d.data(1,:).*gausswin(length(d.data(1,:)),sd)';
d.data(4,:) = d.data(2,:).*gausswin(length(d.data(2,:)),sd)';

% create a more complex signal
for ii = 1:4
    d.data(ii,:) = d.data(ii,:) +noise*randn(size(d.tvec));
end
freqs = [4 32 80 112 440];
for ifreq =1:length(freqs) 
temp(ifreq,:) = sin(2*pi*freqs(ifreq)*d.tvec);
temp_2(ifreq,:) = sin(2*pi*freqs(ifreq)*d.tvec+degtorad(p_off));
end
d.data(1,:) = sum(temp);
d.data(2,:) = sum(temp_2);



figure(888)
subplot(2,2,1)
plot(d.tvec, d.data(1:2,:))
xlim([d.tvec(1) d.tvec(end)])
legend({'orig', ['phase off by ' num2str(p_off)]})
subplot(2,2,3)
plot(d.tvec, d.data(3:4,:))
xlim([d.tvec(1) d.tvec(end)])
legend({'orig', ['phase off by ' num2str(p_off)]})
%% try cohere with mscoh
[c, f] = mscohere(d.data(1,:), d.data(2,:),hanning(wsize),wsize/2,length(d.tvec),d.cfg.hdr{1}.SamplingFrequency);

f_idx = find(f > g_freq(1) & f <= g_freq(2));
ms_coh_mean = mean(c(f_idx));

%% same thing but interpolated by 10
clear d_int
int_ts = 20000;
for ii = 1:4
    d_int.data(ii,:) = interp1(d.tvec, d.data(ii,:), 0:1/int_ts:t_l);
end
d_int.tvec = 0:1/int_ts:t_l;
subplot(2,2,2)
plot(d_int.tvec, d_int.data(1:2,:))
legend({'orig', ['phase off by ' num2str(p_off)]})
subplot(2,2,4)
plot(d_int.tvec, d_int.data(3:4,:))
legend({'orig', ['phase off by ' num2str(p_off)]})

[c_i, f_i] = mscohere(d_int.data(1,:), d_int.data(2,:),hanning(wsize),wsize/2,2*wsize,int_ts);

f_idx = find(f_i > g_freq(1) & f_i <= g_freq(2));
ms_coh_mean_int = mean(c_i(f_idx));

%% try with the fieldtrip method for both the gaussian and the normal
addpath('/Users/jericcarmichael/Documents/GitHub/fieldtrip')
ft_defaults
figure(999)

for jj = 1:2
    % conver the data to fieldtrip format
    d.label = {'c', ['c' num2str(p_off)], 'g', ['g' num2str(p_off)]};
    d_int.label = d.label;
    if jj ==1
        chans = {'c', ['c' num2str(d.cfg.phase_off)]};
        nSub = 2;
        ft_d = TSDtoFT([], d);
    elseif jj ==2
        chans = {'c', ['c' num2str(d.cfg.phase_off)]};
        nSub = 3;
        ft_d = TSDtoFT([], d_int);
    end
    % extract the frequency components
    cfg            = [];
    cfg.output     = 'fourier';
    cfg.method     = 'mtmfft';
    cfg.foi        = 5:0.5:120;
    cfg.tapsmofrq  = 5;
    cfg.keeptrials = 'yes';
    freqfourier    = ft_freqanalysis(cfg, ft_d);
    
    % compute the coherence
    cfg            = [];
    cfg.method     = 'coh';
    cfg.channel    = chans;
    cfg.channelcmb = chans;
    fdfourier      = ft_connectivityanalysis(cfg, freqfourier);
    
    % plot the coherence
    cfg                  = [];
    cfg.parameter        = 'cohspctrm';
    cfg.xlim             = [5 120];
    cfg.refchannel       = chans;
    subplot(1,3,nSub);
    ft_singleplotER(cfg, fdfourier);
    ylim([0.4 1])
    if jj == 1
        f_idx = find(fdfourier.freq > g_freq(1) & fdfourier.freq <= g_freq(2));
        ft_coh_mean = mean(fdfourier.cohspctrm(1,2,f_idx));
    elseif jj == 2
        f_idx = find(fdfourier.freq > g_freq(1) & fdfourier.freq <= g_freq(2));
        ft_coh_mean_int = mean(fdfourier.cohspctrm(1,2,f_idx));
    end
end

%% temp interp of the gen signal
win_size = length(pairs_in.(pairs{iPair}).(main{1}){iEvt}.data);

[Coh_evt_temp.c{iEvt}, Coh_evt_temp.f{iEvt}] = mscohere(pairs_in.(pairs{iPair}).(main{1}){iEvt}.data,...
pairs_in.(pairs{iPair}).(main{2}){iEvt}.data,...
hanning(round(win_size/2)),round(win_size/4),win_size,...
pairs_in.(pairs{iPair}).(main{1}){iEvt}.cfg.hdr{1}.SamplingFrequency);

figure(222)
plot(Coh_evt_temp.f{iEvt}, Coh_evt_temp.c{iEvt})
xlim([0 120])

figure(333)
plot(pairs_in.(pairs{iPair}).(main{1}){iEvt}.tvec, pairs_in.(pairs{iPair}).(main{1}){iEvt}.data, pairs_in.(pairs{iPair}).(main{2}){iEvt}.tvec, pairs_in.(pairs{iPair}).(main{2}){iEvt}.data)

% interp
fs_2 = 20000; 
di.tvec = pairs_in.(pairs{iPair}).(main{1}){iEvt}.tvec(1):1/fs_2:pairs_in.(pairs{iPair}).(main{1}){iEvt}.tvec(end); 
di.data(1,:) = interp1(pairs_in.(pairs{iPair}).(main{1}){iEvt}.tvec, pairs_in.(pairs{iPair}).(main{1}){iEvt}.data, di.tvec);
di.data(2,:) = interp1(pairs_in.(pairs{iPair}).(main{2}){iEvt}.tvec, pairs_in.(pairs{iPair}).(main{2}){iEvt}.data, di.tvec);

win_size = length(di.tvec);
[temp.c, temp.f] = mscohere(di.data(1,:), di.data(2,:),...
hanning(round(win_size/4)),round(win_size/8),win_size,fs_2);

figure(222)
hold on
plot(temp.f, temp.c)
xlim([0 120])

%% print the results
figure(999)
subplot(1,3,1)
plot(f, c,f_i, c_i)
ylim([0 1])
legend({'normal', 'interp'})
xlim([0 120])
% fprintf(['\n Standard mscohere: ' num2str(ms_coh_mean)]);
% fprintf(['\n Interpol mscohere: ' num2str(ms_coh_mean_int)]);
% fprintf(['\n FT_stand mscohere: ' num2str(ft_coh_mean)]);
% fprintf(['\n FT_inter mscohere: ' num2str(ft_coh_mean_int)]);
% fprintf('\n')


%% temp interp of the gen signal
cfg_lfp = [];
cfg_lfp.len = 2000;
cfg_lfp.freq = 50;
cfg_lfp.noise_val = .5;
cfg_lfp.type = 'complex';
cfg_lfp.complex_freq = [55 83 4];
d = Gen_LFP(cfg_lfp);
d.cfg
%%
% win_size = 2^12;
win_size = length(d.tvec);
[Coh_evt_temp.c{iEvt}, Coh_evt_temp.f{iEvt}] = mscohere(d.data(1,:), d.data(2,:),...
    hamming(round(win_size/2)),round(win_size/4),win_size,d.cfg.hdr{1}.SamplingFrequency);

% interp
clear di
fs_2 = 4000; 
di.tvec = d.tvec(1):1/fs_2:d.tvec(end); 
di.data(1,:) = interp1(d.tvec, d.data(1,:), di.tvec);
di.data(2,:) = interp1(d.tvec, d.data(2,:), di.tvec);
di.label{1} = d.label{1};
di.label{2} = d.label{2};

% win_size = length(di.tvec);
[temp.c, temp.f] = mscohere(di.data(1,:), di.data(2,:),...
hamming(round(win_size/2)),round(win_size/4),round(win_size),fs_2);

figure(9999)
plot(Coh_evt_temp.f{iEvt}, Coh_evt_temp.c{iEvt}, temp.f, temp.c)
xlim([0 120])
% try it in FT

d_ft = TSDtoFT([], di);


    cfg            = [];
    cfg.output     = 'fourier';
    cfg.method     = 'mtmfft';
    cfg.foi        = 5:1:120;
    cfg.tapsmofrq  = 5;
    cfg.keeptrials = 'yes';
    freqfourier    = ft_freqanalysis(cfg, d_ft);
    
    % compute the coherence
    cfg            = [];
    cfg.method     = 'coh';
    cfg.channel    = {'CSC1', 'CSC2'};
    cfg.channelcmb = {'CSC1', 'CSC2'};
    fdfourier      = ft_connectivityanalysis(cfg, freqfourier);
    
    % plot the coherence
    cfg                  = [];
    cfg.parameter        = 'cohspctrm';
%     cfg.xlim             = [1 120];
    cfg.refchannel       = {'CSC1', 'CSC2'};
    figure(221)
    ft_singleplotER(cfg, fdfourier);


%% amplitude coherence

% try it with amplitude filter
cfg_filter = [];
cfg_filter.f = [35 45];
d_f = FilterLFP(cfg_filter, d);

% win_size = length(d_f.tvec);
win_size = 2^12;
[Coh_amp_temp.c, Coh_amp_temp.f] = mscohere(d_f.data(1,:), d_f.data(2,:),...
    hamming(round(win_size/2)),round(win_size/4),win_size,d_f.cfg.hdr{1}.SamplingFrequency);
figure(1000)
subplot(2,1,1)
plot(d_f.tvec, d_f.data)
subplot(2,1,2)
plot(Coh_amp_temp.f, Coh_amp_temp.c)
xlim([0 120])