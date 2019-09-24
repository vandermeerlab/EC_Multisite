function [amp_gram_out] = MS_plot_amp_xcorrogram(cfg_in, data)
%%
%
%
%
%
%
%
%
%% set up configurations
global PARAMS
cfg_def = [];
cfg_def.c_ord = linspecer(4); % get the four naris phase colours
cfg_def.moving_avg = '1D';  % use a 1d convelution to smooth the data by cfg.m_avg_nSamples samples
cfg_def.m_avg_nSample = 60;
cfg_def.freq_m_avg_nSample = 1;

cfg_def.cfg_filter = [];
cfg_def.cfg_filter.freq = [3 100];
cfg_def.cfg_filter.freq_step = 1;

cfg_def.linewidth = 3;
cfg_def.this_pair = 'OFC_NAc';
cfg_def.that_pair = 'OFC_CG';
cfg_def.pre_x_lim = [338 339];
cfg_def.ipsi_x_lim = [854 855];
cfg_def.Fs_restrict = 0.01; % dt to be used when computed the insets.
% cfgs for amp-corr-ogram
cfg_def.cfg_amp_cor_ogram = [];
cfg_def.cfg_amp_cor_ogram.dT = .5;

cfg = ProcessConfig2(cfg_def, cfg_in);

%% used for simple plots under the amp xcorr
%switch NaN's to zero before filtering
data.pre.OFC_pot.data(isnan(data.pre.OFC_pot.data)) = 0;
data.pre.NAc_pot.data(isnan(data.pre.NAc_pot.data)) = 0;
data.pre.CG_pot.data(isnan(data.pre.CG_pot.data)) = 0;

cfg_filter = [];
cfg_filter.f = [45 65];

FILT_OFC_p = FilterLFP(cfg_filter, data.pre.OFC_pot);
FILT_NAc_p = FilterLFP(cfg_filter, data.pre.NAc_pot);
FILT_CG_p = FilterLFP(cfg_filter, data.pre.CG_pot);

AMP_OFC_p = FILT_OFC_p;
AMP_NAc_p = FILT_NAc_p;
AMP_CG_p = FILT_CG_p;

AMP_OFC_p.data = abs(hilbert(FILT_OFC_p.data));
AMP_NAc_p.data = abs(hilbert(FILT_NAc_p.data));
AMP_CG_p.data = abs(hilbert(FILT_CG_p.data));


data.ipsi.OFC_pot.data(isnan(data.ipsi.OFC_pot.data)) = 0;
data.ipsi.NAc_pot.data(isnan(data.ipsi.NAc_pot.data)) = 0;
data.ipsi.CG_pot.data(isnan(data.ipsi.CG_pot.data)) = 0;

FILT_OFC_i = FilterLFP(cfg_filter, data.ipsi.OFC_pot);
FILT_NAc_i = FilterLFP(cfg_filter, data.ipsi.NAc_pot);
FILT_CG_i = FilterLFP(cfg_filter, data.ipsi.CG_pot);

AMP_OFC_i = FILT_OFC_i;
AMP_NAc_i = FILT_NAc_i;
AMP_CG_i = FILT_CG_i;

AMP_OFC_i.data = abs(hilbert(FILT_OFC_i.data));
AMP_NAc_i.data = abs(hilbert(FILT_NAc_i.data));
AMP_CG_i.data = abs(hilbert(FILT_CG_i.data));
%% Get the amp xcorr at each frequency

% setup for AMP xcorr across all frequency bands through time


% set up an NaN array
Freq_list = cfg.cfg_filter.freq(1): cfg.cfg_filter.freq_step:cfg.cfg_filter.freq(end);


for iPhase = 1:length(PARAMS.Phases)
    % Define a time vector for the current phase of the naris protocol.
    times = data.(PARAMS.Phases{iPhase}).OFC_pot.tvec(1):cfg.cfg_amp_cor_ogram.dT:data.(PARAMS.Phases{iPhase}).OFC_pot.tvec(end);
    
    for iF = 1:length(Freq_list)
        
        this_F = Freq_list(iF);
        
        if iF>1
            while last_string_l>0
                fprintf('\b'); % delete previous counter display
                last_string_l = last_string_l - 1;
            end
        end
        last_string_l =fprintf('%0.2f', this_F);
        %%
        cfg_filter= [];
        cfg_filter.display_filter = 0;
        cfg_filter.f = [this_F-1 this_F+1];
        cfg_filter.order = 4;
        cfg_filter.verbose = 0;
        cfg_filter.type = 'fdesign';%'cheby1';
        
        FILT_OFC = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).OFC_pot);
        FILT_NAc = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).NAc_pot);
        FILT_CG = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).CG_pot);
        
        
        AMP_OFC = FILT_OFC;
        AMP_NAc = FILT_NAc;
        AMP_CG = FILT_CG;
        
        AMP_OFC.data = abs(hilbert(FILT_OFC.data));
        AMP_NAc.data = abs(hilbert(FILT_NAc.data));
        AMP_CG.data = abs(hilbert(FILT_CG.data));
        
        % create a Amp-xcorr-ogram
        
        for iT = (1/cfg.cfg_amp_cor_ogram.dT):length(times)-(1/cfg.cfg_amp_cor_ogram.dT)
            
            D_OFC_temp = restrict(AMP_OFC, times(iT)-(1/cfg.cfg_amp_cor_ogram.dT), times(iT+1)+(1/cfg.cfg_amp_cor_ogram.dT));
            D_NAc_temp = restrict(AMP_NAc, times(iT)-(1/cfg.cfg_amp_cor_ogram.dT), times(iT+1)+(1/cfg.cfg_amp_cor_ogram.dT));
            D_CG_temp = restrict(AMP_CG, times(iT)-(1/cfg.cfg_amp_cor_ogram.dT), times(iT+1)+(1/cfg.cfg_amp_cor_ogram.dT));
            
            % run it for the OFC-NAc pair
            [ac_ON, lag_ON] = xcov(D_OFC_temp.data,D_NAc_temp.data,200,  'coeff');
            max_ac_ON = ac_ON(lag_ON==0);
            amp_ac.OFC_NAc.(PARAMS.Phases{iPhase})(iF, iT) = max_ac_ON;
            amp_lag.OFC_NAc.(PARAMS.Phases{iPhase})(iF, iT) = times(iT);
            
            % run it for the OFC-CG pair
            [ac_OC, lag_OC] = xcov(D_OFC_temp.data,D_CG_temp.data,200,  'coeff');
            max_ac_OC = ac_OC(lag_OC==0);
            amp_ac.OFC_CG.(PARAMS.Phases{iPhase})(iF, iT) = max_ac_OC;
            amp_lag.OFC_CG.(PARAMS.Phases{iPhase})(iF, iT) = times(iT);
            
            % run it for the NAc-CG pair ( not used in
            % figure but good to see)
            [ac_NC, lag_NC] = xcov(D_NAc_temp.data,D_CG_temp.data,200,  'coeff');
            max_ac_NC = ac_NC(lag_NC==0);
            amp_ac.NAc_CG.(PARAMS.Phases{iPhase})(iF, iT) = max_ac_NC;
            amp_lag.NAc_CG.(PARAMS.Phases{iPhase})(iF, iT) = times(iT);
            
        end
    end
end


%% get the amp xcorr across all frequencies at each time block.

%% append data
% amp_ac = amp_in.amp.ac;
% amp_lag = amp_in.amp.lag;
% Freq_list = amp_in.amp.f;
Fs  = cfg.cfg_amp_cor_ogram.dT;
%%
for iRun = 1:2
    
    if iRun==1
        this_pair = cfg.this_pair;
    else
        this_pair = cfg.that_pair;
    end
    t_pre =  0:Fs:length(amp_ac.(this_pair).pre)*Fs;
    t_pre = t_pre(1:end);
    t_ipsi = t_pre(end):Fs: (length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi))*Fs;
    t_ipsi = t_ipsi(2:end);
    t_contra = t_ipsi(end):Fs: (length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra))*Fs;
    t_contra = t_contra(2:end);
    t_post = t_contra(end):Fs: (length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra)+length(amp_ac.(this_pair).post))*Fs;
    t_post = t_post(2:end-1);
    
    all_time = [t_pre, t_ipsi, t_contra, t_post];
    
    all_sess_og = [amp_ac.(this_pair).pre, amp_ac.(this_pair).ipsi, amp_ac.(this_pair).contra, amp_ac.(this_pair).post];
    all_sess_lag_og = [amp_lag.(this_pair).pre, amp_lag.(this_pair).ipsi, amp_lag.(this_pair).contra, amp_lag.(this_pair).post];
    
    
    if cfg.moving_avg
        if strcmp(cfg.moving_avg, '1D')
            for ii = size(all_sess_og, 1):-1:1
                all_sess(ii,:) = conv(all_sess_og(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                all_sess_lag(ii,:) = conv(all_sess_lag_og(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                
            end
        elseif strcmp(cfg.moving_avg, '1D_freq')
            for ii = size(all_sess_og, 2):-1:1
                all_sess(:,ii) = conv(all_sess_og(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            end
            
            for jj = size(all_sess_lag_og, 2):-1:1
                all_sess_lag(:,jj) = conv(all_sess_lag_og(:,jj), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            end
        elseif strcmp(cfg.moving_avg, '2D')
            
            all_sess = conv2(all_sess_og, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            all_sess_lag = conv2(all_sess_lag_og, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            
        elseif strcmp(cfg.moving_avg, '1D+1D_freq')
            % do the frequency first
            for ii = size(all_sess_og, 2):-1:1
                all_sess_t(:,ii) = conv(all_sess_og(:,ii), ones(cfg.freq_m_avg_nSample,1)/cfg.freq_m_avg_nSample, 'same');
            end
            
            for jj = size(all_sess_lag_og, 2):-1:1
                all_sess_lag_t(:,jj) = conv(all_sess_lag_og(:,jj), ones(cfg.freq_m_avg_nSample,1)/cfg.freq_m_avg_nSample, 'same');
            end
            % then do the time smoothing
            
            for ii = size(all_sess_og, 1):-1:1
                all_sess(ii,:) = conv(all_sess_t(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                all_sess_lag(ii,:) = conv(all_sess_lag_t(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            end
        else
            all_sess = all_sess_og;
            all_sess_lag = all_sess_lag_og;
            
        end
        
    end
    
    
    %% plot the data
    figure(111)
    ax1 = imagesc('XData', all_time, 'CData', all_sess);
    axis xy
    set(gca, 'ytick',[1 size(all_sess,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (0:500:length(all_time)), 'xticklabel', [0:500:2500])
    
    set(ax1, 'AlphaData', ~isnan(all_sess))
    %     title(['Coherogram: ' S{1} ' ' S{2}] )
    xlabel('time (s)')
    ylabel('frequency')
    %         caxis([0.2 1])
    
    set(gca,'FontSize',20);
    xlabel('time (s)'); ylabel('Frequency (Hz)');
    yL = get(gca, 'ylim');
    % add lines to distinguish phases
    ylim([-4.5 yL(2)]);
    x_l = xlim;
    xlim([0 x_l(end)])
    hold on
    h(1) =rectangle('position', [0 -4.5 length(amp_ac.(this_pair).pre)  5.5], 'facecolor',cfg.c_ord(1,:), 'edgecolor', cfg.c_ord(1,:));
    h(2) =rectangle('position', [length(amp_ac.(this_pair).pre) -4.5 length(amp_ac.(this_pair).ipsi) 5.5], 'facecolor',cfg.c_ord(2,:), 'edgecolor', cfg.c_ord(2,:));
    h(3) =rectangle('position', [(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)) -4.5 length(amp_ac.(this_pair).contra) 5.5], 'facecolor',cfg.c_ord(3,:), 'edgecolor', cfg.c_ord(3,:));
    h(4) = rectangle('position', [(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra)) -4.5 length(amp_ac.(this_pair).post) 5.5], 'facecolor',cfg.c_ord(4,:), 'edgecolor', cfg.c_ord(4,:));
    
    text(length(amp_ac.(this_pair).pre)/2, -1.5, 'pre', 'fontsize', 20, 'horizontalAlignment', 'center')
    text((length(amp_ac.(this_pair).pre)+(length(amp_ac.(this_pair).ipsi)/2)), -1.5, 'ipsi', 'fontsize', 20, 'horizontalAlignment', 'center')
    text((length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+(length(amp_ac.(this_pair).contra)/2)), -1.5, 'contra', 'fontsize', 20, 'horizontalAlignment', 'center')
    text((length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).ipsi)+(length(amp_ac.(this_pair).post)/2)), -1.5, 'post', 'fontsize', 20, 'horizontalAlignment', 'center')
    
    
    %     plot(Tp(1):Tp(end), zeros(length(Tp(1):Tp(end)),1)-off_set,'color', cfg.c_ord(1,:), 'linewidth', 1); % pre
    % %     plot(Tp(end):(Tp(end)+Ti(end)), zeros(length(Tp(end):(Tp(end)+Ti(end))),1)-off_set,'color', cfg.c_ord(2,:), 'linewidth', 10); % ipsi
    %     plot((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end)), zeros(length((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))),1)-off_set,'color', cfg.c_ord(3,:), 'linewidth', 10); % contra
    %     plot((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end)), zeros(length((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))),1)-off_set,'color', cfg.c_ord(4,:), 'linewidth', 10); % post
    %     h = vline([length(amp_ac.(this_pair).pre),(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)),(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).post))], {'k', 'k', 'k'} );
    %     if cfg.lfp_on
    %        plot(D1.tvec_norm,(D1.data*10000)+20, 'color', [1 1 1])
    %        plot(D2.tvec_norm, (D2.data*10000)+20, 'color', [0.5 1 1])
    %     end
    ax =gca;
    ax.Clipping = 'off';
    set(h(:), 'linewidth', cfg.linewidth);
    colormap('PARULA');
    
        text(median(cfg.pre_x_lim)*1, 102, 'V', 'fontsize', 28, 'color', cfg.c_ord(1,:))
    text(median(cfg.ipsi_x_lim)*1, 102, 'V', 'fontsize', 28, 'Color', cfg.c_ord(2,:))
    
    
    h = vline([length(amp_ac.(this_pair).pre), (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)),...
        (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra))], {'k','k', 'k', 'k'}) ;
    set(h(:), 'linewidth', 3)
    
    cfg_fig = [];
    cfg_fig.ft_size =30;
    hc = colorbar;
    caxis([-1 1])
    SetFigure(cfg_fig,gcf);
    
    set(hc,'YTick',[-1 1], 'yticklabel', [-1 1])
    set(gcf, 'position', [0 50 1800*.9 480*.9]);
    
    
    %% save
    mkdir(PARAMS.inter_dir, 'Amp_xcorrogram')
    if isunix
        save_dir = [PARAMS.inter_dir, 'Amp_xcorrogram/'];
    else
        save_dir = [PARAMS.inter_dir, 'Amp_xcorrogram\'];
    end
    if iRun == 1
        fname = [cfg.this_pair '_Amp_xcorrogram'];
    else
        fname = [cfg.that_pair '_Amp_xcorrogram'];
    end
    
    
    % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    saveas(gcf, fname, 'png')
    
    %     saveas(gcf, savefig
    pushdir(save_dir)
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir
    close all
    
    %     %% try the amp max ac lag values.
    %     figure(112)
    %     ax1 = imagesc(all_sess_lag);
    %     axis xy
    %     set(gca, 'ytick',[1 size(all_sess,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)), 'xticklabel', [0:500:2500])
    %
    %     %     set(ax1, 'AlphaData', ~isnan(all_sess))
    %     %     title(['Coherogram: ' S{1} ' ' S{2}] )
    %     xlabel('time (s)')
    %     ylabel('frequency')
    %     %         caxis([0.2 1])
    %
    %     set(gca,'FontSize',20);
    %     xlabel('time (s)'); ylabel('Frequency (Hz)');
    %     yL = get(gca, 'ylim');
    %     % add lines to distinguish phases
    %     ylim([-4.5 yL(2)]);
    %     hold on
    %     h(1) =rectangle('position', [0 -4.5 length(amp_ac.(this_pair).pre)  5.5], 'facecolor',cfg.c_ord(1,:), 'edgecolor', cfg.c_ord(1,:));
    %     h(2) =rectangle('position', [length(amp_ac.(this_pair).pre) -4.5 length(amp_ac.(this_pair).ipsi) 5.5], 'facecolor',cfg.c_ord(2,:), 'edgecolor', cfg.c_ord(2,:));
    %     h(3) =rectangle('position', [(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)) -4.5 length(amp_ac.(this_pair).contra) 5.5], 'facecolor',cfg.c_ord(3,:), 'edgecolor', cfg.c_ord(3,:));
    %     h(4) = rectangle('position', [(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra)) -4.5 length(amp_ac.(this_pair).post) 5.5], 'facecolor',cfg.c_ord(4,:), 'edgecolor', cfg.c_ord(4,:));
    %
    %     text(length(amp_ac.(this_pair).pre)/2, -1.5, 'pre', 'fontsize', 20, 'horizontalAlignment', 'center')
    %     text((length(amp_ac.(this_pair).pre)+(length(amp_ac.(this_pair).ipsi)/2)), -1.5, 'ipsi', 'fontsize', 20, 'horizontalAlignment', 'center')
    %     text((length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+(length(amp_ac.(this_pair).contra)/2)), -1.5, 'contra', 'fontsize', 20, 'horizontalAlignment', 'center')
    %     text((length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).ipsi)+(length(amp_ac.(this_pair).post)/2)), -1.5, 'post', 'fontsize', 20, 'horizontalAlignment', 'center')
    %
    %
    %     %     plot(Tp(1):Tp(end), zeros(length(Tp(1):Tp(end)),1)-off_set,'color', cfg.c_ord(1,:), 'linewidth', 1); % pre
    %     % %     plot(Tp(end):(Tp(end)+Ti(end)), zeros(length(Tp(end):(Tp(end)+Ti(end))),1)-off_set,'color', cfg.c_ord(2,:), 'linewidth', 10); % ipsi
    %     %     plot((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end)), zeros(length((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))),1)-off_set,'color', cfg.c_ord(3,:), 'linewidth', 10); % contra
    %     %     plot((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end)), zeros(length((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))),1)-off_set,'color', cfg.c_ord(4,:), 'linewidth', 10); % post
    %     %     h = vline([length(amp_ac.(this_pair).pre),(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)),(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).post))], {'k', 'k', 'k'} );
    %     %     if cfg.lfp_on
    %     %        plot(D1.tvec_norm,(D1.data*10000)+20, 'color', [1 1 1])
    %     %        plot(D2.tvec_norm, (D2.data*10000)+20, 'color', [0.5 1 1])
    %     %     end
    %     ax =gca;
    %     ax.Clipping = 'off';
    %     set(h(:), 'linewidth', cfg.linewidth);
    %     colormap('PARULA');
    %
    %
    %     h = vline([length(amp_ac.(this_pair).pre), (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)),...
    %         (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra))], {'k','k', 'k', 'k'}) ;
    %     set(h(:), 'linewidth', 3)
    %
    %     cfg_fig = [];
    %     cfg_fig.ft_size =30;
    %     hc = colorbar;
    %     caxis([-0.015 0.015])
    %     SetFigure(cfg_fig,gcf);
    %
    %     %         set(hc,'YTick',[-1 1], 'yticklabel', [-1 1])
    %     set(gcf, 'position', [0 50 1800*.9 480*.9]);
    
    %     %% save
    %
    %     mkdir(PARAMS.inter_dir, 'Amp_xcorrogram')
    %     if isunix
    %         save_dir = [PARAMS.inter_dir, 'Amp_xcorrogram/'];
    %     else
    %         save_dir = [PARAMS.inter_dir, 'Amp_xcorrogram\'];
    %     end
    %     if iRun == 1
    %         fname = [cfg.this_pair '_Amp_lag_ogram'];
    %     else
    %
    %         fname = [cfg.that_pair '_Amp_lag_ogram'];
    %     end
    %
    %
    %     % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
    %     h1 = get(gcf);
    %     D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    %     set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %     saveas(gcf, fname, 'png')
    %
    %     %     saveas(gcf, savefig
    %     pushdir(save_dir)
    %     eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    %     popdir
    %     close all
    
    %% zoom in on the some specifed sections
    % run a tighter amp_xcorr verion on this small section.
    
    
    
    OFC_data_cut = restrict(data.pre.OFC_pot, cfg.pre_x_lim(1)-1.5, cfg.pre_x_lim(2)+1.5);
    NAc_data_cut = restrict(data.pre.NAc_pot, cfg.pre_x_lim(1)-1.5, cfg.pre_x_lim(2)+1.5);
    CG_data_cut = restrict(data.pre.CG_pot, cfg.pre_x_lim(1)-1.5, cfg.pre_x_lim(2)+1.5);
    
    % times_pre = OFC_data_cut.tvec;
times_pre =  cfg.pre_x_lim(1)-1:cfg.Fs_restrict:cfg.pre_x_lim(2)+1;    % cycle across frequencies
    for iF = 1:length(Freq_list)
        
        this_F = Freq_list(iF);
        
        
        if iF>1
            while last_string_l>0
                fprintf('\b'); % delete previous counter display
                last_string_l = last_string_l - 1;
            end
        end
        last_string_l =fprintf('%0.2f', this_F);
        %%
        cfg_filter= [];
        cfg_filter.display_filter = 0;
        cfg_filter.f = [this_F-1 this_F+1];
        cfg_filter.order = 4;
        cfg_filter.verbose = 0;
        cfg_filter.type = 'fdesign';%'cheby1';
        
        FILT_OFC = FilterLFP(cfg_filter, OFC_data_cut);
        FILT_NAc = FilterLFP(cfg_filter, NAc_data_cut);
        FILT_CG = FilterLFP(cfg_filter, CG_data_cut);
        
        
        AMP_OFC_restrict = FILT_OFC;
        AMP_NAc_restrict = FILT_NAc;
        AMP_CG_restrict = FILT_CG;
        
        AMP_OFC_restrict.data = abs(hilbert(FILT_OFC.data));
        AMP_NAc_restrict.data = abs(hilbert(FILT_NAc.data));
        AMP_CG_restrict.data = abs(hilbert(FILT_CG.data));
        
        % create a Amp-xcorr-ogram
        
        for iT = 1:length(times_pre)
            
            D_OFC_temp = restrict(AMP_OFC_restrict, times_pre(iT)-1, times_pre(iT)+1);
            D_NAc_temp = restrict(AMP_NAc_restrict, times_pre(iT)-1, times_pre(iT)+1);
            D_CG_temp = restrict(AMP_CG_restrict, times_pre(iT)-1, times_pre(iT)+1);
            
            % run it for the OFC-NAc pair
            [ac_ON, lag_ON] = xcov(D_OFC_temp.data,D_NAc_temp.data,200,  'coeff');
            max_ac_ON = ac_ON(lag_ON==0);
            amp_ac.OFC_NAc_pre_restrict(iF, iT) = max_ac_ON;
            amp_lag.OFC_NAc_pre_restrict(iF, iT) = times_pre(iT);
            
            % run it for the OFC-CG pair
            [ac_OC, lag_OC] = xcov(D_OFC_temp.data,D_CG_temp.data,200,  'coeff');
            max_ac_OC = ac_OC(lag_OC==0);
            amp_ac.OFC_CG_pre_restrict(iF, iT) = max_ac_OC;
            amp_lag.OFC_CG_pre_restrict(iF, iT) = times_pre(iT);
            
            % run it for the NAc-CG pair ( not used in
            % figure but good to see)
            [ac_NC, lag_NC] = xcov(D_NAc_temp.data,D_CG_temp.data,200,  'coeff');
            max_ac_NC = ac_NC(lag_NC==0);
            amp_ac.NAc_CG_pre_restrict(iF, iT) = max_ac_NC;
            amp_lag.NAc_CG_pre_restrict(iF, iT) = times_pre(iT);
            
        end
    end
        disp('Pre inset done')

    % Same thing for the ipsi condition
        
    % correct for the time between sessions
    Ipsi_OFC_pot = data.ipsi.OFC_pot;
    Ipsi_OFC_pot.tvec = (AMP_OFC_i.tvec-AMP_OFC_i.tvec(1))+AMP_OFC.tvec(end);
    Ipsi_NAc_pot = data.ipsi.NAc_pot;
    Ipsi_NAc_pot.tvec = (AMP_NAc_i.tvec-AMP_NAc_i.tvec(1))+AMP_NAc.tvec(end);
    Ipsi_CG_pot = data.ipsi.CG_pot;
    Ipsi_CG_pot.tvec = (AMP_CG_i.tvec-AMP_CG_i.tvec(1))+AMP_CG.tvec(end);
    
    OFC_data_cut = restrict(Ipsi_OFC_pot, cfg.ipsi_x_lim(1)-1.5, cfg.ipsi_x_lim(2)+1.5);
    NAc_data_cut = restrict(Ipsi_NAc_pot, cfg.ipsi_x_lim(1)-1.5, cfg.ipsi_x_lim(2)+1.5);
    CG_data_cut = restrict(Ipsi_CG_pot, cfg.ipsi_x_lim(1)-1.5, cfg.ipsi_x_lim(2)+1.5);
    
    times_ipsi =  cfg.ipsi_x_lim(1)-1:cfg.Fs_restrict:cfg.ipsi_x_lim(2)+1;
    % cycle across frequencies
    for iF = 1:length(Freq_list)
        
        this_F = Freq_list(iF);
        
        
        if iF>1
            while last_string_l>0
                fprintf('\b'); % delete previous counter display
                last_string_l = last_string_l - 1;
            end
        end
        last_string_l =fprintf('%0.2f', this_F);
        %%
        cfg_filter= [];
        cfg_filter.display_filter = 0;
        cfg_filter.f = [this_F-1 this_F+1];
        cfg_filter.order = 4;
        cfg_filter.verbose = 0;
        cfg_filter.type = 'fdesign';%'cheby1';
        
        FILT_OFC = FilterLFP(cfg_filter, OFC_data_cut);
        FILT_NAc = FilterLFP(cfg_filter, NAc_data_cut);
        FILT_CG = FilterLFP(cfg_filter, CG_data_cut);
        
        
        AMP_OFC_restrict = FILT_OFC;
        AMP_NAc_restrict = FILT_NAc;
        AMP_CG_restrict = FILT_CG;
        
        AMP_OFC_restrict.data = abs(hilbert(FILT_OFC.data));
        AMP_NAc_restrict.data = abs(hilbert(FILT_NAc.data));
        AMP_CG_restrict.data = abs(hilbert(FILT_CG.data));
        
        % create a Amp-xcorr-ogram
        
        for iT = 1:length(times_ipsi)
            
            D_OFC_temp = restrict(AMP_OFC_restrict, times_ipsi(iT)-1, times_ipsi(iT)+1);
            D_NAc_temp = restrict(AMP_NAc_restrict, times_ipsi(iT)-1, times_ipsi(iT)+1);
            D_CG_temp = restrict(AMP_CG_restrict, times_ipsi(iT)-1, times_ipsi(iT)+1);
            
            % run it for the OFC-NAc pair
            [ac_ON, lag_ON] = xcov(D_OFC_temp.data,D_NAc_temp.data,200,  'coeff');
            max_ac_ON = ac_ON(lag_ON==0);
            amp_ac.OFC_NAc_ipsi_restrict(iF, iT) = max_ac_ON;
            amp_lag.OFC_NAc_ipsi_restrict(iF, iT) = times_ipsi(iT);
            
            % run it for the OFC-CG pair
            [ac_OC, lag_OC] = xcov(D_OFC_temp.data,D_CG_temp.data,200,  'coeff');
            max_ac_OC = ac_OC(lag_OC==0);
            amp_ac.OFC_CG_ipsi_restrict(iF, iT) = max_ac_OC;
            amp_lag.OFC_CG_ipsi_restrict(iF, iT) = times_ipsi(iT);
            
            % run it for the NAc-CG pair ( not used in
            % figure but good to see)
            [ac_NC, lag_NC] = xcov(D_NAc_temp.data,D_CG_temp.data,200,  'coeff');
            max_ac_NC = ac_NC(lag_NC==0);
            amp_ac.NAc_CG_ipsi_restrict(iF, iT) = max_ac_NC;
            amp_lag.NAc_CG_ipsi_restrict(iF, iT) = times_ipsi(iT);
            
        end
    end
    disp('  Ipsi inset done')
    
    %% run smoothing if required.  [working on this]
     if cfg.moving_avg
        if strcmp(cfg.moving_avg, '1D')
            for ii = size(amp_ac.OFC_NAc_pre_restrict, 1):-1:1
                amp_ac.OFC_NAc_pre_restrict_out(ii,:) = conv(amp_ac.OFC_NAc_pre_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_NAc_pre_restrict_out(ii,:) = conv(amp_lag.OFC_NAc_pre_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_ac.OFC_NAc_ipsi_restrict_out(ii,:) = conv(amp_ac.OFC_NAc_ipsi_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_NAc_ipsi_restrict_out(ii,:) = conv(amp_lag.OFC_NAc_ipsi_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                
                 amp_ac.OFC_CG_pre_restrict_out(ii,:) = conv(amp_ac.OFC_CG_pre_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_CG_pre_restrict_out(ii,:) = conv(amp_lag.OFC_CG_pre_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_ac.OFC_CG_ipsi_restrict_out(ii,:) = conv(amp_ac.OFC_CG_ipsi_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_CG_ipsi_restrict_out(ii,:) = conv(amp_lag.OFC_CG_ipsi_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                
                 amp_ac.NAc_CG_pre_restrict_out(ii,:) = conv(amp_ac.NAc_CG_pre_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.NAc_CG_pre_restrict_out(ii,:) = conv(amp_lag.NAc_CG_pre_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_ac.NAc_CG_ipsi_restrict_out(ii,:) = conv(amp_ac.NAc_CG_ipsi_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.NAc_CG_ipsi_restrict_out(ii,:) = conv(amp_lag.NAc_CG_ipsi_restrict(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            end
        elseif strcmp(cfg.moving_avg, '1D_freq')
            for ii = size(amp_ac.OFC_NAc_pre_restrict, 2):-1:1
                amp_ac.OFC_NAc_pre_restrict_out(:,ii) = conv(amp_ac.OFC_NAc_pre_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_NAc_pre_restrict_out(:,ii) = conv(amp_lag.OFC_NAc_pre_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_ac.OFC_NAc_ipsi_restrict_out(:,ii) = conv(amp_ac.OFC_NAc_ipsi_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_NAc_ipsi_restrict_out(:,ii) = conv(amp_lag.OFC_NAc_ipsi_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');  
                
                amp_ac.OFC_CG_pre_restrict_out(:,ii) = conv(amp_ac.OFC_CG_pre_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_CG_pre_restrict_out(:,ii) = conv(amp_lag.OFC_CG_pre_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_ac.OFC_CG_ipsi_restrict_out(:,ii) = conv(amp_ac.OFC_CG_ipsi_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.OFC_CG_ipsi_restrict_out(:,ii) = conv(amp_lag.OFC_CG_ipsi_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same'); 
                
                amp_ac.NAc_CG_pre_restrict_out(:,ii) = conv(amp_ac.NAc_CG_pre_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.NAc_CG_pre_restrict_out(:,ii) = conv(amp_lag.NAc_CG_pre_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_ac.NAc_CG_ipsi_restrict_out(:,ii) = conv(amp_ac.NAc_CG_ipsi_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
                amp_lag.NAc_CG_ipsi_restrict_out(:,ii) = conv(amp_lag.NAc_CG_ipsi_restrict(:,ii), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same'); 
            end
            
        elseif strcmp(cfg.moving_avg, '2D')
            
            amp_ac.OFC_NAc_pre_restrict_out = conv2(amp_ac.OFC_NAc_pre_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_lag.OFC_NAc_pre_restrict_out = conv2(amp_lag.OFC_NAc_pre_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_ac.OFC_NAc_ipsi_restrict_out = conv2(amp_ac.OFC_NAc_ipsi_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_lag.OFC_NAc_ipsi_restrict_out = conv2(amp_lag.OFC_NAc_ipsi_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            
            amp_ac.OFC_CG_pre_restrict_out = conv2(amp_ac.OFC_CG_pre_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_lag.OFC_CG_pre_restrict_out = conv2(amp_lag.OFC_CG_pre_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_ac.OFC_CG_ipsi_restrict_out = conv2(amp_ac.OFC_CG_ipsi_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_lag.OFC_CG_ipsi_restrict_out = conv2(amp_lag.OFC_CG_ipsi_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            
            amp_ac.NAc_CG_pre_restrict_out = conv2(amp_ac.NAc_CG_pre_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_lag.NAc_CG_pre_restrict_out = conv2(amp_lag.NAc_CG_pre_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_ac.NAc_CG_ipsi_restrict_out = conv2(amp_ac.NAc_CG_ipsi_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
            amp_lag.NAc_CG_ipsi_restrict_out = conv2(amp_lag.NAc_CG_ipsi_restrict, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
       
        else
            amp_ac.OFC_NAc_pre_restrict_out = amp_ac.OFC_NAc_pre_restrict;
            amp_ac.OFC_NAc_ipsi_restrict_out = amp_ac.OFC_NAc_ipsi_restrict;
            
            amp_ac.OFC_CG_pre_restrict_out = amp_ac.OFC_CG_pre_restrict;
            amp_ac.OFC_CG_ipsi_restrict_out = amp_ac.OFC_CG_ipsi_restrict;
            
            amp_ac.NAc_CG_pre_restrict_out = amp_ac.NAc_CG_pre_restrict;
            amp_ac.NAc_CG_ipsi_restrict_out = amp_ac.NAc_CG_ipsi_restrict;
        end
        
    end
    
    %%
    
    figure(222)
    subplot(2,1,1)
    ax1 = imagesc('XData',times_pre,'CData', amp_ac.OFC_NAc_pre_restrict_out);
    axis xy
    set(gca, 'ytick',[1 size(amp_ac.OFC_CG_pre_restrict,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'xticklabel', [])
%     set(gca, 'ytick',[Freq_list(1) Freq_list(end)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)))%, 'xticklabel', [0:500:2500])
    xlim(cfg.pre_x_lim)
    subplot(2,1,2)
    hold on
    plot(AMP_OFC_p.tvec-AMP_OFC_p.tvec(1),AMP_OFC_p.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot(AMP_OFC_p.tvec-AMP_OFC_p.tvec(1),FILT_OFC_p.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
    %         plot(AMP_1.tvec-AMP_1.tvec(1), data.pre.OFC_pot.data)
    plot(AMP_NAc_p.tvec-AMP_NAc_p.tvec(1),AMP_NAc_p.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
    plot(AMP_NAc_p.tvec-AMP_NAc_p.tvec(1),FILT_NAc_p.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
    %         plot(AMP_2.tvec-AMP_2.tvec(1), data.pre.NAc_pot.data)
    
    xlim(cfg.pre_x_lim)
    cfg_fig = [];
    cfg_fig.ft_size =36;
    xlim(cfg.pre_x_lim)
    set(gca, 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'ytick', [])
    SetFigure(cfg_fig,gcf);
    tightfig
    
    saveas(gcf, [save_dir sess '_' S{1} '_' S{2} 'OFC_NAc_pre_sample'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
    saveas(gcf, [save_dir 'OFC_NAc_pre_sample'], 'png')
    
    fname = [save_dir 'OFC_NAc_pre_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    
    
    figure(223)
    subplot(2,1,1)
    ax1 = imagesc('XData',times_ipsi,'CData', amp_ac.OFC_NAc_ipsi_restrict_out);
    axis xy
    set(gca, 'ytick',[1 size(amp_ac.OFC_CG_ipsi_restrict,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', [cfg.ipsi_x_lim(1):0.5:cfg.ipsi_x_lim(2)], 'xticklabel', [])
    xlim(cfg.ipsi_x_lim)
    subplot(2,1,2)
    hold on
    plot((AMP_OFC_i.tvec-AMP_OFC_i.tvec(1))+AMP_OFC.tvec(end),AMP_OFC_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot((AMP_OFC_i.tvec-AMP_OFC_i.tvec(1))+AMP_OFC.tvec(end),FILT_OFC_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
    %         plot(AMP_1_i.tvec-AMP_1_i.tvec(1), data.ipsi.OFC_pot.data)
    plot((AMP_NAc_i.tvec-AMP_NAc_i.tvec(1))+AMP_NAc.tvec(end),AMP_NAc_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
    plot((AMP_NAc_i.tvec-AMP_NAc_i.tvec(1))+AMP_NAc.tvec(end),FILT_NAc_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
    %         plot(AMP_2_i.tvec-AMP_2_i.tvec(1), data.ipsi.NAc_pot.data)
    
    cfg_fig = [];
    cfg_fig.ft_size =36;
    xlim(cfg.ipsi_x_lim)
    set(gca, 'xtick', [cfg.ipsi_x_lim(1):0.5:cfg.ipsi_x_lim(2)], 'ytick', [])
    SetFigure(cfg_fig,gcf);
    tightfig
    
    saveas(gcf, [save_dir sess '_' S{1} '_' S{2} 'OFC_NAc_ipsi_sample'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
    saveas(gcf, [save_dir 'OFC_NAc_ipsi_sample'], 'png')
    
    fname = [save_dir 'OFC_NAc_ipsi_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    %%
    figure(223)
    subplot(2,1,1)
    ax1 = imagesc('XData',times_pre,'CData', amp_ac.OFC_CG_pre_restrict_out);
    axis xy
    set(gca, 'ytick',[1 size(amp_ac.OFC_CG_pre_restrict,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'xticklabel', [])
    %     set(gca, 'ytick',[Freq_list(1) Freq_list(end)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)))%, 'xticklabel', [0:500:2500])
    xlim(cfg.pre_x_lim)
    subplot(2,1,2)
    hold on
    plot(AMP_OFC_p.tvec-AMP_OFC_p.tvec(1),AMP_OFC_p.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot(AMP_OFC_p.tvec-AMP_OFC_p.tvec(1),FILT_OFC_p.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
    %         plot(AMP_1.tvec-AMP_1.tvec(1), data.pre.OFC_pot.data)
    plot(AMP_CG_p.tvec-AMP_CG_p.tvec(1),AMP_CG_p.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
    plot(AMP_CG_p.tvec-AMP_CG_p.tvec(1),FILT_CG_p.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
    %         plot(AMP_2.tvec-AMP_2.tvec(1), data.pre.NAc_pot.data)
    
    xlim(cfg.pre_x_lim)
    cfg_fig = [];
    cfg_fig.ft_size =36;
    xlim(cfg.pre_x_lim)
    set(gca, 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'ytick', [])
    SetFigure(cfg_fig,gcf);
    tightfig
    %% add a scale bar
    
        
        % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
        h1 = get(gcf);
        D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
        set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
        %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
        saveas(gcf, [save_dir 'OFC_CG_pre_sample'], 'png')
        
        fname = [save_dir  'OFC_CG_pre_sample'];
        pushdir(save_dir);
        eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
        popdir;
        close all
        
        

        figure(224)
        subplot(2,1,1)
        ax1 = imagesc('XData',times_ipsi,'CData', amp_ac.OFC_CG_ipsi_restrict_out);
        axis xy
        set(gca, 'ytick',[1 size(amp_ac.OFC_CG_ipsi_restrict,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', [cfg.ipsi_x_lim(1):0.5:cfg.ipsi_x_lim(2)], 'xticklabel', [])
        %     set(gca, 'ytick',[Freq_list(1) Freq_list(end)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)))%, 'xticklabel', [0:500:2500])
        xlim(cfg.ipsi_x_lim)
        subplot(2,1,2)
        hold on
        plot((AMP_OFC_i.tvec-AMP_OFC_i.tvec(1))+AMP_OFC.tvec(end),AMP_OFC_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
        plot((AMP_OFC_i.tvec-AMP_OFC_i.tvec(1))+AMP_OFC.tvec(end),FILT_OFC_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
        %         plot(AMP_1_i.tvec-AMP_1_i.tvec(1), data.ipsi.OFC_pot.data)
        plot((AMP_CG_i.tvec-AMP_CG_i.tvec(1))+AMP_CG.tvec(end),AMP_CG_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
        plot((AMP_CG_i.tvec-AMP_CG_i.tvec(1))+AMP_CG.tvec(end),FILT_CG_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
        xlim(cfg.ipsi_x_lim)
        
        
        cfg_fig = [];
        cfg_fig.ft_size =36;
        xlim(cfg.ipsi_x_lim)
        set(gca, 'xtick', [cfg.ipsi_x_lim(1):0.5:cfg.ipsi_x_lim(2)], 'ytick', [])
        SetFigure(cfg_fig,gcf);
        tightfig
        
        % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
        h1 = get(gcf);
        D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
        set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
        %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
        saveas(gcf, [save_dir 'OFC_CG_ipsi_sample'], 'png')
        
        fname = [save_dir 'OFC_CG_ipsi_sample'];
        pushdir(save_dir);
        eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
        popdir;
        close all
        
        
end

end



