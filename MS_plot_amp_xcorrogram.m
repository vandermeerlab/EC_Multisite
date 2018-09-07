function MS_plot_amp_xcorrogram(cfg_in, amp_in, data)
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

cfg_def.linewidth = 3;
cfg_def.this_pair = 'OFC_NAc';
cfg_def.that_pair = 'OFC_CG';
cfg_def.pre_x_lim = [338 339];
cfg_def.ipsi_x_lim = [854 855];



cfg = ProcessConfig2(cfg_def, cfg_in);

%% 
    %switch NaN's to zero before filtering
    data.pre.OFC_pot.data(isnan(data.pre.OFC_pot.data)) = 0;
    data.pre.NAc_pot.data(isnan(data.pre.NAc_pot.data)) = 0;
    data.pre.CG_pot.data(isnan(data.pre.CG_pot.data)) = 0;
    
    cfg_filter = [];
    cfg_filter.f = [45 65]; 
    
    FILT_1 = FilterLFP(cfg_filter, data.pre.OFC_pot);
    FILT_2 = FilterLFP(cfg_filter, data.pre.NAc_pot);
    FILT_3 = FilterLFP(cfg_filter, data.pre.CG_pot);

    AMP_1 = FILT_1;
    AMP_2 = FILT_2;
    AMP_3 = FILT_3; 
    
    AMP_1.data = abs(hilbert(FILT_1.data));
    AMP_2.data = abs(hilbert(FILT_2.data));
    AMP_3.data = abs(hilbert(FILT_3.data));

    
    data.ipsi.OFC_pot.data(isnan(data.ipsi.OFC_pot.data)) = 0;
    data.ipsi.NAc_pot.data(isnan(data.ipsi.NAc_pot.data)) = 0;
    data.ipsi.CG_pot.data(isnan(data.ipsi.CG_pot.data)) = 0;
    
    FILT_1_i = FilterLFP(cfg_filter, data.ipsi.OFC_pot);
    FILT_2_i = FilterLFP(cfg_filter, data.ipsi.NAc_pot);
    FILT_3_i = FilterLFP(cfg_filter, data.ipsi.CG_pot);

    AMP_1_i = FILT_1_i;
    AMP_2_i = FILT_2_i;
    AMP_3_i = FILT_3_i; 
    
    AMP_1_i.data = abs(hilbert(FILT_1_i.data));
    AMP_2_i.data = abs(hilbert(FILT_2_i.data));
    AMP_3_i.data = abs(hilbert(FILT_3_i.data));
    
%% append data
amp_ac = amp_in.amp_ac;
amp_lag = amp_in.amp_lag; 
Freq_list = amp_in.freq;
Fs  = 0.1; 
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
    ax1 = imagesc( all_sess);
    axis xy
    set(gca, 'ytick',[1 size(all_sess,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)), 'xticklabel', [0:500:2500])
    
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
    
%% try the amp max ac lag values.  
figure(112)
    ax1 = imagesc(all_sess_lag);
    axis xy
    set(gca, 'ytick',[1 size(all_sess,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)), 'xticklabel', [0:500:2500])
    
%     set(ax1, 'AlphaData', ~isnan(all_sess))
    %     title(['Coherogram: ' S{1} ' ' S{2}] )
    xlabel('time (s)')
    ylabel('frequency')
    %         caxis([0.2 1])
    
    set(gca,'FontSize',20);
    xlabel('time (s)'); ylabel('Frequency (Hz)');
    yL = get(gca, 'ylim');
    % add lines to distinguish phases
    ylim([-4.5 yL(2)]);
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
    
    
    h = vline([length(amp_ac.(this_pair).pre), (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)),...
        (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra))], {'k','k', 'k', 'k'}) ;
    set(h(:), 'linewidth', 3)
    
    cfg_fig = [];
    cfg_fig.ft_size =30;
    hc = colorbar;
    caxis([-0.015 0.015])
    SetFigure(cfg_fig,gcf);
    
%         set(hc,'YTick',[-1 1], 'yticklabel', [-1 1])
    set(gcf, 'position', [0 50 1800*.9 480*.9]);
    
    %% save

        mkdir(PARAMS.inter_dir, 'Amp_xcorrogram')
    if isunix
        save_dir = [PARAMS.inter_dir, 'Amp_xcorrogram/'];
    else
        save_dir = [PARAMS.inter_dir, 'Amp_xcorrogram\'];
    end
    if iRun == 1
        fname = [cfg.this_pair '_Amp_lag_ogram'];
    else
        
        fname = [cfg.that_pair '_Amp_lag_ogram'];
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
%% zoom in on the some specifed sections

if iRun == 1
    
    figure(222)
    subplot(2,1,1)
    ax1 = imagesc(all_sess_og);
    axis xy
    set(gca, 'ytick',[1 size(all_sess_og,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)))%, 'xticklabel', [0:500:2500])
    xlim(cfg.pre_x_lim*1/Fs)
    subplot(2,1,2)
    hold on
    plot(AMP_1.tvec-AMP_1.tvec(1),AMP_1.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot(AMP_1.tvec-AMP_1.tvec(1),FILT_1.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
%         plot(AMP_1.tvec-AMP_1.tvec(1), data.pre.OFC_pot.data)
    plot(AMP_2.tvec-AMP_2.tvec(1),AMP_2.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
    plot(AMP_2.tvec-AMP_2.tvec(1),FILT_2.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
%         plot(AMP_2.tvec-AMP_2.tvec(1), data.pre.NAc_pot.data)
    
    xlim(cfg.pre_x_lim)
    cfg_fig = [];
    cfg_fig.ft_size =36;
    xlim(cfg.pre_x_lim)
    set(gca, 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'ytick', [])
    SetFigure(cfg_fig,gcf);
    tightfig
    
    % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
        saveas(gcf, [save_dir this_pair '_pre_sample'], 'png')

    fname = [save_dir this_pair '_pre_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    
    
    figure(223)
    subplot(2,1,1)
    ax1 = imagesc(all_sess_og);
    axis xy
    set(gca, 'ytick',[1 size(all_sess_og,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)))%, 'xticklabel', [0:500:2500])
    xlim(cfg.ipsi_x_lim*1/Fs)
    subplot(2,1,2)
    hold on
    plot((AMP_1_i.tvec-AMP_1_i.tvec(1))+AMP_1.tvec(end),AMP_1_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot((AMP_1_i.tvec-AMP_1_i.tvec(1))+AMP_1.tvec(end),FILT_1_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
%         plot(AMP_1_i.tvec-AMP_1_i.tvec(1), data.ipsi.OFC_pot.data)
    plot((AMP_2_i.tvec-AMP_2_i.tvec(1))+AMP_2.tvec(end),AMP_2_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
    plot((AMP_2_i.tvec-AMP_2_i.tvec(1))+AMP_2.tvec(end),FILT_2_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
%         plot(AMP_2_i.tvec-AMP_2_i.tvec(1), data.ipsi.NAc_pot.data)
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
        saveas(gcf, [save_dir this_pair '_ipsi_sample'], 'png')

    fname = [save_dir this_pair '_ipsi_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    
elseif iRun ==2
     figure(232)
    subplot(2,1,1)
    ax1 = imagesc(all_sess_og);
    axis xy
    set(gca, 'ytick',[1 size(all_sess_og,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)))%, 'xticklabel', [0:500:2500])
    xlim(cfg.pre_x_lim*1/Fs)
    subplot(2,1,2)
    hold on
    plot(AMP_1.tvec-AMP_1.tvec(1),AMP_1.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot(AMP_1.tvec-AMP_1.tvec(1),FILT_1.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
%         plot(AMP_1.tvec-AMP_1.tvec(1), data.pre.OFC_pot.data)
    plot(AMP_3.tvec-AMP_3.tvec(1),AMP_3.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
    plot(AMP_3.tvec-AMP_3.tvec(1),FILT_3.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
%         plot(AMP_2.tvec-AMP_2.tvec(1), data.pre.NAc_pot.data)
    
    xlim(cfg.pre_x_lim)
    cfg_fig = [];
    cfg_fig.ft_size =36;
    xlim(cfg.pre_x_lim)
    set(gca, 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'ytick', [])
    SetFigure(cfg_fig,gcf);
    tightfig
    
    % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
        saveas(gcf, [save_dir this_pair '_pre_sample'], 'png')

    fname = [save_dir this_pair '_pre_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    
    
    figure(233)
    subplot(2,1,1)
    ax1 = imagesc(all_sess_og);
    axis xy
    set(gca, 'ytick',[1 size(all_sess_og,1)], 'yticklabel',[Freq_list(1) Freq_list(end)], 'xtick', (1:5000:length(all_time)))%, 'xticklabel', [0:500:2500])
    xlim(cfg.ipsi_x_lim*1/Fs)
    subplot(2,1,2)
    hold on
    plot((AMP_1_i.tvec-AMP_1_i.tvec(1))+AMP_1.tvec(end),AMP_1_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot((AMP_1_i.tvec-AMP_1_i.tvec(1))+AMP_1.tvec(end),FILT_1_i.data+0.00025, 'color', [0.4 0.4 0.4],'Linewidth', 1)
%         plot(AMP_1_i.tvec-AMP_1_i.tvec(1), data.ipsi.OFC_pot.data)
    plot((AMP_3_i.tvec-AMP_3_i.tvec(1))+AMP_3.tvec(end),AMP_3_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 2)
    plot((AMP_3_i.tvec-AMP_3_i.tvec(1))+AMP_3.tvec(end),FILT_3_i.data, 'color', [0.6 0.6 0.6],'Linewidth', 1)
%         plot(AMP_2_i.tvec-AMP_2_i.tvec(1), data.ipsi.NAc_pot.data)
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
        saveas(gcf, [save_dir this_pair '_ipsi_sample'], 'png')

    fname = [save_dir this_pair '_ipsi_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    
    
end

end


