function MS_coherogram_fig(cfg_in, data)
%% MS_coherogram_fig: conputes a coherogram for two inputs across time.
% This makes use of the Chronux Toolbox (v2.11, a the time of writing this
% function).
%
%
%          Inputs:
%           -
%           -
%           -
%          Outputs:
%           -
%           -
%           -
%
% EC - 2017-01-08

%% set up variables
cfg_def = [];
cfg_def.linewidth = 3;
cfg_def.c_ord = linspecer(4);
cfg_def.Fs=2000; % sampling frequency
cfg_def.fpass=[0 100]; % frequencies of interest
cfg_def.tapers=[5 9]; % tapers
cfg_def.trialave=0; % average over trials
cfg_def.err=0; % no error computation
cfg_def.pad = -1;
cfg_def.pre_x_lim = [338 339];
cfg_def.ipsi_x_lim = [854 855];

movingwin=[1 .05]; % set the moving window dimensions

cfg = ProcessConfig2(cfg_def, cfg_in);
c_ord = linspecer(4);
global PARAMS
mkdir(PARAMS.inter_dir, 'Coherogram')
if isunix
    save_dir = [PARAMS.inter_dir, 'Coherogram/'];
else
    save_dir = [PARAMS.inter_dir, 'Coherogram\'];
end

addpath(genpath(PARAMS.Chronux_dir))

%% append the data
sites = fieldnames(data.pre);
for iSite =1:length(sites)
    all_data = [];
    if strcmp(sites{iSite}, 'ExpKeys') || strcmp(sites{iSite}, 'pos')
        continue
    else
        
        %% remove the high amplitude and chewing artifacts
        for ii = 1:4
            csc_artif = data.(PARAMS.Phases{ii}).(sites{iSite});
            csc_artif.data = abs(csc_artif.data); % detect artifacts both ways
            
            cfg_artif_det = [];
            cfg_artif_det.method = 'raw';
            cfg_artif_det.threshold = std(csc_artif.data)*7;
            cfg_artif_det.minlen = 0;
            evt_artif = TSDtoIV(cfg_artif_det,csc_artif);
            
            cfg_temp = []; cfg_temp.d = [-0.05 0.05];
            evt_artif = ResizeIV(cfg_temp,evt_artif);
            artif_idx = TSD_getidx2(data.(PARAMS.Phases{ii}).(sites{iSite}),evt_artif); % if error, try TSD_getidx (slower)
            data.(PARAMS.Phases{ii}).(sites{iSite}).data(artif_idx) = NaN;
        end
        %%
        all_data.hdr.Fs = data.pre.(sites{iSite}).cfg.hdr{1}.SamplingFrequency;
        [x,y] = size(data.ipsi.(sites{iSite}).data);
        if x>y
            all_data.data = [data.pre.(sites{iSite}).data ; data.ipsi.(sites{iSite}).data ; data.contra.(sites{iSite}).data; data.post.(sites{iSite}).data];
        else
            all_data.data = [data.pre.(sites{iSite}).data , data.ipsi.(sites{iSite}).data , data.contra.(sites{iSite}).data, data.post.(sites{iSite}).data];
        end
        all_data.tvec = [data.pre.(sites{iSite}).tvec', data.ipsi.(sites{iSite}).tvec', data.contra.(sites{iSite}).tvec',data.post.(sites{iSite}).tvec'];
        pre_norm = (data.pre.(sites{iSite}).tvec - data.pre.(sites{iSite}).tvec(1))';
        ipsi_norm = ((data.ipsi.(sites{iSite}).tvec - data.ipsi.(sites{iSite}).tvec(1))+(pre_norm(end)+(pre_norm(end)-pre_norm(end-1))))';
        contra_norm = ((data.contra.(sites{iSite}).tvec - data.contra.(sites{iSite}).tvec(1))+(ipsi_norm(end)+(ipsi_norm(end)-ipsi_norm(end-1))))';
        post_norm = ((data.post.(sites{iSite}).tvec - data.post.(sites{iSite}).tvec(1))+(contra_norm(end)+(contra_norm(end)-contra_norm(end-1))))';
        
        all_data.tvec_norm = [pre_norm, ipsi_norm, contra_norm, post_norm];
        all_data.type = 'AMPX';
        all_data.labels = sites{iSite};
        
        Sites.(sites{iSite}) = all_data;
    end
end
clear csc_artif

%%  loop through the pairs of good cites to make the coherograms

for iPair  =1:length(data.pre.ExpKeys.GoodPairs)
    S =strsplit(data.pre.ExpKeys.GoodPairs{iPair}, '_');
    
    if strcmp(S{1}, 'PiriO')
        S{1} = 'Piri_O';
    end
    if strcmp(S{2}, 'PiriO')
        S{2} = 'Piri_O';
    end
    if strcmp(S{1}, 'PiriN')
        S{1} = 'Piri_N';
    end
    if strcmp(S{2}, 'PiriN')
        S{2} = 'Piri_N';
    end
    
    S{1} = [S{1} '_pot'];
    S{2} = [S{2} '_pot'];
    
    D1 = Sites.(S{1});
    D2 = Sites.(S{2});
    % get the time range for each naris phase for plotting the color bars.
    % there is probably a more computationally frendy way to do this...
    [~,~,~,~,~,Tp,~]=cohgramc(data.pre.(S{1}).data',data.pre.(S{2}).data',movingwin, cfg);
    [~,~,~,~,~,Ti,~]=cohgramc(data.ipsi.(S{1}).data',data.ipsi.(S{2}).data',movingwin, cfg);
    [~,~,~,~,~,Tc,~]=cohgramc(data.contra.(S{1}).data',data.contra.(S{2}).data',movingwin, cfg);
    [~,~,~,~,~,Tpo,~]=cohgramc(data.post.(S{1}).data',data.post.(S{2}).data',movingwin, cfg);
    
    %    CC = [Cp; Ci; Cc; Cpo];
    %    ff = [fp fi fc fpo];
    %    tt = [Tp Ti Tc Tpo];
    [C,~,~,~,~,t,f]=cohgramc(D1.data',D2.data',movingwin, cfg);
    
    %%
    figure(1)
    ax1 =imagesc(t,f,C'); axis xy
    
    set(ax1, 'AlphaData', ~isnan(C'))
    %     title(['Coherogram: ' S{1} ' ' S{2}] )
    xlabel('time (s)')
    ylabel('frequency')
    %         caxis([0.2 1])
    
    set(gca,'FontSize',20);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    % add lines to distinguish phases
    ylim([-6 100]);
    hold on
    rectangle('position', [Tp(1) -6 (Tp(end)-Tp(1))  5.5], 'facecolor',cfg.c_ord(1,:), 'edgecolor', cfg.c_ord(1,:))
    rectangle('position', [Tp(end) -6 ((Tp(end)+Ti(end))-Tp(end))  5.5], 'facecolor',cfg.c_ord(2,:), 'edgecolor', cfg.c_ord(2,:))
    rectangle('position', [(Tp(end)+Ti(end)) -6 ((Tp(end)+Ti(end)+Tc(end))-(Tp(end)+Ti(end)))  5.5], 'facecolor',cfg.c_ord(3,:), 'edgecolor', cfg.c_ord(3,:))
    rectangle('position', [(Tp(end)+Ti(end)+Tc(end)) -6 ((Tp(end)+Ti(end)+Tc(end)+Tpo(end))-(Tp(end)+Ti(end)+Tc(end)))  5.5], 'facecolor',cfg.c_ord(4,:), 'edgecolor', cfg.c_ord(4,:))

    text((Tp(end)-Tp(1))/2, -3.5, 'pre', 'fontsize', 20, 'horizontalAlignment', 'center')
    text(median(Tp(end):(Tp(end)+Ti(end))), -3.5, 'ipsi', 'fontsize', 20, 'horizontalAlignment', 'center')
    text(median((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))), -3.5, 'contra', 'fontsize', 20, 'horizontalAlignment', 'center')
    text(median((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))), -3.5, 'post', 'fontsize', 20, 'horizontalAlignment', 'center')

    
%     plot(Tp(1):Tp(end), zeros(length(Tp(1):Tp(end)),1)-off_set,'color', cfg.c_ord(1,:), 'linewidth', 1); % pre
% %     plot(Tp(end):(Tp(end)+Ti(end)), zeros(length(Tp(end):(Tp(end)+Ti(end))),1)-off_set,'color', cfg.c_ord(2,:), 'linewidth', 10); % ipsi
%     plot((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end)), zeros(length((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))),1)-off_set,'color', cfg.c_ord(3,:), 'linewidth', 10); % contra
%     plot((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end)), zeros(length((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))),1)-off_set,'color', cfg.c_ord(4,:), 'linewidth', 10); % post
    h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );
    %     if cfg.lfp_on
    %        plot(D1.tvec_norm,(D1.data*10000)+20, 'color', [1 1 1])
    %        plot(D2.tvec_norm, (D2.data*10000)+20, 'color', [0.5 1 1])
    %     end

    set(h(:), 'linewidth', cfg.linewidth);
    colormap('PARULA');
    %     drawArrow = @(x,y,L) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, 'LineWidth', L);  % make an arrow.  Matlab builtin function is just terrible.
    %     drawArrow([250 254], [0 10], 5)
    text(median(cfg.pre_x_lim), 102, 'V', 'fontsize', 28, 'color', c_ord(1,:))
    text(median(cfg.ipsi_x_lim), 102, 'V', 'fontsize', 28, 'Color', c_ord(2,:))
%     tightfig
    cfg_fig = [];
    cfg_fig.ft_size =30;
    SetFigure(cfg_fig,gcf);
        hc = colorbar;
%     caxis(hc, [0 1])
%     set(hc,'Ticks',[0 1], 'Tickslabels', [0 1])
    set(hc, 'Limits', [0 1], 'Ticks', [0 1])
    set(gcf, 'position', [0 50 1800*.9 480*.9]);
%     tightfig;
%     ax =gca;
%     ax.Clipping = 'off';
    %     sess = strrep(data.pre.ExpKeys.date, '-', '_');
    if isfield(data.pre.ExpKeys, 'Subject')
        sess = [data.contra.ExpKeys.Subject '_' strrep(data.contra.ExpKeys.date, '-', '-')];
        subject = data.contra.ExpKeys.Subject;
    else
        sess = [data.contra.ExpKeys.ratID '_' strrep(data.contra.ExpKeys.date, '-', '-')];
        subject = data.contra.ExpKeys.ratID;
    end
    sess = strrep(sess, '-', '_');
    
    saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'png')
    saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
    fname = [sess '_' S{1} '_' S{2} '_Coherogram'];
    pushdir(save_dir)
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir
    close all
    
    
    
    %% get a zoom in
    movingwin = [0.5 0.05];
    [C,~,~,~,~,t,f]=cohgramc(D1.data(1:floor(.5*length(D1.data)))',D2.data(1:floor(.5*length(D2.data)))',movingwin, cfg);
    
    
    figure(111)
    subplot(2,1,1)
    ax1 =imagesc(t,f,C'); axis xy
    set(ax1, 'AlphaData', ~isnan(C'))
    xlim(cfg.pre_x_lim)
    set(gca, 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'xticklabel', [])
    
    
    subplot(2,1,2)
    hold on
    plot(D1.tvec_norm(1:floor(.5*length(D1.tvec_norm))),D1.data(1:floor(.5*length(D1.data)))+0.00050, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot(D2.tvec_norm(1:floor(.5*length(D2.tvec_norm))), D2.data(1:floor(.5*length(D2.data))), 'color', [0.6 0.6 0.6], 'Linewidth', 2)
%     [axL, icons] = legend(S, 'box', 'off', 'fontsize', 18, 'location');
%     axL.FontSize = 30;
%     set(icons, 'LineWidth', 2)
    cfg_fig = [];
    cfg_fig.ft_size =36;
    xlim(cfg.pre_x_lim)
    set(gca, 'xtick', [cfg.pre_x_lim(1):0.5:cfg.pre_x_lim(2)], 'ytick', [])
    SetFigure(cfg_fig,gcf);
    tightfig
    saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_pre_sample'], 'png')
    % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
    fname = [sess '_' S{1} '_' S{2} '_pre_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    
    figure(112)
    subplot(2,1,1)
    ax1 =imagesc(t,f,C'); axis xy
    set(ax1, 'AlphaData', ~isnan(C'))
    xlim(cfg.ipsi_x_lim)
    set(gca, 'xtick', [cfg.ipsi_x_lim(1):0.5:cfg.ipsi_x_lim(2)], 'xticklabel', [])
    
    subplot(2,1,2)
    hold on
    plot(D1.tvec_norm,(D1.data)+0.00050, 'color', [0.4 0.4 0.4],'Linewidth', 2)
    plot(D2.tvec_norm, D2.data, 'color', [0.6 0.6 0.6], 'Linewidth', 2)
    xlim(cfg.ipsi_x_lim)
    set(gca, 'xtick', [cfg.ipsi_x_lim(1):0.5:cfg.ipsi_x_lim(2)], 'ytick', [])
%     [axL, icons] = legend(S, 'box', 'off', 'fontsize', 12);
%     axL.FontSize = 30;
%     set(icons, 'LineWidth', 2)
    line([cfg.ipsi_x_lim(1)-0.01 cfg.ipsi_x_lim(1)-0.01],[0 0.0002],'color',  'k', 'linewidth', 0.5)
    box off
    cfg_fig = [];
    cfg_fig.ft_size =36;
    SetFigure(cfg_fig,gcf);
    tightfig
    
    saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_ipsi_sample'], 'png');
    % saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
    h1 = get(gcf);
    D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
    set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
    %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
    fname = [sess '_' S{1} '_' S{2} '_ipsi_sample'];
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
    popdir;
    close all
    
    % figure(111)
    % subplot(1,2,2)
    % ax1 =imagesc(t,f,C'); axis xy
    % xlim([254 255])
    %
    % subplot(3,2,4:5)
    % hold on
    % plot(D1.tvec_norm,(D1.data)+0.00020, 'color', [0.4 0.4 0.4])
    % plot(D2.tvec_norm, D2.data, 'color', [0.6 0.6 0.6])
    % xlim([254 255])
    clearvars('C', 'D1', 'D2');
end

%%


% xlim([0 100])



% %% make a spectrogram of the combined data
% [~,Fp,Tp,Pp] = spectrogram(data.pre.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
% [~,Fi,Ti,Pi] = spectrogram(data.ipsi.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
% [~,Fc,Tc,Pc] = spectrogram(data.contra.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
% [~,Fpo,Tpo,Ppo] = spectrogram(data.post.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
%
%
%
% [~,F,T,P] = spectrogram(all_data.data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
% figure(100);
% P(P==inf) = NaN;
% imagesc(T,F,10*log10(P)); % converting to dB as usual
%
% caxis([-150 -80])
% h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );
% % add lines to distinguish phases
% ylim([0 120]);
% hold on
% plot(0:Tp(end), zeros(length(0:Tp(end)),1),'color', cfg.c_ord(1,:), 'linewidth', 5); % pre
% plot(Tp(end):(Tp(end)+Ti(end)), zeros(length(Tp(end):(Tp(end)+Ti(end))),1),'color', cfg.c_ord(2,:), 'linewidth', 5); % ipsi
% plot((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end)), zeros(length((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))),1),'color', cfg.c_ord(3,:), 'linewidth', 5); % contra
% plot((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end)), zeros(length((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))),1),'color', cfg.c_ord(4,:), 'linewidth', 5); % post
%
% set(h(:), 'linewidth', cfg.linewidth);
% colormap('PARULA');
%
% SetFigure([],gcf);
% set(gcf, 'position', [0 50 1600*.9 420*.9]);
% tightfig;
% sess = strrep(data.pre.ExpKeys.date, '-', '_');
% if isfield(data.pre.ExpKeys, 'Subject')
%     sess = [data.contra.ExpKeys.Subject '_' strrep(data.contra.ExpKeys.date, '-', '-')];
%     subject = data.contra.ExpKeys.Subject;
% else
%     sess = [data.contra.ExpKeys.ratID '_' strrep(data.contra.ExpKeys.date, '-', '-')];
%     subject = data.contra.ExpKeys.ratID;
% end
% %% save the session wide spectrogram
% saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Spectrogram'], 'png')
% saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Spectrogram'], 'fig')
% h = get(gcf);
% D = h.PaperPosition; % Returns 1x4 vector [left right width height]
% set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
% %         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
% fname = [subject '_' sess '_' sites{iSite} '_Spectrogram'];
% pushdir(save_dir)
% eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
% popdir
% %         print -depsc2 -tiff -r300 -painters test.eps
%
%
% end
% close all
% end
% end
