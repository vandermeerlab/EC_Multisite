function MS_spec_fig(cfg_in, data)
%%        :
%
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

cfg_def = [];
cfg_def.win = 512;
cfg_def.noverlap = cfg_def.win/4;
cfg_def.linewidth = 3;
cfg_def.c_ord = linspecer(4);
cfg = ProcessConfig2(cfg_def, cfg_in);

global PARAMS
mkdir(PARAMS.inter_dir, 'Spec')
if isunix
    save_dir = [PARAMS.inter_dir, '/Spec/'];
else
    save_dir = [PARAMS.inter_dir, '\Spec\'];
end

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
            cfg_artif_det.threshold = std(csc_artif.data)*5;
            cfg_artif_det.minlen = 0;
            evt_artif = TSDtoIV(cfg_artif_det,csc_artif);
            
            cfg_temp = []; cfg_temp.d = [-0.25 0.25];
            evt_artif = ResizeIV(cfg_temp,evt_artif);
            artif_idx = TSD_getidx2(data.(PARAMS.Phases{ii}).(sites{iSite}),evt_artif); % if error, try TSD_getidx (slower)
            data.(PARAMS.Phases{ii}).(sites{iSite}).data(artif_idx) = 0;
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
        all_data.type = 'AMPX';
        all_data.labels = sites{iSite};
        
        
        %% make a spectrogram of the combined data
        [~,Fp,Tp,Pp] = spectrogram(data.pre.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
        [~,Fi,Ti,Pi] = spectrogram(data.ipsi.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
        [~,Fc,Tc,Pc] = spectrogram(data.contra.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
        [~,Fpo,Tpo,Ppo] = spectrogram(data.post.(sites{iSite}).data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
        
        [~,F,T,P] = spectrogram(all_data.data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);
%         t_data = data.post.(sites{iSite}).data;
%         t_data(t_data==0) =NaN;
%         [~,Fp_c,Tp_c,Pp_c] = spectrogram(t_data,rectwin(cfg.win),cfg.noverlap,1:120,all_data.hdr.Fs);

        figure(200);
        
        %         P(P==inf) = NaN;
        ax1 = imagesc(T,F,10*log10(P)); % converting to dB as usual
        set(ax1, 'AlphaData', ~isinf(10*log10(P)))
        set(gca,'FontSize',28);
        axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
        ax = gca; 
        ax.YTickLabels = {0 [], [], 60, [], [], 120}; % set
        set(gca, 'tickdir', 'out')
        
        caxis([-150 -80])
        h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );
        for jj = 1:length(h)
            h(jj).LineWidth = cfg.linewidth;
        end
        % add lines to distinguish phases
        yL = get(gca, 'ylim');
        % add lines to distinguish phases
        ylim([-4.5 yL(2)]);        
        hold on
        %         plot(0:Tp(end), zeros(length(0:Tp(end)),1),'color', cfg.c_ord(1,:), 'linewidth', 5); % pre
        %         plot(Tp(end):(Tp(end)+Ti(end)), zeros(length(Tp(end):(Tp(end)+Ti(end))),1),'color', cfg.c_ord(2,:), 'linewidth', 5); % ipsi
        %         plot((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end)), zeros(length((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))),1),'color', cfg.c_ord(3,:), 'linewidth', 5); % contra
        %         plot((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end)), zeros(length((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))),1),'color', cfg.c_ord(4,:), 'linewidth', 5); % post
        
        h(1) =rectangle('position', [0 -4.5 length(0:Tp(end)) 5.5], 'facecolor',cfg.c_ord(1,:), 'edgecolor', cfg.c_ord(1,:));
        h(2) =rectangle('position', [length(0:Tp(end)) -4.5 length(Tp(end):(Tp(end)+Ti(end))) 5.5], 'facecolor',cfg.c_ord(2,:), 'edgecolor', cfg.c_ord(2,:));
        h(3) =rectangle('position', [length(0:(Tp(end)+Ti(end))) -4.5 length(Ti(end):(Tp(end)+Ti(end)+Tc(end))) 5.5], 'facecolor',cfg.c_ord(3,:), 'edgecolor', cfg.c_ord(3,:));
        h(4) = rectangle('position',[length(0:(Tp(end)+Ti(end)+Tc(end))) -4.5 length(Tc(end):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))) 5.5], 'facecolor',cfg.c_ord(4,:), 'edgecolor', cfg.c_ord(4,:));
        
        text(length(0:Tp(end))/2, -1.5, 'pre', 'fontsize', 28, 'horizontalAlignment', 'center')
        text(length(0:Tp(end))+(length(Tp(end):(Tp(end)+Ti(end)))/2), -1.5, 'ipsi', 'fontsize', 28, 'horizontalAlignment', 'center')
        text(length(0:Tp(end)+Ti(end))+(length(Ti(end):(Ti(end)+Tc(end)))/2), -1.5, 'contra', 'fontsize', 28, 'horizontalAlignment', 'center')
        text(length(0:Tp(end)+Ti(end)+Tc(end))+(length(Tc(end):(Tc(end)+Tpo(end)))/2), -1.5, 'post', 'fontsize', 28, 'horizontalAlignment', 'center')
        
        
        
        set(h(:), 'linewidth', cfg.linewidth);
        colormap('PARULA');
        cfg_fig.ft_size = 28;
        SetFigure(cfg_fig,gcf);
        set(gcf, 'position', [0 50 1600*.9 420*.9]);
        tightfig;
        sess = strrep(data.pre.ExpKeys.date, '-', '_');
        if isfield(data.pre.ExpKeys, 'Subject')
            sess = [data.contra.ExpKeys.Subject '_' strrep(data.contra.ExpKeys.date, '-', '-')];
            subject = data.contra.ExpKeys.Subject;
        else
            sess = [data.contra.ExpKeys.ratID '_' strrep(data.contra.ExpKeys.date, '-', '-')];
            subject = data.contra.ExpKeys.ratID;
        end
        %% save the session wide spectrogram
        saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'png')
                saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'fig')
        h = get(gcf);
        D = h.PaperPosition; % Returns 1x4 vector [left right width height]
        set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
        
        saveas_eps([subject '_' sess '_' sites{iSite} '_Spectrogram'], save_dir);
        close(200)
        %% make the cross-spectral correlation matrix
        for ii  =1:4
            if ii == 1
                Px = Pp; Fx = Fp; label = 'pre';
                Color =  cfg.c_ord(1,:);
            elseif ii == 2
                Px = Pi; Fx = Fi; label = 'ipsi';
                Color =  cfg.c_ord(2,:);
            elseif ii == 3
                Px = Pc; Fx = Fc; label = 'contra';
                Color =  cfg.c_ord(3,:);
            elseif ii == 4
                Px = Ppo; Fx = Fpo; label = 'post';
                Color =  cfg.c_ord(4,:);
            end
            figure(200)
            [x_cor,p_cor] = corrcoef(Px');
            imagesc(Fx,Fx,x_cor)
            axis xy;
            SetFigure([], gcf);
            pbaspect([1 1 1])
            set(gcf, 'position', [358   167   642   631])
            set(gca, 'xtick', [], 'ytick', [], 'xColor', Color, 'yColor', Color, 'linewidth', 6);
            
            rectangle('position', [0.5, 0.5, length(Fx), length(Fx)], 'EdgeColor', Color, 'linewidth', 6)
            tightfig
            saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_xcorr' label], 'fig')
            
            saveas_eps([subject '_' sess '_' sites{iSite} '_xcorr' label], save_dir)
            %             saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_xcorr' label], 'epsc')
            close(200)
            %             %% try it with correlation
            %
            %             figure(201)
            %             imagesc(Fx,Fx,p_cor)
            %             axis xy;
            %             colormap(flipud(parula))
            %             caxis([0 0.06])
            %             SetFigure([], gcf);
            %             pbaspect([1 1 1])
            %             set(gcf, 'position', [358   167   642   631]);
            %             tightfig
            %
            % make a combined figure as a png for quick inspection.
            figure(400)
            subplot(2,4,ii)
            x_cor = corrcoef(Px');
            imagesc(Fx,Fx,x_cor)
            pbaspect([1 1 1])
            axis xy;
            xlabel(label)
            set(gca, 'xtick', [], 'ytick', []);

        end
        subplot(4,2,5:8)
        ax2 =imagesc(T,F,10*log10(P));
        title(sites{iSite})
        set(ax2, 'AlphaData', ~isinf(10*log10(P)))
        set(gca,'FontSize',20);
        axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');
        caxis([-150 -80])
        h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );
    
         yL = get(gca, 'ylim');
        % add lines to distinguish phases
        ylim([-4.5 yL(2)]);        hold on

        h(1) =rectangle('position', [0 -4.5 length(0:Tp(end)) 5.5], 'facecolor',cfg.c_ord(1,:), 'edgecolor', cfg.c_ord(1,:));
        h(2) =rectangle('position', [length(0:Tp(end)) -4.5 length(Tp(end):(Tp(end)+Ti(end))) 5.5], 'facecolor',cfg.c_ord(2,:), 'edgecolor', cfg.c_ord(2,:));
        h(3) =rectangle('position', [length(0:(Tp(end)+Ti(end))) -4.5 length(Ti(end):(Tp(end)+Ti(end)+Tc(end))) 5.5], 'facecolor',cfg.c_ord(3,:), 'edgecolor', cfg.c_ord(3,:));
        h(4) = rectangle('position',[length(0:(Tp(end)+Ti(end)+Tc(end))) -4.5 length(Tc(end):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))) 5.5], 'facecolor',cfg.c_ord(4,:), 'edgecolor', cfg.c_ord(4,:));
        
        text(length(0:Tp(end))/2, -1.5, 'pre', 'fontsize', 20, 'horizontalAlignment', 'center')
        text(length(0:Tp(end))+(length(Tp(end):(Tp(end)+Ti(end)))/2), -1.5, 'ipsi', 'fontsize', 20, 'horizontalAlignment', 'center')
        text(length(0:Tp(end)+Ti(end))+(length(Ti(end):(Ti(end)+Tc(end)))/2), -1.5, 'contra', 'fontsize', 20, 'horizontalAlignment', 'center')
        text(length(0:Tp(end)+Ti(end)+Tc(end))+(length(Tc(end):(Tc(end)+Tpo(end)))/2), -1.5, 'post', 'fontsize', 20, 'horizontalAlignment', 'center')
        
        set(h(:), 'linewidth', cfg.linewidth);
        colormap('PARULA');
        set(gcf, 'position', [327   248   995   543]);
     
        saveas_eps([subject '_' sess '_' sites{iSite} '_all'], save_dir)

        close(400)
    end
    close all
end
end