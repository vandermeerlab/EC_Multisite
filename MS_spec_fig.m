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
cfg = ProcessConfig2(cfg_def, cfg_in);

global PARAMS
%% append the data
sites = fieldnames(data.pre);
for iSite =1:length(sites)
    all_data = [];
    if strcmp(sites{iSite}, 'ExpKeys') || strcmp(sites{iSite}, 'pos')
        continue
    else
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
        figure(100);
        imagesc(T,F,10*log10(P)); % converting to dB as usual
        set(gca,'FontSize',20);
        axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');
        caxis([-150 -80])
        h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );
        set(h(:), 'linewidth', cfg.linewidth)
        colormap('PARULA')
        
        SetFigure([],gcf)
        set(gcf, 'position', [600 50 1600 780]);
        sess = strrep(data.pre.ExpKeys.date, '-', '_');
        subject  = data.pre.ExpKeys.ratID;
        saveas(gcf, [PARAMS.inter_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'png')
        saveas(gcf, [PARAMS.inter_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'fig')
        saveas(gcf, [PARAMS.inter_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
        
        
        %% make the cross-spectral correlation matrix
        for ii  =1:4
            if ii == 1
                Px = Pp; Fx = Fp; label = 'pre';
            elseif ii == 2
                Px = Pi; Fx = Fi; label = 'ipsi';
            elseif ii == 3
                Px = Pc; Fx = Fc; label = 'contra';
            elseif ii == 4
                Px = Ppo; Fx = Fpo; label = 'post';
            end
            figure(200)
            x_cor = corrcoef(Px');
            imagesc(Fx,Fx,x_cor)
            axis xy;
            SetFigure([], gcf);
            pbaspect([1 1 1])
            set(gcf, 'position', [358   167   642   631]);
            saveas(gcf, [PARAMS.inter_dir subject '_' sess '_' sites{iSite} '_xcorr' label], 'fig')
            saveas(gcf, [PARAMS.inter_dir subject '_' sess '_' sites{iSite} '_xcorr' label], 'epsc')
            close(200)
            % make a combined figure as a png for quick inspection.
            figure(400)
            subplot(2,4,ii)
            x_cor = corrcoef(Px');
            imagesc(Fx,Fx,x_cor)
            pbaspect([1 1 1])
            axis xy;
            xlabel(label)
        end
        subplot(4,2,5:8)
        title(sites{iSite})
        imagesc(T,F,10*log10(P))
        set(gca,'FontSize',20);
        axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');
        caxis([-150 -80])
        h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );
        set(h(:), 'linewidth', cfg.linewidth)
        colormap('PARULA')
        saveas(gcf, [PARAMS.inter_dir subject '_' sess '_' sites{iSite} '_all'], 'png')
        close(400)
    end
    close all
end
end