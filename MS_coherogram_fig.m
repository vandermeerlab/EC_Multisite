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
cfg_def.fpass=[0 120]; % frequencies of interest
cfg_def.tapers=[5 9]; % tapers
cfg_def.trialave=0; % average over trials
cfg_def.err=0; % no error computation
cfg_def.pad = -1;

movingwin=[10 0.5]; % set the moving window dimensions

cfg = ProcessConfig2(cfg_def, cfg_in);

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
            cfg_artif_det.threshold = std(csc_artif.data)*6;
            cfg_artif_det.minlen = 0;
            evt_artif = TSDtoIV(cfg_artif_det,csc_artif);
            
            cfg_temp = []; cfg_temp.d = [-0.5 0.5];
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
    
    imagesc(t,f,C'); axis xy
    
    title(['Coherogram: ' S{1} ' ' S{2}] )
    xlabel('time (s)')
    ylabel('frequency')
    %         caxis([0.2 1])
    
    set(gca,'FontSize',20);
    xlabel('time (s)'); ylabel('Frequency (Hz)');
    % add lines to distinguish phases
    ylim([0 120]);
    hold on
    plot(Tp(1):Tp(end), zeros(length(Tp(1):Tp(end)),1),'color', cfg.c_ord(1,:), 'linewidth', 5); % pre
    plot(Tp(end):(Tp(end)+Ti(end)), zeros(length(Tp(end):(Tp(end)+Ti(end))),1),'color', cfg.c_ord(2,:), 'linewidth', 5); % ipsi
    plot((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end)), zeros(length((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))),1),'color', cfg.c_ord(3,:), 'linewidth', 5); % contra
    plot((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end)), zeros(length((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))),1),'color', cfg.c_ord(4,:), 'linewidth', 5); % post
        h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );

    
    set(h(:), 'linewidth', cfg.linewidth);
    colormap('PARULA');
    
    SetFigure([],gcf);
    set(gcf, 'position', [0 50 1600*.9 420*.9]);
    tightfig;
    hc = colorbar;
    caxis([0 1])
    set(hc,'YTick',[0 1])
    
%     sess = strrep(data.pre.ExpKeys.date, '-', '_');
if isfield(data.pre.ExpKeys, 'Subject')
    sess = [data.contra.ExpKeys.Subject '_' strrep(data.contra.ExpKeys.date, '-', '-')];
    subject = data.contra.ExpKeys.Subject;
else
    sess = [data.contra.ExpKeys.ratID '_' strrep(data.contra.ExpKeys.date, '-', '-')];
    subject = data.contra.ExpKeys.ratID;
end
sess = strrep(sess, '-', '_');
%% save the session wide spectrogram
saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'png')
%saveas(gcf, [save_dir sess '_' S{1} '_' S{2} '_Coherogram'], 'fig')
h1= get(gcf);
D = h1.PaperPosition; % Returns 1x4 vector [left right width height]
set(gcf, 'PaperSize', [D(3) D(4)]); %default PaperSize is [8.5 11]
%         saveas(gcf, [save_dir subject '_' sess '_' sites{iSite} '_Spectrogram'], 'epsc')
fname = [sess '_' S{1} '_' S{2} '_Coherogram'];
pushdir(save_dir)
eval(sprintf('print -depsc2 -tiff -r300 -painters %s',fname));
popdir
close all
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
