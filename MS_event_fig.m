function MS_event_fig(cfg_in, events, data)
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
cfg_def.c_ord = linspecer(7);
cfg_def.raw_plot_boost = 0.5*10^-3; % how to offset the LFP traces
cfg_def.width = 0.5; % width of sample event window (in seconds)
cfg_def.line_end  =.5;
cfg_def.width = .5;
cfg_def.ft_size = 20;
cfg_def.linewidth = 2;
cfg_def.mrk_off = -5;
cfg_def.markers = {'#', '+', 'x', 'o'};
cfg_def.type = 'low';
cfg_def.nEvt = 6;
cfg = ProcessConfig2(cfg_def, cfg_in);

global PARAMS

cfg.raw_plot_boost = (round(length(fieldnames(data.post))/1.5) *.1) * 10^-3;
%% cycle trough the sites
sites = fieldnames(events);
for iSite =1:length(sites)
    all_data = [];
    if strcmp(sites{iSite}(end-2:end), 'trk') % only use the "pot" phase
        continue
    else
        if length(events.(sites{iSite}).post.(cfg.type).tstart) < 6 % ensure there are at least some events
            nEvt = length(events.(sites{iSite}).post.(cfg.type).tstart);
        else
            nEvt = cfg.nEvt;
        end
        for iEvt = randperm(length(events.(sites{iSite}).post.(cfg.type).tstart(1:end)),nEvt) % loop through 4 random examples in the post session. 
            if iEvt == 1 || iEvt == length(events.(sites{iSite}).post.(cfg.type).tstart); % first and last events can cause errors if they fo beyond the plot range.  just skip them
                continue
            else
                ctr = mean(cat(2,events.(sites{iSite}).post.(cfg.type).tstart(iEvt),events.(sites{iSite}).post.(cfg.type).tend(iEvt)),2);
                bg_tstart_idx = nearest_idx3(ctr-cfg.width/2,data.post.(sites{iSite}).tvec);
                bg_tend_idx = nearest_idx3(ctr+cfg.width/2,data.post.(sites{iSite}).tvec);
                
                % prepare the figure paramters
                %%
                h1 = figure(101);
                set(gcf,'PaperPositionMode','auto')
                
                % plot the raw traces for the four corner sites
                set(gca, 'ColorOrder', cfg.c_ord)
                hold all
                plot_loop = 1;
                labels = {};
                y_val = [];
                for iSite2 =length(sites):-1:1
                    if strcmp(sites{iSite2}(end-2:end), 'trk')
                        continue
                    else
                        h.LFP = plot(data.post.(sites{iSite2}).tvec(bg_tstart_idx:bg_tend_idx),data.post.(sites{iSite2}).data(bg_tstart_idx:bg_tend_idx)+plot_loop*cfg.raw_plot_boost, 'linewidth', cfg.linewidth, 'color', cfg.c_ord(plot_loop, :));
                        y_val = [y_val, median(data.post.(sites{iSite2}).data(bg_tstart_idx:bg_tend_idx))+plot_loop*cfg.raw_plot_boost];
                        plot_loop = plot_loop + 1;
                        labels = [labels, sites{iSite2}(1:end-4)];
                        xlim([data.post.(sites{iSite2}).tvec(bg_tstart_idx), data.post.(sites{iSite2}).tvec(bg_tend_idx)])
                        y_lim = ylim;
                        set(gca, 'ytick', sort(y_val), 'yticklabel', labels,'TickDir','out')
                        set(gca, 'xtick', data.post.(sites{iSite2}).tvec(bg_tstart_idx):.5:data.post.(sites{iSite2}).tvec(bg_tend_idx),'xTickLabel', 0:500:500, 'fontsize', cfg.ft_size, 'fontname', 'helvetica','fontweight','bold')
                        text(data.post.(sites{iSite2}).tvec(bg_tstart_idx+round(.4*length(data.post.(sites{iSite2}).tvec(bg_tstart_idx:bg_tend_idx)))), y_lim(1)-0.5*10^4, 'Time (ms)','fontsize',cfg.ft_size, 'fontweight','bold', 'fontname', 'helvetica')
                    end
                end
                rectangle('position',[events.(sites{iSite}).post.(cfg.type).tstart(iEvt), y_lim(1)+0.1*10^-4, events.(sites{iSite}).post.(cfg.type).tend(iEvt) - events.(sites{iSite}).post.(cfg.type).tstart(iEvt), 0.5*10^-4], 'facecolor', [.6 .6 .6], 'edgecolor', [.6 .6 .6])
                %% save the figure
                if isfield(data.pre.ExpKeys, 'Subject')
                    sess = [data.contra.ExpKeys.Subject '_' strrep(data.contra.ExpKeys.date, '-', '-')];
                    subject = data.contra.ExpKeys.Subject;
                else
                    sess = [data.contra.ExpKeys.ratID '_' strrep(data.contra.ExpKeys.date, '-', '-')];
                    subject = data.contra.ExpKeys.ratID;
                end
                title([sess '-' num2str(iEvt) '-' sites{iSite}]);
                sess = strrep(sess, '-', '-');
                cfg_fig = [];
                cfg_fig.ft_size  = 24;
                SetFigure(cfg_fig, gcf)
                mkdir(PARAMS.inter_dir, 'Events');
                if isunix
                    saveas(gcf, [PARAMS.inter_dir 'Events/' subject '_' sess '_' sites{iSite} '_' num2str(iEvt) '_' cfg.type], 'png')
                    saveas(gcf, [PARAMS.inter_dir 'Events/' subject '_' sess '_' sites{iSite} '_' num2str(iEvt) '_' cfg.type ], 'fig')
                            saveas_eps([subject '_' sess '_' sites{iSite} '_' num2str(iEvt) '_' cfg.type], [PARAMS.inter_dir 'Events/'])

                else
                    saveas(gcf, [PARAMS.inter_dir 'Events\' subject '_' sess '_' sites{iSite} '_' num2str(iEvt) '_' cfg.type ], 'png')
                    saveas(gcf, [PARAMS.inter_dir 'Events\' subject '_' sess '_' sites{iSite} '_' num2str(iEvt)  '_' cfg.type], 'fig')
                            saveas_eps([subject '_' sess '_' sites{iSite} '_' num2str(iEvt) '_' cfg.type], [PARAMS.inter_dir 'Events\'])

                    
                end
            end
        end
    end
    close all
end
end