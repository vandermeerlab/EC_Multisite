function MS_plot_session_phase(cfg_in, Naris)
%% MS_plot_session_phase: ;loops over subjects and plots session wide phase
% measures (amplitude cross-corrolation "amp"; phase coherence "coh").
%
%
%
%
%
% INPUTS:
%      - cfg_in: [struct] contains all the configuration paramaters
%                    - cfg.measure: 'amp', 'coh', or 'both';
%                    determines which phase measure to plot
%
%      - Naris: [struct] contains the outputs from MS_amp_xcorr_session for
%      "amp" or MS_coh_session for "coh"


%% setup parameters
global PARAMS
cfg_def = [];
cfg_def.measure = 'coh';
cfg_def.plot_type = 'all'; % can be "all", 'no_piri'
cfg_def.color.blue = double([158,202,225])/255;
cfg_def.color.green = double([168,221,181])/255;
cfg_def.linewidth = 4;
cfg_def.ylim = [0 1];
cfg_def.legend = 'off';
cfg = ProcessConfig2(cfg_def, cfg_in);

%% master legend with all combinations
if strcmp(cfg.legend, 'on')
    figure(1)
    for iS = 1:length(PARAMS.all_pairs)
        hold on
        h(iS)=plot(1:10 , NaN(1,10), 'color', PARAMS.pair_c_ord(iS,:), 'linewidth', 8);
    end
    set(gca, 'visible', 'off')
    legend(PARAMS.all_pairs)
    legend boxoff
    cfg.set_fig.ft_size = 22;
    SetFigure(cfg.set_fig, gcf)
    saveas(gcf, [PARAMS.inter_dir '/sess/' 'legend_tall'], 'fig')
    %     saveas(gcf, [PARAMS.inter_dir '/sess/' 'legend_tall'], 'epsc')
    saveas_eps('legend_tall',[PARAMS.inter_dir 'sess/'])
    
    close(1)
    figure(2)
    subplot(3,1,1)
    for iS = 1:7
        hold on
        h(iS)=plot(1:10 , NaN(1,10), 'color', PARAMS.pair_c_ord(iS,:), 'linewidth', 8);
    end
    set(gca, 'visible', 'off')
    legend(PARAMS.all_pairs{1:7}, 'orientation', 'horizontal', 'units', 'normalized');
    legend boxoff
    
    subplot(3,1,2)
    for iS = 8:14
        hold on
        h(iS)=plot(1:10 , NaN(1,10), 'color', PARAMS.pair_c_ord(iS,:), 'linewidth', 8);
    end
    set(gca, 'visible', 'off')
    legend(PARAMS.all_pairs{8:14}, 'orientation', 'horizontal', 'units', 'normalized');
    legend boxoff
    subplot(3,1,3)
    for iS = 15:21
        hold on
        h(iS)=plot(1:10 , NaN(1,10), 'color', PARAMS.pair_c_ord(iS,:), 'linewidth', 8);
    end
    set(gca, 'visible', 'off')
    legend(PARAMS.all_pairs{15:21}, 'orientation', 'horizontal', 'units', 'normalized');
    legend boxoff
    cfg.set_fig.ft_size = 22;
    cfg.set_fig.resize = 1;
    cfg.set_fig.pos = [600 50 560*1.8 420*1];
    SetFigure(cfg.set_fig, gcf)
    saveas(gcf, [PARAMS.inter_dir '/sess/' 'legend_long'], 'fig')
    %     saveas(gcf, [PARAMS.inter_dir '/sess/' 'legend_long'], 'epsc')
    saveas_eps('legend_long',[PARAMS.inter_dir 'sess/'])
    
    close(2)
end
%% loop through session and average the values
Subjects = fieldnames(Naris);
for iSub = 1:length(Subjects)
    sess_list = fieldnames(Naris.(Subjects{iSub}));
    % loop over sessions to get the average
    
    switch cfg.measure
        %% get the coherence across session %%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'coh'
            pairs =fieldnames(Naris.(Subjects{iSub}).(sess_list{1}).coh.cxx);
            for iS = 1:length(pairs)
                for iPhase = 1:length(PARAMS.Phases)
                    all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
                    all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
                    for iSess = 1:length(sess_list)
                        all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase})');
                        all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase})');
                    end
                    all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = mean(all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}));
                    all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = mean(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}));
                end
            end
            %% actually plot everything
            figure(iSub)
            %             c_ord = linspecer(length(pairs));
            labels = {};
            rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
            rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
            for iS = 1:length(pairs)
                if strcmp(cfg.plot_type, 'no_piri')
                    if isempty(strfind(pairs{iS}, 'Piri'))
                        c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                        hold on
                        plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).contra, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).contra, 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                        xlim([0 100])
                        labels{end+1} = pairs{iS};
                    end
                else
                    hold on
                    c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                    plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).contra, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).contra, 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                    xlim([0 100])
                    labels{end+1} = pairs{iS};
                    
                end
            end
            
            %             legend(labels, 'location', 'Northwest')
            %             legend boxoff
            
            for iS = 1:length(pairs)
                if strcmp(cfg.plot_type, 'no_piri')
                    if isempty(strfind(pairs{iS}, 'Piri'))
                        c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                        hold on
                        plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).ipsi, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).ipsi,'--', 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                        xlim([0 100])
                    end
                else
                    c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                    hold on
                    plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).ipsi, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).ipsi,'--', 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                    xlim([0 100])
                end
            end
            xlabel('Frequency (Hz)');
            ylabel('Coherence')
            ylim(cfg.ylim)
            set(gca, 'ytick', [0:.5:1], 'xtick',[0 100])
            cfg.set_fig =[];
            cfg.set_fig.ft_size = 28;
            SetFigure(cfg.set_fig, gcf)
            
            if strcmp(cfg.plot_type, 'no_piri')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_coh_no_piri'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_coh_no_piri'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_coh_no_piri'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_coh_no_piri']',[PARAMS.inter_dir 'sess/'])
                
            else
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_coh'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_coh'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_coh'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_coh']',[PARAMS.inter_dir 'sess/'])
                
            end
            
            
            %% create a figure for the OFC_NAc as the example for Figure
            if ismember('OFC_NAc', pairs)
                figure(iSub+1)
                c_ord = linspecer(length(pairs));
                rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3])
                rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
                
                hold on
                plot(all_naris.(Subjects{iSub}).coh.fxx.OFC_NAc.contra, all_naris.(Subjects{iSub}).coh.cxx.OFC_NAc.contra, 'color', [0.3639    0.5755    0.7484], 'linewidth', cfg.linewidth)
                xlim([0 100])
                legend({'OFC - NAc'}, 'location', 'Northwest')
                legend boxoff
                plot(all_naris.(Subjects{iSub}).coh.fxx.OFC_NAc.ipsi, all_naris.(Subjects{iSub}).coh.cxx.OFC_NAc.ipsi, '--', 'color',[0.3639    0.5755    0.7484] , 'linewidth', cfg.linewidth)
                xlim([0 100])
                
                xlabel('Frequency (Hz)');
                ylabel('Coherence')
                ylim(cfg.ylim)
                set(gca, 'ytick', [0:0.5:1], 'xtick',[0 100])
                cfg.set_fig.ft_size = 28;
                SetFigure(cfg.set_fig, gcf)
                
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_NAc'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_NAc'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_NAc'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_OFC_NAc']',[PARAMS.inter_dir 'sess/'])
                
                
                set(findall(gca, 'Type', 'Line'),'LineWidth',6)
                ylabel([]); xlabel([]);
                set(gca, 'ytick', [0 1], 'xtick',[0 100])
                cfg.set_fig.ft_size = 36;
                SetFigure(cfg.set_fig, gcf)
                legend off
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_NAc_small'], 'epsc')
                close all
            end
            
            
            if ismember('OFC_CG', pairs)
                figure(iSub+2)
                c_ord = linspecer(length(pairs));
                rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3])
                rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
                
                hold on
                plot(all_naris.(Subjects{iSub}).coh.fxx.OFC_CG.contra, all_naris.(Subjects{iSub}).coh.cxx.OFC_CG.contra, 'color', [0.9153    0.2816    0.2878], 'linewidth', cfg.linewidth)
                xlim([0 100])
                legend({'OFC - CG'}, 'location', 'Northwest')
                legend boxoff
                plot(all_naris.(Subjects{iSub}).coh.fxx.OFC_CG.ipsi, all_naris.(Subjects{iSub}).coh.cxx.OFC_CG.ipsi, '--', 'color',[0.9153    0.2816    0.2878], 'linewidth', cfg.linewidth)
                xlim([0 100])
                
                
                xlabel('Frequency (Hz)');
                ylabel('Coherence')
                ylim(cfg.ylim)
                set(gca, 'ytick', [0:0.5:1], 'xtick',[0 100])
                cfg.set_fig.ft_size = 36;
                SetFigure(cfg.set_fig, gcf)
                
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_CA'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_CA'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_CA'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_OFC_CA']',[PARAMS.inter_dir 'sess/'])
                
                
                set(findall(gca, 'Type', 'Line'),'LineWidth',6)
                set(gca, 'ytick', [0 1], 'xtick',[0 100])
                ylabel([]); xlabel([]);
                cfg.set_fig.ft_size = 36;
                SetFigure(cfg.set_fig, gcf)
                legend off
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_CG_small'], 'epsc')
                close all
            end
            %% same for amp    %%%%%%%%%%%%%%%%%%%%%%
        case 'amp'
            pairs =fieldnames(Naris.(Subjects{iSub}).(sess_list{1}).amp.ac);
            for iS = 1:length(pairs)
                for iPhase = 1:length(PARAMS.Phases)
                    all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
                    all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
                    for iSess = 1:length(sess_list)
                        all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}));
                        all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}));
                    end
                    all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}) = nanmean(all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}));
                    all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}) = nanmean(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}));
                end
            end
            
            %% actually plot everything
            figure(iSub)
            %             c_ord = linspecer(length(pairs));
            labels = {};
            rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
            rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
            for iS = 1:length(pairs)
                if strcmp(cfg.plot_type, 'no_piri')
                    if isempty(strfind(pairs{iS}, 'Piri'))
                        c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                        hold on
                        plot(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).contra, all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).contra, 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                        xlim([0 100])
                        labels{end+1} = pairs{iS};
                    end
                else
                    hold on
                    c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                    plot(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).contra, all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).contra, 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                    xlim([0 100])
                    labels{end+1} = pairs{iS};
                    
                end
            end
            
            %
            %                legend(labels, 'location', 'Northwest')
            %             legend boxoff
            
            for iS = 1:length(pairs)
                if strcmp(cfg.plot_type, 'no_piri')
                    if isempty(strfind(pairs{iS}, 'Piri'))
                        c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                        hold on
                        plot(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).ipsi, all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).ipsi,'--', 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                        xlim([0 100])
                    end
                else
                    c_idx = find(strcmp(PARAMS.all_pairs,pairs{iS}));
                    hold on
                    plot(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).ipsi, all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).ipsi,'--', 'color', PARAMS.pair_c_ord(c_idx, :), 'linewidth', cfg.linewidth)
                    xlim([0 100])
                end
            end
            xlabel('Frequency (Hz)');
            ylabel('Amplitude xcorr')
            ylim(cfg.ylim)
            set(gca, 'ytick', [0:0.5:1], 'xtick',[0 100])

            cfg.set_fig =[];
            cfg.set_fig.ft_size = 28;
            SetFigure(cfg.set_fig, gcf)
            
            if strcmp(cfg.plot_type, 'no_piri')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_no_piri'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_no_piri'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_no_piri'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_amp_no_piri']',[PARAMS.inter_dir 'sess/'])
                
            else
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_amp']',[PARAMS.inter_dir 'sess/'])
                
            end
            %% create a figure for the OFC_NAc as the example for Figure
            if ismember('OFC_NAc', pairs)
                figure(iSub+1)
                c_ord = linspecer(length(pairs));
                rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3])
                rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
                
                hold on
                plot(all_naris.(Subjects{iSub}).amp.f.OFC_NAc.contra, all_naris.(Subjects{iSub}).amp.ac.OFC_NAc.contra, 'color', [0.3639    0.5755    0.7484], 'linewidth', cfg.linewidth)
                xlim([0 100])
                legend({'OFC - NAc'}, 'location', 'Northwest')
                legend boxoff
                plot(all_naris.(Subjects{iSub}).amp.f.OFC_NAc.ipsi, all_naris.(Subjects{iSub}).amp.ac.OFC_NAc.ipsi, '--', 'color',[0.3639    0.5755    0.7484] , 'linewidth', cfg.linewidth)
                xlim([0 100])
                
                xlabel('Frequency (Hz)');
                ylabel('Amplitude xcorr')
                ylim(cfg.ylim)
                set(gca, 'ytick', [0:0.5:1], 'xtick',[0 100])
                cfg.set_fig.ft_size = 28;
                SetFigure(cfg.set_fig, gcf)
                
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_amp_OFC_NAc']',[PARAMS.inter_dir 'sess/'])
                
                
                set(findall(gca, 'Type', 'Line'),'LineWidth',6)
                ylabel([]); xlabel([]);
                set(gca, 'ytick', [0 1], 'xtick',[0 100])
                cfg.set_fig.ft_size = 36;
                SetFigure(cfg.set_fig, gcf)
                legend off
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc_small'], 'epsc')
                close all
            end
            
            
            if ismember('OFC_CG', pairs)
                figure(iSub+2)
                c_ord = linspecer(length(pairs));
                rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3])
                rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
                
                hold on
                plot(all_naris.(Subjects{iSub}).amp.f.OFC_CG.contra, all_naris.(Subjects{iSub}).amp.ac.OFC_CG.contra, 'color', [0.9153    0.2816    0.2878], 'linewidth', cfg.linewidth)
                xlim([0 100])
                legend({'OFC - CG'}, 'location', 'Northwest')
                legend boxoff
                plot(all_naris.(Subjects{iSub}).amp.f.OFC_CG.ipsi, all_naris.(Subjects{iSub}).amp.ac.OFC_CG.ipsi, '--', 'color',[0.9153    0.2816    0.2878], 'linewidth', cfg.linewidth)
                xlim([0 100])
                
                
                xlabel('Frequency (Hz)');
                ylabel('Amplitude xcorr')
                ylim(cfg.ylim)
                set(gca, 'ytick', [0:0.5:1], 'xtick',[0 100])
                cfg.set_fig.ft_size = 36;
                SetFigure(cfg.set_fig, gcf)
                
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_CG'], 'png')
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_CG'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_CG'], 'epsc')
                saveas_eps([Subjects{iSub} '_sess_amp_OFC_CG']',[PARAMS.inter_dir 'sess/'])
                
                
                set(findall(gca, 'Type', 'Line'),'LineWidth',6)
                set(gca, 'ytick', [0 1], 'xtick',[0 100])
                ylabel([]); xlabel([]);
                cfg.set_fig.ft_size = 36;
                SetFigure(cfg.set_fig, gcf)
                legend off
                saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_OFC_CG_small'], 'epsc')
                close all
            end
            
            
            
            
            
            
            
            
            
            
            
            %             %% do both %%%%%%%%%%%%%%%%%%%%%%%%
            %             case 'both'
            %             pairs =fieldnames(Naris.(Subjects{iSub}).(sess_list{1}).coh.cxx);
            %             for iS = 1:length(pairs)
            %                 for iPhase = 1:length(PARAMS.Phases)
            %                     all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
            %                     all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
            %                     for iSess = 1:length(sess_list)
            %                         all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase})');
            %                         all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase})');
            %                     end
            %                     all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = mean(all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).(PARAMS.Phases{iPhase}));
            %                     all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}) = mean(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).(PARAMS.Phases{iPhase}));
            %                 end
            %             end
            %
            %             %% actually plot everything
            %             figure(iSub)
            %             c_ord = linspecer(length(pairs));
            %             rectangle('position', [45, 0, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
            %             rectangle('position', [70, 0, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
            %             for iS = 1:length(pairs)
            %                 hold on
            %                 plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).contra, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).contra, 'color', c_ord(iS, :), 'linewidth', cfg.linewidth)
            %                 xlim([0 100])
            %             end
            %             legend(pairs)
            %             for iS = 1:length(pairs)
            %                 hold on
            %                 plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).ipsi, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).ipsi,'--', 'color', c_ord(iS, :), 'linewidth', cfg.linewidth)
            %                 xlim([0 100])
            %             end
            %             xlabel('Frequency (Hz)');
            %             ylabel('Coherence')
            %             ylim(cfg.ylim)
            %             cfg.set_fig.ft_size = 22;
            %             SetFigure(cfg.set_fig, gcf)
            %
            %             saveas(gcf, [PARAMS.inter_dir Subject{iS} '_sess_coh'], 'png')
            %             saveas(gcf, [PARAMS.inter_dir Subject{iS} '_sess_coh'], 'fig')
            %             saveas(gcf, [PARAMS.inter_dir Subject{iS} '_sess_coh'], 'epsc')
    end
end
end
