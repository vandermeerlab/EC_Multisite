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
cfg_def.color.blue = double([158,202,225])/255;
cfg_def.color.green = double([168,221,181])/255;
cfg_def.linewidth = 4;
cfg_def.ylim = [0 1];
cfg = ProcessConfig2(cfg_def, cfg_in);


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
            c_ord = linspecer(length(pairs));
            rectangle('position', [45, 0, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
            rectangle('position', [70, 0, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
            if cfg.pairs
                for iS = 1:length(pairs)
                    if ismember(pairs{iS}, cfg.pairs)
                        hold on
                        plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).contra, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).contra, 'color', c_ord(iS, :), 'linewidth', cfg.linewidth)
                        xlim([0 100])
                    end
                end
            else
                for iS = 1:length(pairs)
                    hold on
                    plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).contra, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).contra, 'color', c_ord(iS, :), 'linewidth', cfg.linewidth)
                    xlim([0 100])
                end
                
            end
                
                legend(pairs, 'location', 'Northwest')
                for iS = 1:length(pairs)
                    hold on
                    plot(all_naris.(Subjects{iSub}).coh.fxx.(pairs{iS}).ipsi, all_naris.(Subjects{iSub}).coh.cxx.(pairs{iS}).ipsi,'--', 'color', c_ord(iS, :), 'linewidth', cfg.linewidth)
                    xlim([0 100])
                end
                xlabel('Frequency (Hz)');
                ylabel('Coherence')
                ylim(cfg.ylim)
                cfg.set_fig.ft_size = 22;
                SetFigure(cfg.set_fig, gcf)
                
                saveas(gcf, [PARAMS.inter_dir Subjects{iSub} '_sess_coh'], 'png')
                saveas(gcf, [PARAMS.inter_dir Subjects{iSub} '_sess_coh'], 'fig')
                saveas(gcf, [PARAMS.inter_dir Subjects{iSub} '_sess_coh'], 'epsc')
                
                
                
                %% same for amp    %%%%%%%%%%%%%%%%%%%%%%
                case 'amp'
                    pairs =fieldnames(Naris.(Subjects{iSub}).(sess_list{1}).amp.ac);
                    for iS = 1:length(pairs)
                        for iPhase = 1:length(PARAMS.Phases)
                            all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
                            all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}) = [];
                            for iSess = 1:length(sess_list)
                                all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase})');
                                all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}) = cat(1,all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}), Naris.(Subjects{iSub}).(sess_list{iSess}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase})');
                            end
                            all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}) = nanmean(all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).(PARAMS.Phases{iPhase}));
                            all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}) = nanmean(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).(PARAMS.Phases{iPhase}));
                        end
                    end
                    
                    %% actually plot everything
                    figure(iSub)
                    c_ord = linspecer(length(pairs));
                    rectangle('position', [45, 0, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
                    rectangle('position', [70, 0, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
                    for iS = 1:length(pairs)
                        hold on
                        plot(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).contra, all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).contra, 'color', c_ord(iS, :), 'linewidth', cfg.linewidth)
                        xlim([0 100])
                    end
                    legend(pairs)
                    for iS = 1:length(pairs)
                        hold on
                        plot(all_naris.(Subjects{iSub}).amp.f.(pairs{iS}).ipsi, all_naris.(Subjects{iSub}).amp.ac.(pairs{iS}).ipsi,'--', 'color', c_ord(iS, :), 'linewidth', cfg.linewidth)
                        xlim([0 100])
                    end
                    xlabel('Frequency (Hz)');
                    ylabel('Amplitude xcorr')
                    ylim(cfg.ylim)
                    cfg.set_fig.ft_size = 22;
                    SetFigure(cfg.set_fig, gcf)
                    
                    
                    saveas(gcf, [PARAMS.inter_dir Subjects{iSub} '_sess_amp'], 'png')
                    saveas(gcf, [PARAMS.inter_dir Subjects{iSub} '_sess_amp'], 'fig')
                    saveas(gcf, [PARAMS.inter_dir Subjects{iSub} '_sess_amp'], 'epsc')
                    
                    
                    
                    
                    
                    
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
