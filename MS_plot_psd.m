function MS_plot_psd(cfg_in, Naris)
%% MS_plot_psd: plots multiple power spectral densities for the data files
% in the "Naris" structure (output from MS_collect_psd)
%
% inputs:
%    -cfg_in: [struct] contains configuration paramters
%    -Naris: [struct] contains power and frequency values for each channel
%    for each subject/session/phase
%
%    this script currently uses a global parameter set to determine where
%    to save the output figures

%% set up defaults

cfg_def = [];
cfg_def.type = 'both'; % whether to output the 'standard' or "white" filtered PSD
cfg_def.linewidth = 4;
cfg_def.color.blue = double([158,202,225])/255;
cfg_def.color.green = double([168,221,181])/255;
cfg = ProcessConfig2(cfg_def, cfg_in);
global PARAMS
c_ord = linspecer(length(PARAMS.Phases));
mkdir(PARAMS.inter_dir, 'PSD')

%% cycle through data from each subject, session, and channel to give an
%  overview figure (raster) and individual PSDs (vector)
switch cfg.type
    case {'standard'} % make use the standard PSD
        %         sub_list = fieldnames(Naris);
        %                 sub_list = fieldnames(Naris);
        %                 for iSub = 1:length(fieldnames(Naris))
        sess_list = fieldnames(Naris);
        for iSess = 1:length(sess_list);
            site_list = fieldnames(Naris.(sess_list{iSess}).(PARAMS.Phases{1}));
            S1 = ceil(length(site_list)/2);
            h.(['n' num2str(iSess)]) = figure((iSess));
            for iSite = 1:length(site_list)
                h_site.(['n' num2str(iSite)]) = figure((iSess)*10 +(iSite));
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth);
                end
                
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-115 -85])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                clear H
                %% move site text
                text(10, -110, site_list{iSite}(1:end-4), 'fontsize', 48)
                %                     legend(PARAMS.Phases)
                cfg_f.ft_size = 24;
                SetFigure(cfg_f, h_site.(['n' num2str(iSite)]))
                
                if ~isunix
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_' site_list{iSite} ], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} ],[PARAMS.inter_dir 'PSD\'])
                    
                else
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_' site_list{iSite} ], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} ],[PARAMS.inter_dir 'PSD/'])
                    
                end
                close(h_site.(['n' num2str(iSite)]));
                
                % hold all plots for the summary
                h_all = figure(iSess);
                subplot(2,S1,iSite)
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth/2)
                end
                
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-115 -85])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                %% move site text
                text(10, -110, site_list{iSite}(1:end-4), 'fontsize', 12)
                %                     legend(PARAMS.Phases, 'location', 'SouthEast')
                
            end
            
            %Square_subplots();
            cfg_f.ft_size = 24;
            SetFigure(cfg_f, h_all)
            if ~isunix
                saveas(h_all, [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_all'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all'],[PARAMS.inter_dir 'PSD\'])
            else
                saveas(h_all, [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_all'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all'],[PARAMS.inter_dir 'PSD/'])
            end
            close all
        end
    case{'white'};
        %         sub_list = fieldnames(Naris);
        %         for iSub = 1:length(fieldnames(Naris))
        sess_list = fieldnames(Naris);
        for iSess = 1:length(sess_list);
            site_list = fieldnames(Naris.(sess_list{iSess}).(PARAMS.Phases{1}));
            S1 = ceil(length(site_list)/2);
            h.(['n' num2str(iSess)]) = figure((iSess));
            for iSite = 1:length(site_list)
                h_site.(['n' num2str(iSite)]) = figure((iSess)*10 +(iSite));
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth);
                end
                
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-135 -110])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                clear H
                %% move site text
                text(10, -115, site_list{iSite}(1:end-4), 'fontsize', 48)
                %                     legend(PARAMS.Phases)
                cfg_f.ft_size = 24;
                SetFigure(cfg_f, h_site.(['n' num2str(iSite)]))
                
                if ~isunix
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_' site_list{iSite} '_white'], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_white'],[PARAMS.inter_dir 'PSD\'])
                    
                else
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_' site_list{iSite} '_white'], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_white'],[PARAMS.inter_dir 'PSD/'])
                    
                end
                close(h_site.(['n' num2str(iSite)]));
                
                % hold all plots for the summary
                h_all = figure(iSess);
                subplot(2,S1,iSite)
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth/2)
                end
                
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-135 -110])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                %% move site text
                text(10, -115, site_list{iSite}(1:end-4), 'fontsize', 12)
                %                     legend(PARAMS.Phases, 'location', 'SouthEast')
                
            end
            
            %Square_subplots();
            cfg_f.ft_size = 24;
            SetFigure(cfg_f, h_all)
            if ~isunix
                saveas(h_all, [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_all_white'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all_white'],[PARAMS.inter_dir 'PSD\'])
            else
                saveas(h_all, [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_all_white'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all_white'],[PARAMS.inter_dir 'PSD/'])
            end
            close all
        end
        
        
    otherwise % b/c Matlab can't deal with "or" in a case switch...
        sess_list = fieldnames(Naris);
        for iSess = 1:length(sess_list);
            site_list = fieldnames(Naris.(sess_list{iSess}).(PARAMS.Phases{1}));
            S1 = ceil(length(site_list)/2);
            h.(['n' num2str(iSess)]) = figure((iSess));
            for iSite = 1:length(site_list)
                h_site.(['n' num2str(iSite)]) = figure((iSess)*10 +(iSite));
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth);
                end
                
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-115 -85])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                clear H
                %% move site text
                text(10, -110, site_list{iSite}(1:end-4), 'fontsize', 48)
                %                     legend(PARAMS.Phases)
                cfg_f.ft_size = 24;
                SetFigure(cfg_f, h_site.(['n' num2str(iSite)]))
                
                if ~isunix
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_' site_list{iSite} ], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} ],[PARAMS.inter_dir 'PSD\'])
                    
                else
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_' site_list{iSite} ], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} ],[PARAMS.inter_dir 'PSD/'])
                    
                end
                close(h_site.(['n' num2str(iSite)]));
                
                % hold all plots for the summary
                h_all = figure(iSess);
                subplot(2,S1,iSite)
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth/2)
                end
                
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-115 -85])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                %% move site text
                text(10, -110, site_list{iSite}(1:end-4), 'fontsize', 12)
                %                     legend(PARAMS.Phases, 'location', 'SouthEast')
                
            end
            
            %Square_subplots();
            cfg_f.ft_size = 24;
            SetFigure(cfg_f, h_all)
            if ~isunix
                saveas(h_all, [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_all'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all'],[PARAMS.inter_dir 'PSD\'])
            else
                saveas(h_all, [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_all'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all'],[PARAMS.inter_dir 'PSD/'])
            end
            close all
        end
        %         sub_list = fieldnames(Naris);
        %         for iSub = 1:length(fieldnames(Naris))
        sess_list = fieldnames(Naris);
        for iSess = 1:length(sess_list);
            site_list = fieldnames(Naris.(sess_list{iSess}).(PARAMS.Phases{1}));
            S1 = ceil(length(site_list)/2);
            h.(['n' num2str(iSess)]) = figure((iSess));
            for iSite = 1:length(site_list)
                h_site.(['n' num2str(iSite)]) = figure((iSess)*10 +(iSite));
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth);
                end
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-135 -110])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                clear H
                %% move site text
                text(10, -115, site_list{iSite}(1:end-4), 'fontsize', 48)
                %                     legend(PARAMS.Phases)
                cfg_f.ft_size = 24;
                SetFigure(cfg_f, h_site.(['n' num2str(iSite)]))
                
                if ~isunix
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_' site_list{iSite} '_white'], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_white'],[PARAMS.inter_dir 'PSD\'])
                    
                else
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_' site_list{iSite} '_white'], 'png')
                    saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_white'],[PARAMS.inter_dir 'PSD/'])
                    
                end
                close(h_site.(['n' num2str(iSite)]));
                
                % hold all plots for the summary
                h_all = figure(iSess);
                subplot(2,S1,iSite)
                for iPhase = 1:length(PARAMS.Phases)
                    hold on
                    plot(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                        10*log10(Naris.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                        'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth/2)
                end
                xlim([0 100])
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-135 -110])
                %                     xlabel('Frequency (Hz)')
                %                     ylabel('Power')
                %                     title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                y_val = ylim;
                rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
                rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
                H = get(gca, 'Children');
                set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
                %% move site text
                text(10, -115, site_list{iSite}(1:end-4), 'fontsize', 12)
                %                     legend(PARAMS.Phases, 'location', 'SouthEast')
            end
            
            %Square_subplots();
            cfg_f.ft_size = 24;
            SetFigure(cfg_f, h_all)
            if ~isunix
                saveas(h_all, [PARAMS.inter_dir 'PSD\' sess_list{iSess} '_all_white'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all_white'],[PARAMS.inter_dir 'PSD\'])
            else
                saveas(h_all, [PARAMS.inter_dir 'PSD/' sess_list{iSess} '_all_white'], 'png')
                saveas_eps([ sess_list{iSess} '_' site_list{iSite} '_all_white'],[PARAMS.inter_dir 'PSD/'])
            end
            close all
        end        
end
