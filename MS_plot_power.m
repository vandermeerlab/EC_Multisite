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
cfg_def.linewidth = 2;
cfg = ProcessConfig2(cfg_def, cfg_in);
global PARAMS
c_ord = linspecer(length(PARAMS.Phases));
%% cycle through data from each subject, session, and channel to give an
%  overview figure (raster) and individual PSDs (vector)
switch cfg.type
    case {'standard'} % make use the standard PSD
        sub_list = fieldnames(Naris);
        for iSub = 1:length(fieldnames(Naris))
            sess_list = fieldnames(Naris.(sub_list{iSub}));
            for iSess = 1:length(sess_list);
                site_list = fieldnames(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{1}));
                S1 = ceil(length(site_list)/2);
                h.(['n' num2str(iSess)]) = figure((iSub)*10 + (iSess));
                iSite_pot= 0; iSite_trk = 0;
                for iSite = 1:length(site_list)
                    h_site.(['n' num2str(iSite)]) = figure((iSub)*100 + (iSess)*10 +(iSite));
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel('Frequency (Hz)')
                    ylabel('Power')
                    title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                    %legend(PARAMS.Phases)
                    SetFigure([], h_site.(['n' num2str(iSite)]))
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite}], 'png')
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite}], 'eps')
                    close(h_site.(['n' num2str(iSite)]));
                    if strcmp(site_list{iSite}(end-2:end), 'trk')
                        h_all_trk= figure((iSub)*1000 + (iSess));
                        iSite_trk = iSite_trk+1;
                        subplot(3,3,iSite_trk)
                    else
                        h_all_pot = figure((iSub)*100 + (iSess));
                        iSite_pot = iSite_pot+1;
                        subplot(3,3,iSite_pot)
                    end
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel(site_list{iSite})
                    ylabel('Power')
                    %                     legend(PARAMS.Phases)
                    
                end
                fig_pos = get(h_all_pot, 'position');
                tightfig(h_all_pot)
                set(h_all_pot, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                fig_pos = get(h_all_trk, 'position');
                tightfig(h_all_trk)
                set(h_all_trk, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                % SetFigure([], h_all)
                if isunix
                    sum_dir = '/Summary/';
                    mkdir(PARAMS.inter_dir, sum_dir)
                else
                    sum_dir = '\Summary\';
                    mkdir(PARAMS.inter_dir, sum_dir)
                end
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_pot'], 'png')
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_pot'], 'eps')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_trk'], 'png')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_trk'], 'eps')
                close all
            end
        end
    case{'white'};
        sub_list = fieldnames(Naris);
        for iSub = 1:length(fieldnames(Naris))
            sess_list = fieldnames(Naris.(sub_list{iSub}));
            for iSess = 1:length(sess_list);
                site_list = fieldnames(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{1}));
                S1 = ceil(length(site_list)/2);
                h.(['n' num2str(iSess)]) = figure((iSub)*10 + (iSess));
                iSite_pot= 0; iSite_trk = 0;
                for iSite = 1:length(site_list)
                    h_site.(['n' num2str(iSite)]) = figure((iSub)*100 + (iSess)*10 +(iSite));
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel('Frequency (Hz)')
                    ylabel('Power')
                    title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                    %legend(PARAMS.Phases)
                    SetFigure([], h_site.(['n' num2str(iSite)]))
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite} '_white'], 'png')
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite} '_white'], 'eps')
                    close(h_site.(['n' num2str(iSite)]));
                    if strcmp(site_list{iSite}(end-2:end), 'trk')
                        h_all_trk= figure((iSub)*1000 + (iSess));
                        iSite_trk = iSite_trk+1;
                        subplot(3,3,iSite_trk)
                    else
                        h_all_pot = figure((iSub)*100 + (iSess));
                        iSite_pot = iSite_pot+1;
                        subplot(3,3,iSite_pot)
                    end
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel(site_list{iSite})
                    ylabel('Power')
                    %                     legend(PARAMS.Phases, 'location', 'SouthEast')
                    
                end
                fig_pos = get(h_all_pot, 'position');
                tightfig(h_all_pot)
                set(h_all_pot, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                fig_pos = get(h_all_trk, 'position');
                tightfig(h_all_trk)
                set(h_all_trk, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                % SetFigure([], h_all)
                if isunix
                    sum_dir = '/Summary/';
                    mkdir(PARAMS.inter_dir, sum_dir)
                else
                    sum_dir = '\Summary\';
                    mkdir(PARAMS.inter_dir, sum_dir)
                end
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_pot'], 'png')
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_pot'], 'eps')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_trk'], 'png')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_trk'], 'eps')
                close all
            end
        end
        
    otherwise % b/c Matlab can't deal with "or" in a case switch...
        sub_list = fieldnames(Naris);
        for iSub = 1:length(fieldnames(Naris))
            sess_list = fieldnames(Naris.(sub_list{iSub}));
            for iSess = 1:length(sess_list);
                site_list = fieldnames(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{1}));
                S1 = ceil(length(site_list)/2);
                h.(['n' num2str(iSess)]) = figure((iSub)*10 + (iSess));
                iSite_pot= 0; iSite_trk = 0;
                for iSite = 1:length(site_list)
                    h_site.(['n' num2str(iSite)]) = figure((iSub)*100 + (iSess)*10 +(iSite));
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel('Frequency (Hz)')
                    ylabel('Power')
                    title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                    %legend(PARAMS.Phases)
                    SetFigure([], h_site.(['n' num2str(iSite)]))
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite}], 'png')
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite}], 'eps')
                    close(h_site.(['n' num2str(iSite)]));
                    if strcmp(site_list{iSite}(end-2:end), 'trk')
                        h_all_trk= figure((iSub)*1000 + (iSess));
                        iSite_trk = iSite_trk+1;
                        subplot(3,3,iSite_trk)
                    else
                        h_all_pot = figure((iSub)*100 + (iSess));
                        iSite_pot = iSite_pot+1;
                        subplot(3,3,iSite_pot)
                    end
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel(site_list{iSite})
                    ylabel('Power')
                    %                     legend(PARAMS.Phases)
                    
                end
                fig_pos = get(h_all_pot, 'position');
                tightfig(h_all_pot)
                set(h_all_pot, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                fig_pos = get(h_all_trk, 'position');
                tightfig(h_all_trk)
                set(h_all_trk, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                % SetFigure([], h_all)
                if isunix
                    sum_dir = '/Summary/';
                    mkdir(PARAMS.inter_dir, sum_dir)
                else
                    sum_dir = '\Summary\';
                    mkdir(PARAMS.inter_dir, sum_dir)
                end
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_pot'], 'png')
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_pot'], 'eps')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_trk'], 'png')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_trk'], 'eps')
                close all
                
            end
        end
        sub_list = fieldnames(Naris);
        for iSub = 1:length(fieldnames(Naris))
            sess_list = fieldnames(Naris.(sub_list{iSub}));
            for iSess = 1:length(sess_list);
                site_list = fieldnames(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{1}));
                S1 = ceil(length(site_list)/2);
                h.(['n' num2str(iSess)]) = figure((iSub)*10 + (iSess));
                iSite_pot= 0; iSite_trk = 0;
                for iSite = 1:length(site_list)
                    h_site.(['n' num2str(iSite)]) = figure((iSub)*100 + (iSess)*10 +(iSite));
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel('Frequency (Hz)')
                    ylabel('Power')
                    title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                    %legend(PARAMS.Phases)
                    SetFigure([], h_site.(['n' num2str(iSite)]))
                    Square_subplots
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite} '_white'], 'png')
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir sess_list{iSess} '_' site_list{iSite} '_white'], 'eps')
                    close(h_site.(['n' num2str(iSite)]));
                    if strcmp(site_list{iSite}(end-2:end), 'trk')
                        h_all_trk= figure((iSub)*1000 + (iSess));
                        iSite_trk = iSite_trk+1;
                        subplot(3,3,iSite_trk)
                    else
                        h_all_pot = figure((iSub)*100 + (iSess));
                        iSite_pot = iSite_pot+1;
                        subplot(3,3,iSite_pot)
                    end
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.White_Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    xlabel(site_list{iSite})
                    ylabel('Power')
                    %legend(PARAMS.Phases, 'location', 'SouthEast')
                    
                end
                fig_pos = get(h_all_pot, 'position');
                tightfig(h_all_pot)
                set(h_all_pot, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                fig_pos = get(h_all_trk, 'position');
                tightfig(h_all_trk)
                set(h_all_trk, 'position', [fig_pos(1), fig_pos(2), 700, 700])
                % SetFigure([], h_all)
                if isunix
                    sum_dir = '/Summary/';
                    mkdir(PARAMS.inter_dir, sum_dir)
                else
                    sum_dir = '\Summary\';
                    mkdir(PARAMS.inter_dir, sum_dir)
                end
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_pot'], 'png')
                saveas(h_all_pot, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_pot'], 'eps')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_trk'], 'png')
                saveas(h_all_trk, [PARAMS.inter_dir sum_dir sess_list{iSess} '_all_white_trk'], 'eps')
                close all
            end
        end
        
end
