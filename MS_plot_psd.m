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
    case {'standard', 'both'} % make use the standard PSD
        sub_list = fieldnames(Naris);
        for iSub = 1:length(fieldnames(Naris))
            sess_list = fieldnames(Naris.(sub_list{iSub}));
            for iSess = 1:length(sess_list);
                h.(['n' num2str(iSess)]) = figure((iSub)*10 + (iSess));
                site_list = fieldnames(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{1}));
                for iSite = 1:length(site_list)
                    h_site.(['n' num2str(iSite)]) = figure((iSub)*100 + (iSess)*10 +(iSite));
                    for iPhase = 1:length(PARAMS.Phases)
                        hold on
                        plot(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.F,...
                            10*log10(Naris.(sub_list{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).(site_list{iSite}).psd.Pxx),...
                            'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth)
                    end
                    xlim([0 120])
                    title(strrep([sess_list{iSess} '  ' site_list{iSite}], '_', '-'))
                    legend(PARAMS.Phases)
                    SetFigure([], h_site.(['n' num2str(iSite)]))
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir '_' sess_list{iSess} '_' site_list{iSite}], 'png')
                    saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir '_' sess_list{iSess} '_' site_list{iSite}], 'eps')
                    close(h_site.(['n' num2str(iSite)]));
                end
            end
        end
        
    case{'white', 'both'};
        
        
        
        
        
end
