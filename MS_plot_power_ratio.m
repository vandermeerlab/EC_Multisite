function MS_plot_power_ratio(cfg_in, Naris_in)
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
cfg_def.plot_type = 'raw';
cfg_def.ylims = [-75 75];
cfg_def.ylims_norm = [0 3];
cfg_def.pot_trk = 'pot';
cfg = ProcessConfig2(cfg_def, cfg_in);
global PARAMS
c_ord = [linspecer(length(PARAMS.Phases)); [.6 .6 .6]];

%% Collect all the power ratios across the different sites per subject/session
AOC_low = []; AOC_high = [];
AOC_con_low = []; AOC_con_high = [];
types = {'Pxx', 'White_Pxx'};
for iType = 1:length(types)
    all_AOC_low.(types{iType}) = [];
    all_AOC_high.(types{iType}) = [];
    all_AOC_con_low.(types{iType}) = [];
    all_AOC_con_high.(types{iType}) = [];
end
subjects = fieldnames(Naris_in);
for iType = 1:length(types)
    it_log = {};
    for iSub = 1:length(subjects)
        sess_list = fieldnames(Naris_in.(subjects{iSub}));
        sites = {'PL'    'IL'    'OFC'    'PiriO'    'NAc'    'Piri_N'    'CG'};
        % create the empty array
        AOC_low.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
        AOC_high.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
        AOC_con_low.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
        AOC_con_high.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
        for iSess = 1:length(sess_list);
            %
            for iSite = 1:length(sites)
                site_idx = strcmp(Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_labels, [sites{iSite} '_' cfg.pot_trk]);
                if sum(site_idx) ==1
                    site_idx = find(site_idx ==1);
                    AOC_low.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio.(types{iType}).low(:,site_idx);
                    AOC_high.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio.(types{iType}).high(:,site_idx);
                    con_types = fieldnames(Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}));
                    AOC_con_low.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}).(con_types{1})(:,site_idx);
                    AOC_con_high.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}).(con_types{1})(:,site_idx);
                end
            end
            it_log{end+1,1} =  sess_list{iSess};
        end
        all_AOC_low.(types{iType}) = cat(3,all_AOC_low.(types{iType}), AOC_low.(types{iType}));
        all_AOC_high.(types{iType}) = cat(3,all_AOC_high.(types{iType}), AOC_high.(types{iType}));
        all_AOC_con_low.(types{iType}) = cat(3,all_AOC_con_low.(types{iType}), AOC_con_low.(types{iType}));
        all_AOC_con_high.(types{iType}) = cat(3,all_AOC_con_high.(types{iType}), AOC_con_high.(types{iType}));
    end
end

%% get the stats
for iType = 1:length(types)
    AOC_low.(types{iType}) = nanmean(all_AOC_low.(types{iType}), 3)';
    AOC_high.(types{iType}) = nanmean(all_AOC_high.(types{iType}),3)';
    AOC_con_low.(types{iType}) =  nanmean(all_AOC_con_low.(types{iType}),3)';
    AOC_con_high.(types{iType}) =  nanmean(all_AOC_con_high.(types{iType}),3)';
    % try it relative ot the control condition
    
    %     rm_idx = strfind(sites, 'Piri'); Index = find(not(cellfun('isempty', rm_idx)));
    %     AOC_low.(types{iType})([4, 6],:) = [];
    %     AOC_high.(types{iType})([4, 6],:) = [];
    %     AOC_con_low.(types{iType})([4, 6],:) = [];
    %     AOC_con_high.(types{iType})([4, 6],:) = [];
    %     sites{1,6} = [];sites{1,4} = [];
    %     sites = sites(~cellfun('isempty',sites));
    for iPhase = 1:size(AOC_low.(types{iType}),2)
        norm_low.(types{iType})(:,iPhase) = AOC_low.(types{iType})(:,iPhase) ./ AOC_low.(types{iType})(:,5);
        norm_high.(types{iType})(:,iPhase) = AOC_high.(types{iType})(:,iPhase) ./ AOC_high.(types{iType})(:,5);
        norm_con_low.(types{iType})(:,iPhase) = AOC_con_low.(types{iType})(:,iPhase) ./ AOC_con_low.(types{iType})(:,5);
        norm_con_high.(types{iType})(:,iPhase) = AOC_con_high.(types{iType})(:,iPhase) ./ AOC_con_high.(types{iType})(:,5);
    end
    %% switch between normalized to control and raw
    %     switch cfg.plot_type
    %
    %         case 'raw'
    %             % plot
    %             figure(iType)
    %             subplot(4,1,1)
    %             b= bar(AOC_low.(types{iType}));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             title(types{iType})
    %             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
    %             leg_val = PARAMS.Phases; leg_val{5} = 'Control';
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims])
    %
    %             subplot(4,1,2)
    %             b = bar(AOC_high.(types{iType}));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims])
    %
    %             subplot(4,1,3)
    %             b = bar(AOC_con_low.(types{iType}));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):25:cfg.ylims(2)])
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims])
    %
    %             subplot(4,1,4)
    %             b = bar(AOC_con_high.(types{iType}));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):25:cfg.ylims(2)])
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims])
    %
    %         case 'norm'
    %             % plot
    %             figure(iType)
    %             subplot(4,1,1)
    %             b= bar(norm_low.(types{iType})(:,1:4));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             title(types{iType})
    %             set(gca, 'xticklabel', sites)
    %             leg_val = PARAMS.Phases; leg_val{5} = 'Control';
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims_norm])
    %
    %             subplot(4,1,2)
    %             b = bar(norm_high.(types{iType})(:,1:4));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             set(gca, 'xticklabel', sites)
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims_norm])
    %
    %             subplot(4,1,3)
    %             b = bar(norm_con_low.(types{iType})(:,1:4));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             set(gca, 'xticklabel', sites)
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims_norm])
    %
    %             subplot(4,1,4)
    %             b = bar(norm_con_high.(types{iType})(:,1:4));
    %             for iPhase = 1:5
    %                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
    %             end
    %             set(gca, 'xticklabel', sites)
    %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
    %             ylim([cfg.ylims_norm])
    %
    %             %% save the figure
    %             if isunix
    %                 saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type])
    %                 saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'png')
    %             else
    %                 saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type])
    %                 saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type], 'png')
    %             end
    close all
    %%%%%%%%%%%%%%%%%%% Multisite 4-site Vs cross Piriform figure %%%%%%%%%%%%%
    sites{4} = 'PC_O_F_C'; sites{6} = 'PC_N_A_c';
    cfg_plt1.pos = [600 50 560*1.4 560*1.8];
    cfg_plt1.ft_size = 18;
    
    for iFig = 1:2
        if iFig ==1
            F_id = 'Four';
            s_idx = [1,2,3,5,7]; % corresponds to the PL, OFC, NAc, and CG
        elseif iFig == 2
            F_id = 'Piri';
            s_idx = [3:6]; % corresponds to OFC, OFC_Piri, NAc, and NAc_Piri
        end
        %% Gerenate the summary figure for the PL, OFC, NAc, CG
        switch cfg.plot_type
            
            case 'raw'
                % plot
                figure(iType)
                subtightplot(9,3,[2,3,5,6],0.1) % strange subplots make for nice figures.  Ignore white space...
                b= bar(AOC_low.(types{iType})(s_idx,:));
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['low gamma AOC (' num2str(cfg.power_ratio.gamma_freq(1,1)) '-' num2str(cfg.power_ratio.gamma_freq(1,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
                ylim([cfg.ylims])
                SetFigure(cfg_plt1, gcf)
                
                subtightplot(9,3,[8,9,11,12],0.1)
                b = bar(AOC_high.(types{iType})(s_idx,:));
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['high gamma AOC (' num2str(cfg.power_ratio.gamma_freq(2,1)) '-' num2str(cfg.power_ratio.gamma_freq(2,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
                ylim([cfg.ylims])
                SetFigure(cfg_plt1, gcf)
                
                subtightplot(9,3,[17,18,20,21],0.1)
                b = bar(AOC_con_low.(types{iType})(s_idx,:));
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['low control AOC (' num2str(cfg.power_ratio.contrast(1,1)) '-' num2str(cfg.power_ratio.contrast(1,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
                %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
                ylim([cfg.ylims])
                SetFigure(cfg_plt1, gcf)
                
                subtightplot(9,3,[23,24,26,27],0.1)
                b = bar(AOC_con_high.(types{iType})(s_idx,:));
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['high control AOC (' num2str(cfg.power_ratio.contrast(2,1)) '-' num2str(cfg.power_ratio.contrast(2,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
                %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
                ylim([cfg.ylims])
                
                SetFigure(cfg_plt1, gcf)
                %%
            case 'norm'
                % plot
                figure(iType)
                subtightplot(12,1,1:3,0.1)
                b= bar(norm_low.(types{iType})(s_idx,:), 'BaseValue', 1);
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['low gamma AOC (' num2str(cfg.power_ratio.gamma_freq(1,1)) '-' num2str(cfg.power_ratio.gamma_freq(1,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx))
                ylim([cfg.ylims_norm])
                
                subtightplot(12,1,4:6,0.1)
                b = bar(norm_high.(types{iType})(s_idx,:), 'BaseValue', 1);
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['high gamma AOC (' num2str(cfg.power_ratio.gamma_freq(2,1)) '-' num2str(cfg.power_ratio.gamma_freq(2,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx))
                ylim([cfg.ylims_norm])
                
                
                subtightplot(12,1,7:9,0.1)
                b = bar(norm_con_low.(types{iType})(s_idx,:), 'BaseValue', 1);
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['low control AOC (' num2str(cfg.power_ratio.contrast(1,1)) '-' num2str(cfg.power_ratio.contrast(1,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx))
                %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
                ylim([cfg.ylims_norm])
                
                subtightplot(12,1,10:12,0.1)
                b = bar(norm_con_high.(types{iType})(s_idx,:), 'BaseValue', 1);
                for iPhase = 1:5
                    set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
                end
                title(['high control AOC (' num2str(cfg.power_ratio.contrast(2,1)) '-' num2str(cfg.power_ratio.contrast(2,2)) 'Hz)'], 'fontweight', 'normal');
                set(gca, 'xticklabel', sites(s_idx))
                %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
                ylim([cfg.ylims_norm])
                
%                 cfg_plt1.pos = [600 50 560*1.4 560*1.8];
                cfg_plt1.ft_size = 18;
                SetFigure(cfg_plt1, gcf)
        end
        %%
        mkdir(PARAMS.inter_dir, 'AOC_fit')
        if isunix
            saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type])
            saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'png')
            saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'epsc')
            
        else
            saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type])
            saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type], 'png')
            saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type], 'epsc')
        end
        close all
        
    end
end
% create the legend values and add some space
figure(9999)
leg_val = PARAMS.Phases; leg_val{5} = 'control';
b = bar(norm_high.(types{iType})(s_idx,:), 'BaseValue', 1);
for iPhase = 1:5
    set(b(iPhase), 'FaceColor', c_ord(iPhase,:), 'visible', 'off')
end
set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
l = legend(leg_val, 'location', 'south', 'orientation', 'horizontal');
legend boxoff
axis off
cfg_plt1.pos = [600 50 560*1.4 560*1.8];
cfg_plt1.ft_size = 18;
SetFigure(cfg_plt1, gcf)
if isunix
    saveas(gcf, [PARAMS.inter_dir '/AOC_fit/legend'], 'epsc')
else
    saveas(gcf, [PARAMS.inter_dir '\AOC_fit\legend'], 'epsc')
end
end
