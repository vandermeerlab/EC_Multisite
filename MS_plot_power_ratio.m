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
cfg = ProcessConfig2(cfg_def, cfg_in);
global PARAMS
c_ord = linspecer(length(PARAMS.Phases));


%% Collect all the power ratios across the different sites per subject/session
AOC_low = []; AOC_high = [];
AOC_con_low = []; AOC_con_high = [];
types = {'Pxx', 'White_Pxx'};
for iType = 1:length(types)
    all_AOC_low.(types{iType}) = [];
    all_AOC_high.(types{iType}) = [];
    all_AOC_con_low.(types{iType}) = [];
    all_AOC_con_h.(types{iType}) = [];
end
subjects = fieldnames(Naris_in);
for iType = 1:length(types)
    for iSub = 1:length(subjects)
        sess_list = fieldnames(Naris_in.(subjects{iSub}));
        for iSess = 1:length(sess_list);
            sites = PARAMS.all_sites;
            for iSite = 1:length(sites)
                site_idx = strcmp(Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_labels, [sites{iSite} '_pot']);
                if sum(site_idx) ==1
                    site_idx = find(site_idx ==1);
                    AOC_low.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio.(types{iType}).low(:,site_idx);
                    AOC_high.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio.(types{iType}).high(:,site_idx);
                    con_types = fieldnames(Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}));
                    AOC_con_low.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}).(con_types{1})(:,site_idx);
                    AOC_con_high.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}).(con_types{1})(:,site_idx);
                end
            end
        end
        all_AOC_low.(types{iType}) = cat(3,all_AOC_low.(types{iType}), AOC_low.(types{iType}));
        all_AOC_high.(types{iType}) = cat(3,all_AOC_high.(types{iType}), AOC_high.(types{iType}));
        all_AOC_con_low.(types{iType}) = cat(3,all_AOC_con_low.(types{iType}), AOC_con_low.(types{iType}));
        all_AOC_con_high.(types{iType}) = cat(3,all_AOC_con_high.(types{iType}), AOC_con_high.(types{iType}));
    end
end

%% get the stats
for iType = 1:length(types)
    AOC_low.(types{iType}) = mean(all_AOC_low.(types{iType}), 3)';
    AOC_high.(types{iType}) = mean(all_AOC_high.(types{iType}),3)';
    AOC_con_low.(types{iType}) =  mean(all_AOC_con_low.(types{iType}),3)';
    AOC_con_high.(types{iType}) =  mean(all_AOC_con_high.(types{iType}),3)';
    
    figure(iType)
    subplot(4,1,1)
    bar(AOC_low.(types{iType}))
    title(types{iType})
    set(gca, 'xticklabel', sites)
    subplot(4,1,2)
    bar(AOC_high.(types{iType}))
    set(gca, 'xticklabel', sites)
    
    subplot(4,1,3)
    bar(AOC_con_low.(types{iType}))
    set(gca, 'xticklabel', sites)
    
    subplot(4,1,4)
    bar(AOC_con_high.(types{iType}))
    set(gca, 'xticklabel', sites)
    
end
