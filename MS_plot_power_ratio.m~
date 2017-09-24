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
subjects = fieldnames(Naris_in);
for iSub = 1:length(subjects)
    sess_list = fieldnames(Naris_in.(subjects{iSub}));
    for iSess = 1:length(sess_list); 
        for iSite = 1:length(PARAMS.all_sites)
            for iType = 1:length(types)
                if strcmp(
                AOC_low(iSite,:) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio.(types{iType})
            
            
        
        
        end
    end
end


