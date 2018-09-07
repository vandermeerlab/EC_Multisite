%% temp psd processing

global PARAMS
fprintf(PARAMS.log, date);
PARAMS.inter_dir = '/Volumes/Fenrir/MS_temp/';

% Extract the data from each recroding phase within each session and separate pot vs track sections
cd(PARAMS.inter_dir)


for iSub = 4:length(PARAMS.Subjects)
    
    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'])
    
    %%%%%%%%%%%%%%%% Power Measures %%%%%%%%%%%%%%%%%
    
    
    %% generate the PSD for each session for each site
        fprintf(PARAMS.log,'\n\nExtracting Power Metrics');
        sess_list = fieldnames(data);
        for iSess  = 1:length(sess_list)
            fprintf(['Session ' sess_list{iSess} '\n'])
            fprintf(PARAMS.log,['\nGetting Power ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_collect_psd([],data.(strrep(sess_list{iSess}, '-', '_')));
            fprintf(PARAMS.log, '...complete');
        end
% end
    %% plot the PSD for each session and each site
    cfg_psd = [];
    %         cfg_psd.type = 'white';
    MS_plot_psd(cfg_psd, Naris.(PARAMS.Subjects{iSub}))
    
    
    sess_list = fieldnames(data);
    for iSess = 1:length(sess_list)
        fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        fprintf(['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        MS_spec_fig([], data.(strrep(sess_list{iSess}, '-', '_')));
        fprintf('...complete\n');
        fprintf(PARAMS.log, '...complete');
    end
    
    clear data
    %% plot the gamma band power ratios
    cfg_pow_ratio_plot = [];
% cfg_pow_ratio_plot.ylims = [-100 100];
% cfg_pow_ratio_plot.plot_type = 'raw';
% cfg_pow_ratio_plot.ylims_norm = [-2 2];
% temporary
% cfg_pow_ratio_plot.power_ratio.contrast = [25 45; 90 110];
% cfg_pow_ratio_plot.power_ratio.gamma_freq = [45 65; 70 90];

% cfg_pow_ratio_plot.pot_trk = 'pot'; 
MS_plot_power_ratio(cfg_pow_ratio_plot, Naris)
    
    
    
end

%% get the average PSD for each site

MS_plot_psd_avg([], Naris)