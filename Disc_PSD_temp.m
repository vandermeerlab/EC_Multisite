%% temp psd processing

global PARAMS
fprintf(PARAMS.log, date);
% Extract the data from each recroding phase within each session and separate pot vs track sections
cd(PARAMS.inter_dir)


for iSub = 1%:length(PARAMS.Subjects)
    
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
    
    
    
    
end