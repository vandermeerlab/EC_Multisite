function MS_single_analysis(analyses)
%% single analyses lets you run individiual analyses

%% make a log
% Extract the data from each recroding phase within each session and separate pot vs track sections
%% load the intermediate files
global PARAMS
fprintf(PARAMS.log, date);
fprintf(PARAMS.log,'\n\nLoading intermediates');
load([PARAMS.inter_dir 'MS_data_R104.mat'])

% load([PARAMS.inter_dir 'MS_data.mat'])
% load([PARAMS.inter_dir 'MS_naris.mat'])
load([PARAMS.inter_dir 'MS_events.mat'])

%% get the ratio of the power in multiple bands relative to the exponential f curve
if ismember('power_ratio', analyses)
    
    fprintf(PARAMS.log,'\n\nExtracting Power Metrics');
    for iSub = 1:length(PARAMS.Subjects)
        sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
        for iSess  = 1:length(sess_list)
            fprintf(['Session ' sess_list{iSess} '\n'])
            fprintf(PARAMS.log,['\nGetting Power ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_collect_psd([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
            fprintf(PARAMS.log, '...complete');
        end
    end
    
    save([PARAMS.inter_dir 'MS_narispre.mat'], 'Naris', '-v7.3')
    
    fprintf(PARAMS.log,'\n\nExtracting Power Ratio');
    for iSub = 1:length(PARAMS.Subjects)
        sess_list = fieldnames(Naris.(PARAMS.Subjects{iSub}));
        for iSess  = 1:length(sess_list)
            fprintf(PARAMS.log,['\nGetting ratio ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            cfg_pow_ratio.id = sess_list{iSess};
            [Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), cfg_p_ratio] = MS_get_power_ratio(cfg_pow_ratio,Naris.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
            %         Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')) = MS_get_power_ratio(cfg_pow_ratio,Naris_trk.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
            fprintf(PARAMS.log, '...complete');
        end
    end
    save([PARAMS.inter_dir 'MS_naris.mat'], 'Naris', '-v7.3')
    
    % plot the gamma band power ratios
    cfg_pow_ratio_plot.ylims = [-100 100];
    cfg_pow_ratio_plot.plot_type = 'raw';
    cfg_pow_ratio_plot.ylims_norm = [-2 2];
    % temporary
    cfg_pow_ratio_plot.power_ratio.contrast = [25 45; 90 110];
    cfg_pow_ratio_plot.power_ratio.gamma_freq = [45 65; 70 90];
    
    cfg_pow_ratio_plot.pot_trk = 'pot';
    MS_plot_power_ratio(cfg_pow_ratio_plot, Naris)
    
    % plot the gamma band power ratios
    cfg_pow_ratio_plot = [];
    cfg_pow_ratio_plot.plot_type = 'norm';
    cfg_pow_ratio_plot.ylims_norm = [-2 3];
    % temporary
    cfg_pow_ratio_plot.power_ratio.contrast = [25 45; 90 110];
    cfg_pow_ratio_plot.power_ratio.gamma_freq = [45 65; 70 90];
    
    cfg_pow_ratio_plot.pot_trk = 'pot';
    MS_plot_power_ratio(cfg_pow_ratio_plot, Naris)
end



%% split pot vs trk
% fprintf(PARAMS.log,'\n\nSplitting the data into pot and trk');
% [Naris_pot, Naris_trk]  = MS_pot_trk_split(Naris);
% [data_pot, data_trk]  = MS_pot_trk_split(data);
% % [Events_pot, Events_trk]  = MS_pot_trk_split(Events);


%% plot the PSDs
if ismember('plot_psd', analyses)
    cfg_psd.type = 'white';
    MS_plot_psd(cfg_psd, Naris);
    
    %% count the events
    cfg_evt_plot =[];
    cfg_evt_plot.sites = {'PL_pot', 'OFC_pot', 'NAc_pot', 'CG_pot'};
    
    MS_plot_event_stats(cfg_evt_plot, Events)
end

%% get an example event from each session and plot all sites together for the same event.
if ismember('event_fig', analyses)
    for iSub = 1:length(PARAMS.Subjects)
        sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
        for iSess = 1:length(sess_list)
            fprintf(PARAMS.log,['\nPlotting Events ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            MS_event_fig([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
            fprintf(PARAMS.log, '...complete');
        end
    end
    
end
%% generate a spectrogram across each session for each site.
if ismember('spectrogram', analyses) || ismember('spec', analyses)
    for iSub = 1:length(PARAMS.Subjects)
        sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
        for iSess = 1:length(sess_list)
            fprintf(PARAMS.log,['\nPlotting Spec ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            MS_spec_fig([], data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
            fprintf(PARAMS.log, '...complete');
        end
    end
end
%% get the phase slope values across all subjects, sessions, pairs, events
if ismember('phase', analyses)
    
    %for iSub = 1:length(PARAMS.Subjects)
iSub = 2;
        sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
%       for iSess = 1:length(sess_list)
 iSess =3; 
%if iSess ~= 2
    %fprintf(PARAMS.log,['\nExtracting phase slope ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            mat_all{iSub,iSess} = MS_get_phase_metrics_serial([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
            
%             fprintf(PARAMS.log, '...complete');

 %       end
end
   % end
    
    
    save([PARAMS.inter_dir 'MS_mat23.mat'], 'mat_all', '-v7.3')
end


%% get the phase slope values across all subjects, sessions, pairs, events
% if ismember('coh', analyses)
    
    %for iSub = 1:length(PARAMS.Subjects)
iSub = 2;
        sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
    %    for iSess = 1:length(sess_list)
 iSess =2; 
    %fprintf(PARAMS.log,['\nExtracting phase slope ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
            Coh_matl{iSub,iSess} = MS_coherence([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
            
%             fprintf(PARAMS.log, '...complete');
       end
   end
    
    
%     save([PARAMS.inter_dir 'MS_matS2.mat'], 'mat_all', '-v7.3')
% end
%%
% MS_plot_power([], Naris);

% %% Get the phase coherence metrics
% % create pairs of channels for detected events.
%
% for iSub = 1:length(PARAMS.Subjects)
%     sess_list = fieldnames(Events.(PARAMS.Subjects{iSub}));
%     for iSess = 1:length(sess_list)
%         [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), Coh_mat.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))]  = MS_event_pairs([], Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
%     end
% end
%
% %% plot the COH metrics
%
% stats_coh =  MS_Coh_plot_stats(Coh_mat);
%
%
%
% %% get the coordinates from the Expkeys
%
% % stats_subjects = MS_get_subject_info(data);
% fclose(PARAMS.log);
