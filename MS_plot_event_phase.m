function [all_mean] =  MS_plot_event_phase(cfg_in)
% %% for debugging
% run('/Users/jericcarmichael/Documents/GitHub/EC_Multisite/MS_initialize.m')
% global PARAMS
%
% % for USB data
% PARAMS.inter_dir = '/Volumes/Fenrir/MS_temp/';
% cd(PARAMS.inter_dir)

% Subjects = {'R102','R104','R107', 'R108', 'R112','R122','R123'};

%% default configurations
global PARAMS

cfg_def = [];
cfg_def.Subjects = {'R102','R104','R107', 'R108', 'R112','R122','R123'}; % can also be {'all'}
cfg_def.filter = [];
cfg_def.filter1.f = [45 65];
cfg_def.filter2.f = [70 90];
cfg_def.color.blue = double([158,202,225])/255;
cfg_def.color.green = double([168,221,181])/255;
cfg_def.color.OFC_NAc = double([123,225,160])/255;
% cfg_def.color.OFC_CG= double([255,168,213])/255;% old too bright
cfg_def.color.OFC_CG= double([158 1 66])/255;
cfg_def.lgd_size = 24;
cfg_def.measures = [3 10]; % which measures to use. Default is amplitude xcorr and phase slope. 
cfg_def.amp_xcorr_lags = -0.05:0.0005:0.05;
% If there is a specified measure or set of measures find the index

cfg = ProcessConfig2(cfg_def, cfg_in);

cfg.stats_dir = fopen([PARAMS.stats_dir  'Phase_Event_stats.txt'], 'w'); % can be used to append to a text file



%%
for iSub = cfg.Subjects
    fprintf('\nAnalyzing subject %s\n', iSub{1})
    if isunix
        load([PARAMS.inter_dir 'Phase_outputs/' iSub{1} '_phase_out.mat']);
    else
        load([PARAMS.inter_dir 'Phase_outputs\' iSub{1} '_phase_out.mat']);
    end
    if strcmp(iSub{1}, 'all')
        mat_in = mat_out;
        clear mat_out;
    else
        mat_in = Phase_mat;
        clear Phase_mat;
    end
    
    %% pre allocate all the sufields of the different analyses;
    sess_list = fieldnames(mat_in);
    
    labels = mat_in.(sess_list{1}).labels;
    bands = {'low', 'high'};
    all_mean = []; all_std = [];
    measures = fieldnames(mat_in.(sess_list{1}).(PARAMS.Phases{1}));
    
    
    for ii = size(labels,1):-1:1
        for jj = size(labels,2):-1:1
            Idx_mat(ii, jj) = str2double([num2str(ii) num2str(jj)]);
            
        end
        S = strsplit(labels{ii, jj},'_');
        single_labels{ii} = S{1};
    end
    
    %% move everything to the all_out structure
    sess_list = fieldnames(mat_in);
    
    for iPhase =1:length(PARAMS.Phases)
        for iMs = length(measures):-1:1 % go backwards to not overwrite AMP_AC_max and AMP_LAG_MAx
            if strcmp(measures{iMs}, 'EVT_COUNT')
                continue
            else
                for iBand = 1:2
                    if strcmp(measures{iMs}, 'AMP_LAG_max') || strcmp(measures{iMs}, 'AMP_AC_max')
                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}) = cell(size(labels));
                        continue
                    else
                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}) = cell(size(labels));
                        all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}) = cell(size(labels));
                        all_out = {};
                        for iSess = 1:length(sess_list)
                            all_out = cat(3,all_out, mat_in.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(measures{iMs}).(bands{iBand}));
                        end
                        for ii = size(all_out,1):-1:1
                            for jj = size(all_out,2):-1:1
                                % get the mean and median values
                                
                                temp = all_out(ii, jj,:);
                                %                     if isempty(temp{1,1,1})
                                %                         all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  [];
                                %                         all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  [];
                                %                     else
                                %                         if size(temp{1,1,1},1) > size(temp{1,1,1},2)
                                %                             temp(cellfun(@(temp) any(isnan(temp)),temp)) = [];
                                %                             temp = cell2mat(squeeze(temp)');
                                %                         else
                                temp(cellfun(@(temp) any(isnan(temp)),temp)) = [];
                                temp = cell2mat(squeeze(temp));
                                if isempty(temp) || size(temp,1) <10
                                    all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = [];
                                    all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  [];
                                    all_collect.AMP_LAG.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = [];
                                    all_collect.AMP_AC.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = [];
                                else
                                    if strcmp(measures{iMs}, 'AMP_AC') % collect all of the amp xcorr values
                                        for iEvt = size(temp,1):-1:1
                                            [all_collect.AMP_AC.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj}(iEvt),max_idx] = max(temp(iEvt,:));
                                            all_collect.AMP_LAG.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj}(iEvt) = cfg.amp_xcorr_lags(max_idx); % uses standard -0.05:0.0005:0.05
                                        end
                                    end
                                    
                                    %                         end
                                    if strcmp(measures{iMs}, 'Phase_lag_cxy')
                                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj}=  circ_mean(temp);
                                        all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  circ_std(temp);
                                        
                                    elseif strcmp(measures{iMs}, 'Phase_lag_mean')
                                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj}=  circ_mean(temp);
                                        all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  circ_std(temp);
                                        all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  temp;
                                    elseif strcmp(measures{iMs}, 'Phase_lag_F')  || strcmp(measures{iMs}, 'COH_fxx') || strcmp(measures{iMs}, 'PS_D')...
                                            || strcmp(measures{iMs}, 'AMP_AC') || strcmp(measures{iMs}, 'AMP_LAG');
                                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanmean(temp, 1);
                                        all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanstd(temp, [],1);
                                        
                                    elseif strcmp(measures{iMs}, 'EVT_COUNT')
                                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nansum(temp, 1);
                                        
                                    elseif strcmp(measures{iMs}, 'PS_slope')
                                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj}=  nanmean((temp./(2*pi))*1000); % ./(2*pi))*1000 convert to time rather than deg/Hz
                                        all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanstd((temp./(2*pi))*1000);
                                        
                                    else
                                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanmedian(temp,1);
                                        all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanstd(temp,[],1);
                                    end
                                    
                                    % work around for amp_lag max values
                                    if strcmp(measures{iMs}, 'AMP_AC')
                                        if isempty(temp)
                                            all_mean.AMP_AC_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = NaN;
                                            all_mean.AMP_LAG_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = NaN;
                                            all_std.AMP_AC_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = NaN;
                                            all_std.AMP_LAG_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = NaN;
                                        else
                                            lag = -0.05:0.0005:0.05;
                                            for iEvt = size(temp,1):-1:1
                                                [t_max(iEvt), idx] = max(temp(iEvt,:));
                                                lags_out(iEvt) = lag(idx);
                                            end
                                            all_mean.AMP_AC_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = nanmean(t_max);
                                            all_std.AMP_AC_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = nanmean(t_max);
                                            
                                            all_mean.AMP_LAG_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = nanstd(lags_out);
                                            all_std.AMP_LAG_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = nanstd(lags_out);
                                            clear lags_out t_max
                                        end
                                    end
                                end
                            end
                        end
                        clear temp
                    end
                end
            end
        end
    end
    clear all_out
end
%% get stats for AMP_xcorr max and lag [only used descriptives
stats_id = [PARAMS.stats_dir 'EVENT_AMP_stats.txt'];
stats_file = fopen(stats_id, 'w');


fprintf(stats_file,['**************** EVENT_AMP_stats ****************\n']);
fprintf(stats_file,['**************** ' date ' ****************\n']);
for iBand = 1:2
temp_con = all_collect.AMP_AC.contra.(bands{iBand});
temp_ipsi = all_collect.AMP_AC.ipsi.(bands{iBand});
sig_mat.AMP_AC.(bands{iBand}) = cell(size(temp_con));
for ii = size(temp_con,1):-1:1
    for jj = size(temp_con,2):-1:1
        if isempty(temp_con{ii,jj})
            temp_con{ii,jj} = NaN;
            desc_con(ii,jj) = NaN; desc_ipsi(ii,jj) = NaN;

        else
            % get mean for descriptives
            desc_con(ii,jj) = nanmean(temp_con{ii,jj});
            desc_ipsi(ii,jj) = nanmean(temp_ipsi{ii,jj});
            % get KS p values if needed. might use to add an astrix to the
            % sig plots or something
            if ~isempty(temp_con{ii,jj}) && ~isempty(temp_ipsi{ii,jj})
            [~,p] = kstest2(temp_con{ii,jj}, temp_ipsi{ii,jj});
            ks_mat.AMP_AC.(bands{iBand})(ii,jj) = p; 
            if p <0.001
                sig_mat.AMP_AC.(bands{iBand}){ii,jj} = '***';
            elseif p >=0.001 && p< 0.01
                sig_mat.AMP_AC.(bands{iBand}){ii,jj} = '**';
            elseif p >=0.01 && p<0.05
                sig_mat.AMP_AC.(bands{iBand}){ii,jj} = '*';
            else
                sig_mat.AMP_AC.(bands{iBand}){ii,jj} = 'n/s';
            end
            else
                sig_mat.AMP_AC.(bands{iBand}){ii,jj} = 'N/A';
                ks_mat.AMP_AC.(bands{iBand})(ii,jj) = NaN;
            end
        end
    end
end
fprintf(stats_file,['\n**************** CONTRA MEAN: ' bands{iBand} ' ****************\n']);
fprintf_table_EC(stats_file, desc_con, PARAMS.all_sites, PARAMS.all_sites);

fprintf(stats_file,['\n**************** IPSI MEAN: ' bands{iBand} ' ****************\n']);
fprintf_table_EC(stats_file, desc_con, PARAMS.all_sites, PARAMS.all_sites);

fprintf(stats_file,['\n**************** KS: ' bands{iBand} ' ****************\n']);
fprintf_table_EC(stats_file, ks_mat.AMP_AC.(bands{iBand}), PARAMS.all_sites, PARAMS.all_sites);

fprintf(stats_file,['\n**************** KS SIG: ' bands{iBand} ' ****************\n']);
fprintf_table_EC(stats_file, sig_mat.AMP_AC.(bands{iBand}), PARAMS.all_sites, PARAMS.all_sites);

end
fclose(stats_file);

%% get the average within the frequency bands

% EC this needs to be cleaned up.
for iPhase = 1:length(PARAMS.Phases)
    for iBand = 1:length(bands)
        for ii = size(labels,1):-1:1
            for jj = size(labels,2):-1:1
                
                % Coh
                if isempty(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    Coh_mean.(bands{iBand})(ii, jj) = NaN;
                    Coh_std.(bands{iBand})(ii, jj) = NaN;
                    % if it is in the lower triangular put the low
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    this_cxx = all_mean.COH_cxx.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    this_f = all_mean.COH_fxx.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    Coh_mean.(bands{iBand})(ii, jj) = nanmean(this_cxx(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                    Coh_std.(bands{iBand})(ii, jj) = nanstd(this_cxx(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                end
                
                % same for the phase lag
                if isempty(all_mean.Phase_lag_mean.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    Phase_diff.(bands{iBand})(ii, jj) = NaN;
                    Phase_diff_std.(bands{iBand})(ii, jj) = NaN;
                    % if it is in the lower triangular put the low
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    Phase_diff.(bands{iBand})(ii, jj) = all_mean.Phase_lag_mean.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    Phase_diff_std.(bands{iBand})(ii, jj) = all_std.Phase_lag_mean.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    
                end
                
                
                %                     if isempty(all_mean.EVT_COUNT.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                %                         Evt_count.(bands{iBand})(ii, jj) = NaN;
                %                         % if it is in the lower triangular put the low
                %                     elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                %                         Evt_count.(bands{iBand})(ii, jj) = nansum(all_mean.EVT_COUNT.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj});
                %                     end
                
                % PS
                if isempty(all_mean.PS_slope.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    PS_slope_mean.(bands{iBand})(ii, jj) = NaN;
                    PS_slope_std.(bands{iBand})(ii, jj) = NaN;
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    this_ps = all_mean.PS_slope.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    this_f = all_mean.PS_F.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    PS_slope_mean.(bands{iBand})(ii, jj) = nanmean(this_ps(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                    PS_slope_std.(bands{iBand})(ii, jj) = nanstd(this_ps(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                    
                end
                
                
                
                % AMP_Ac_MAx
                if isempty(all_mean.AMP_AC_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    AMP_AC_max_mean.(bands{iBand})(ii, jj) = NaN;
                    AMP_AC_max_std.(bands{iBand})(ii, jj) = NaN;
                    % if it is in the lower triangular put the low
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    AMP_AC_max_mean.(bands{iBand})(ii, jj) = all_mean.AMP_AC_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    AMP_AC_max_std.(bands{iBand})(ii, jj) = all_std.AMP_AC_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                end
                
                % amp lag
                if isempty(all_mean.AMP_LAG_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    AMP_LAG_max_mean.(bands{iBand})(ii, jj) = NaN;
                    AMP_LAG_max_std.(bands{iBand})(ii, jj) = NaN;
                    % if it is in the lower triangular put the low
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    AMP_LAG_max_mean.(bands{iBand})(ii, jj) = all_mean.AMP_LAG_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    AMP_LAG_max_std.(bands{iBand})(ii, jj) = all_std.AMP_LAG_max.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                end
            end
        end
    end
    mean_out.Coh_out.(PARAMS.Phases{iPhase}) = tril(Coh_mean.low) + tril(Coh_mean.high,-1)';
    mean_out.Phase_diff_out.(PARAMS.Phases{iPhase}) = tril(Phase_diff.low) + tril(Phase_diff.high,-1)';
    %         mean_out.Evt_count_out.(PARAMS.Phases{iPhase})  = tril(Evt_count.low) + tril(Evt_count.high,-1)';
    mean_out.AMP_AC_max_out.(PARAMS.Phases{iPhase})  = tril(AMP_AC_max_mean.low) + tril(AMP_AC_max_mean.high,-1)';
    mean_out.AMP_LAG_max_out.(PARAMS.Phases{iPhase})  = tril(AMP_LAG_max_mean.low) + tril(AMP_LAG_max_mean.high,-1)';
    mean_out.PS_slope_out.(PARAMS.Phases{iPhase})  = tril(PS_slope_mean.low) + tril(PS_slope_mean.high,-1)';
    
    % same for STD
    std_out.Coh_out.(PARAMS.Phases{iPhase}) = tril(Coh_std.low) + tril(Coh_std.high,-1)';
    std_out.Phase_diff_out.(PARAMS.Phases{iPhase}) = tril(Phase_diff_std.low) + tril(Phase_diff_std.high,-1)';
    %         mean_out.Evt_count_out.(PARAMS.Phases{iPhase})  = tril(Evt_count.low) + tril(Evt_count.high,-1)';
    std_out.AMP_AC_max_out.(PARAMS.Phases{iPhase})  = tril(AMP_AC_max_std.low) + tril(AMP_AC_max_std.high,-1)';
    std_out.AMP_LAG_max_out.(PARAMS.Phases{iPhase})  = tril(AMP_LAG_max_std.low) + tril(AMP_LAG_max_std.high,-1)';
    std_out.PS_slope_out.(PARAMS.Phases{iPhase})  = tril(PS_slope_std.low) + tril(PS_slope_std.high,-1)';
    
end
%% set up the colors to keep them consistent
low_mat = tril(ones(size(labels)),-1);
high_mat = triu(ones(size(labels)),1);

c_labels ={};
for ii = 1:length(labels)
    for jj = 1:length(labels)
        if low_mat(ii, jj)
            c_labels = cat(2,c_labels,labels{ii, jj});
        end
    end
end
c_ord = linspecer(length(c_labels));
%% cycle through measures to get plots (used for checks)
for iMs = cfg.measures
    % get the line plots for certain veriables
    if strcmp(measures{iMs}, 'Phase_lag_cxy') || strcmp(measures{iMs}, 'PS_slope') || strcmp(measures{iMs}, 'COH_cxx') || strcmp(measures{iMs}, 'AMP_AC')
        Phase_c = linspecer(4);
        low_mat = tril(ones(size(labels)),-1);
        figure(1111);
        for iBand = 1:length(bands)
            for iPhase =1:length(PARAMS.Phases)
                % plot the coherence
                subplot(2,2,iPhase)
                labs = {};
                handles_p = {};
                for ii = 1:size(labels,1)
                    for jj = 1:size(labels,2)
                        % select the bands and corresponding idx matrix.
                        if strcmp(bands{iBand}, 'low')
                            comp_mat = low_mat;
                        elseif strcmp(bands{iBand}, 'high')
                            comp_mat = high_mat';
                        end
                        
                        if isempty(strfind(labels{ii,jj}, 'Piri'))
                            if ~isempty(all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}) && comp_mat(ii, jj)
                                iC = find(ismember(c_labels, labels{ii,jj}));
                                hold on
                                % for the coherence
                                if strcmp(measures{iMs}, 'COH_cxx')
                                    h =shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}, all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}, all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj})));
                                    h.mainLine.Color = c_ord(iC,:);
                                    h.patch.FaceColor = c_ord(iC,:);
                                    h.patch.FaceAlpha = .5;
                                    handles_p = [handles_p, h.patch];
                                end
                                % for the phase lag
                                if strcmp(measures{iMs}, 'Phase_lag_cxy')
                                    sem =  rad2deg(all_std.Phase_lag_cxy.(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj})/sqrt(length(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}));
                                    h =shadedErrorBar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}),sem);
                                    h.mainLine.Color = c_ord(iC,:);
                                    h.patch.FaceColor = c_ord(iC,:);
                                    h.patch.FaceAlpha = .5;
                                    
                                    %                                 hE(iC) = plot(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), 'o', 'color', c_ord(iC,:), 'markerfacecolor', c_ord(iC,:));
                                    %                                 hE(iC) = errorbar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), sem, 'o', 'color', c_ord(iC,:))%, 'CapSize',2);
                                    handles_p = [handles_p, h.patch];
                                end
                                % for the phase slope
                                if strcmp(measures{iMs}, 'PS_slope')
                                    h =shadedErrorBar(all_mean.PS_F.(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}, all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}, all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj})));
                                    h.mainLine.Color = c_ord(iC,:);
                                    h.patch.FaceColor = c_ord(iC,:);
                                    h.patch.FaceAlpha = .5;
                                    handles_p = [handles_p, h.patch];
                                end
                                
                                if strcmp(measures{iMs}, 'AMP_AC')
                                    h =shadedErrorBar(all_mean.AMP_LAG.(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}, all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}, all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii, jj})));
                                    h.mainLine.Color = c_ord(iC,:);
                                    h.patch.FaceColor = c_ord(iC,:);
                                    h.patch.FaceAlpha = .5;
                                    handles_p = [handles_p, h.patch];
                                end
                                
                                
                                labs = [labs, labels{ii, jj}];
                                iC = iC+1;
                            end
                        end
                    end
                end
                if strcmp(measures{iMs}, 'PS_slope')
                    y_lim = [-20 20];
                    ylabel('Phase offset (ms)');
                    xlim([30 100])
                    box_pos =  [0, y_lim(1), 100, sum(abs(y_lim))];
                    gamma_box_h = 2;
                    title_pos = [80  8 1.00011];
                    
                elseif strcmp(measures{iMs}, 'Phase_lag_cxy')
                    y_lim = [-200 200];
                    ylabel('Phase lag (deg)');
                    xlim([0 100])
                    box_pos =  [0, -200, 100, 400];
                    gamma_box_h = 10;
                    title_pos = [80  160 1.00011];
                elseif strcmp(measures{iMs}, 'AMP_AC')
                    y_lim =[0 1];
                    ylabel('Amplitude xcor');
                    xlim([-0.05 0.05])
                    box_pos =  [-0.05, 0, .1, 1];
                    gamma_box_h = 0.05;
                    title_pos = [0.04  0.9 1.00011];
                else
                    y_lim = [0 1];
                    ylabel('Coherence');
                    xlim([0 100])
                    box_pos =  [0, 0, 100, 1];
                    gamma_box_h = 0.05;
                    title_pos = [80  0.9 1.00011];
                    
                end
                title(PARAMS.Phases{iPhase}, 'color', Phase_c(iPhase,:),'fontsize', 32,'Position',title_pos) % add the phase title
                
                if iPhase ==4
                    hL = legend(handles_p, labs);
                    set(hL,'Position', [0.45 0.5 0.05 0.1], 'Units', 'normalized')
                end
                ylim(y_lim);
                %                 text(.5,1.02,PARAMS.Phases{iPhase},'units','normalized','fontsize', 32,'horizontalalignment','center','verticalalignment','top', 'color', Phase_c(iPhase,:));
                xlabel('Frequency (Hz)');
                % put boxes for the gamma bands
                rectangle('position', [cfg.filter1.f(1), y_lim(1), (cfg.filter1.f(2)-cfg.filter1.f(1)), gamma_box_h], 'FaceColor',[0.5 0.5 0.5], 'edgecolor', [0.5 0.5 0.5])
                rectangle('position', [cfg.filter2.f(1), y_lim(1), (cfg.filter2.f(2)-cfg.filter2.f(1)), gamma_box_h], 'FaceColor',[0.7 0.7 0.7], 'edgecolor', [0.7 0.7 0.7])
                rectangle('position',box_pos, 'EdgeColor', Phase_c(iPhase,:), 'linewidth', 2)
            end
            Square_subplots
            tightfig
            cfg_count_fig.resize = 1;
            cfg_count_fig.pos = [100 40 1190 765];
            SetFigure(cfg_count_fig, gcf)
            %             for  iPhase = 1:4
            %                 subplot(2,2,iPhase)
            %                 hold on
            %                 text(.5,1.02,PARAMS.Phases{iPhase},'units','normalized','fontsize', 32,'horizontalalignment','center','verticalalignment','top', 'color', Phase_c(iPhase,:));
            %             end
            mkdir(PARAMS.inter_dir, 'Phase_plots')
            saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '_' bands{iBand} '.fig']);
            saveas_eps([iSub{1} '_' measures{iMs} '_' bands{iBand}], [PARAMS.inter_dir 'Phase_plots/'])
            %             saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs}], 'epsc')
            print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '_' bands{iBand} '.png'],'-r300')
            close all
        end
    end
end

%% same thing but each site-pair gets a single plot.  (used for figure 3)
%% cycle through measures to get plots
for iMs = cfg.measures%[3 10]%1:length(measures)
    % get the line plots for certain veriables
    if strcmp(measures{iMs}, 'Phase_lag_cxy') || strcmp(measures{iMs}, 'PS_slope') || strcmp(measures{iMs}, 'COH_cxx') || strcmp(measures{iMs}, 'AMP_AC')
        
        for iBand = 1:length(bands)
            low_mat = tril(ones(size(labels)),-1);
            for ii = 1:size(labels,1)
                for jj = 1:size(labels,2)
                    % select the bands and corresponding idx matrix.
                    
                    if strcmp(bands{iBand}, 'low')
                        comp_mat = low_mat;
                        mean_ii = jj; % used for getting the mean and STD values.  Set up as a upper/lower matrix earlier. Pain to dig out.
                        mean_jj = ii;
                    elseif strcmp(bands{iBand}, 'high')
                        comp_mat = high_mat';
                        mean_ii = ii; % swap so that it is using the upper of the matrix which is for high gamma.
                        mean_jj = jj;
                    end
                    
%                     if isempty(strfind(labels{ii,jj}, 'Piri'))
                        if ~isempty(all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj}) %&& comp_mat(ii, jj)
                            % for the coherence
                            h_curr= figure();
                            hold on
                            if strcmp(labels{ii,jj}, 'CG_OFC') || strcmp(labels{ii,jj}, 'OFC_CG')
                                this_color = cfg_def.color.OFC_CG;
                            elseif strcmp(labels{ii,jj}, 'NAc_OFC') || strcmp(labels{ii,jj}, 'OFC_NAc')
                                this_color = cfg_def.color.OFC_NAc;
                            else
                                if sum(ismember(c_labels, labels{ii, jj}))
                                    this_color = c_ord(ismember(c_labels, labels{ii, jj}),:);
                                else
                                    this_color = c_ord(ismember(c_labels, labels{jj, ii}),:);  % get the colour for that pair by swapping order and finding.
                                end
                            end
                            % Coherence (not used)
%                             if strcmp(measures{iMs}, 'COH_cxx')
%                                 h =shadedErrorBar(all_mean.COH_fxx.contra.(bands{iBand}){ii, jj}, all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj}, all_std.(measures{iMs}).contra.(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj})));
%                                 h.mainLine.Color = this_color;
%                                 h.patch.FaceColor = this_color;
%                                 h.patch.FaceAlpha = .5;
%                                 h.patch.EdgeColor = this_color;
%                                 h.edge(1).Color = this_color;
%                                 h.edge(2).Color = this_color;
%                                 if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
%                                     h2 =shadedErrorBar(all_mean.COH_fxx.ipsi.(bands{iBand}){ii, jj}, all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj}, all_std.(measures{iMs}).ipsi.(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})));
%                                     h2.mainLine.Color = [.6 .6 .6];
%                                     h2.patch.FaceColor = [.6 .6 .6];
%                                     h2.patch.FaceAlpha = .5;
%                                     h2.patch.EdgeColor = [.6 .6 .6];
%                                     h2.edge(1).Color = [.6 .6 .6];
%                                     h2.edge(2).Color = [.6 .6 .6];
%                                 end
%                                 y_lim = [0 1];
%                                 ylabel('Coherence'); xlabel('Frequency')
%                                 x_lim=([0 100]);
%                                 box_pos =  [0, 0, 100, 1];
%                                 gamma_box_h = 1;
%                                 
%                             end
%                             % for the phase lag (not used)
%                             if strcmp(measures{iMs}, 'Phase_lag_cxy')
%                                 sem =  rad2deg(all_std.Phase_lag_cxy.contra.(bands{iBand}){ii, jj})/sqrt(length(all_mean.Phase_lag_cxy.contra.(bands{iBand}){ii, jj}));
%                                 h =shadedErrorBar(all_mean.Phase_lag_F.contra.(bands{iBand}){ii, jj}, rad2deg(all_mean.Phase_lag_cxy.contra.(bands{iBand}){ii, jj}),sem);
%                                 h.mainLine.Color = this_color;
%                                 h.patch.FaceColor = this_color;
%                                 h.patch.FaceAlpha = .5;
%                                 h.patch.EdgeColor = this_color;
%                                 h.edge(1).Color = this_color;
%                                 h.edge(2).Color = this_color;
%                                 
%                                 if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
%                                     sem =  rad2deg(all_std.Phase_lag_cxy.ipsi.(bands{iBand}){ii, jj})/sqrt(length(all_mean.Phase_lag_cxy.ipsi.(bands{iBand}){ii, jj}));
%                                     h2 =shadedErrorBar(all_mean.Phase_lag_F.ipsi.(bands{iBand}){ii, jj}, rad2deg(all_mean.Phase_lag_cxy.ipsi.(bands{iBand}){ii, jj}),sem);
%                                     h2.mainLine.Color = [.6 .6 .6];
%                                     h2.patch.FaceColor = [.6 .6 .6];
%                                     h2.patch.FaceAlpha = .5;
%                                     h2.patch.EdgeColor = [.6 .6 .6];
%                                     h2.edge(1).Color = [.6 .6 .6];
%                                     h2.edge(2).Color = [.6 .6 .6];
%                                 end
%                                 
%                                 %                                 hE(iC) = plot(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), 'o', 'color', c_ord(iC,:), 'markerfacecolor', c_ord(iC,:));
%                                 %                                 hE(iC) = errorbar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), sem, 'o', 'color', c_ord(iC,:))%, 'CapSize',2);
%                                 
%                                 y_lim = [-200 0 200];
%                                 ylabel('Phase lag (deg)');
%                                 x_lim=([0 100]);
%                                 box_pos =  [0, -200, 100, 400];
%                                 gamma_box_h = 400;
%                                 
%                                 
%                             end
                            % for the phase slope
                            if strcmp(measures{iMs}, 'PS_slope')
                                
                                hold on
                                if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
                                    h2 =shadedErrorBar(all_mean.PS_F.ipsi.(bands{iBand}){ii, jj}, all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj}, all_std.(measures{iMs}).ipsi.(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})));
                                    h2.mainLine.Color = [.8 .8 .8];
                                    h2.mainLine.LineWidth = 3;
                                    h2.patch.FaceColor = [.6 .6 .6];
                                    h2.patch.FaceAlpha = .5;
                                    h2.patch.EdgeColor = [.6 .6 .6];
                                    h2.edge(1).Color = [.6 .6 .6];
                                    h2.edge(2).Color = [.6 .6 .6];
                                end
                                
                                % swap order or else it is hard to see
                                h =shadedErrorBar(all_mean.PS_F.contra.(bands{iBand}){ii, jj}, all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj}, all_std.(measures{iMs}).contra.(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj})));
                                h.mainLine.Color = this_color;
                                h.mainLine.LineWidth = 3;
                                h.patch.FaceColor = this_color;
                                h.patch.FaceAlpha = .7;
                                h.patch.EdgeColor = this_color;
                                h.edge(1).Color = this_color;
                                h.edge(2).Color = this_color;
                                
                                y_lim = [-15, 0, 15];
                                ylabel('Phase offset (ms)'); xlabel('Frequency');
                                x_lim=([30 100]);
                                box_pos =  [30, -5, 100, 40];
                                gamma_box_h = 2;
                                
                                if ~isempty(all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj})
                                    % add text for mean PS value for each gamma
                                    this_ps = all_mean.PS_slope.contra.(bands{iBand}){ii,jj};
                                    this_f = all_mean.PS_F.contra.(bands{iBand}){ii,jj};
                                    contra_mean = nanmean(this_ps(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                                    contra_std = nanstd(this_ps(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                                    % add text
                                    mean_txt = text(60, 13, ['mean:' num2str(contra_mean,'%.2f') ' ' char(177) ' ' ...
                                        num2str(contra_std,'%.2f')], 'color', 'k', 'fontsize', 48, 'FontName', 'helvetica');
                                    
                                    % just for text in paper (get 55-70Hz) mean for CG-OFC
                                    if (iMs ==10 && iBand ==1 && strcmp(labels{ii,jj}, 'CG_OFC')) || (iMs ==10 && iBand ==1 &&strcmp(labels{ii,jj}, 'OFC_CG')) 
                                        contra_mid_mean = nanmean(this_ps(nearest_idx(55, this_f):nearest_idx(70, this_f)));
                                        contra_mid_std = nanstd(this_ps(nearest_idx(55, this_f):nearest_idx(70, this_f)));
                                        fprintf('\nOFC-CG 55-70Hz mean %.2fms sd %.2fms \n', contra_mid_mean, contra_mid_std)
                                    
                                    end
                                % find the ipsi
                                if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
                                    this_ps = all_mean.PS_slope.ipsi.(bands{iBand}){ii,jj};
                                    this_f = all_mean.PS_F.ipsi.(bands{iBand}){ii,jj};
                                    ipsi_mean = nanmean(this_ps(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                                    ipsi_std = nanstd(this_ps(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                                    % add text
                                    mean_txt_ip = text(60, 10, ['mean:' num2str(ipsi_mean,'%.2f') ' ' char(177) ' ' ...
                                        num2str(ipsi_std,'%.2f')], 'color', [.7 .7 .7], 'fontsize', 48, 'FontName', 'helvetica');
                                end
                            end
                            
                            if strcmp(measures{iMs}, 'AMP_AC')
                                
                                hold on
                                
                                if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
                                    h2 =shadedErrorBar(all_mean.AMP_LAG.ipsi.(bands{iBand}){ii, jj}*1000, all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj}, all_std.(measures{iMs}).ipsi.(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})));
                                    h2.mainLine.Color = [.6 .6 .6];
                                    h2.patch.FaceColor = [.6 .6 .6];
                                    h2.patch.FaceAlpha = .5;
                                    h2.patch.EdgeColor = [.6 .6 .6];
                                    h2.edge(1).Color = [.6 .6 .6];
                                    h2.edge(2).Color = [.6 .6 .6];
                                    [ip_val, ip_idx] = max(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj});
                                end
                                h =shadedErrorBar(all_mean.AMP_LAG.contra.(bands{iBand}){ii, jj}*1000, all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj}, all_std.(measures{iMs}).contra.(bands{iBand}){ii, jj}/sqrt(length(all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj})));
                                h.mainLine.Color = this_color;
                                h.patch.FaceColor = this_color;
                                h.patch.FaceAlpha = .5;
                                h.patch.EdgeColor = this_color;
                                h.edge(1).Color = this_color;
                                h.edge(2).Color = this_color;
                                [con_val, con_idx] = max(all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj});
                                
                                
                                y_lim =[0 1];
                                ylabel('Amplitude xcor');
                                x_lim=([-50 0 50]);
                                xlabel('time (ms)')
                                box_pos =  [-0.05, 0, .1, 1];
                                gamma_box_h = 1;
                                ylim([y_lim(1) y_lim(end)]);xlim([x_lim(1) x_lim(end)])
                                h_con = vline(all_mean.AMP_LAG.contra.(bands{iBand}){ii, jj}(con_idx));
                                h_con.Color = this_color; h_con.LineStyle = '- -'; h_con.LineWidth = 3;
                                
                                
                                if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
                                    h_ip = [];
                                    h_ip = vline(all_mean.AMP_LAG.contra.(bands{iBand}){ii, jj}(ip_idx));
                                    h_ip.Color = [.6 .6 .6]; h_ip.LineStyle =  ':'; h_ip.LineWidth = 3;
                                    lgd = legend([h_con, h_ip],['Contra lag: ' num2str(all_mean.AMP_LAG.contra.(bands{iBand}){ii, jj}(con_idx)*1000,2) 'ms'],...
                                        [ 'Ipsi lag:' num2str(all_mean.AMP_LAG.ipsi.(bands{iBand}){ii, jj}(ip_idx)*1000) 'ms'], 'location', 'northwest');
                                else
                                    lgd = legend(h_con,['Contra lag: ' num2str(all_mean.AMP_LAG.contra.(bands{iBand}){ii, jj}(con_idx)*1000,2) 'ms'],'location', 'northwest');
                                end
                                
                                lgd.FontSize = cfg.lgd_size ; lgd.Box = 'off';
                                
                                % add text for means SD
                                %   % for contra
                                [max_val, max_id] = max(all_mean.(measures{iMs}).contra.(bands{iBand}){ii, jj});
                                mean_txt = text(2, .94, ['mean:' num2str(max_val,'%.2f') ' ' char(177) ' ' ...
                                    num2str(all_std.(measures{iMs}).contra.(bands{iBand}){ii, jj}(max_id),'%.2f')], 'color', 'k', 'fontsize', cfg.lgd_size, 'FontName', 'helvetica');
                                % for ipsi
                                [max_val, max_id] = max(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj});
                                mean_txt_ip = text(2, .85, ['mean:' num2str(max_val,'%.2f') ' ' char(177) ' ' ...
                                    num2str(all_std.(measures{iMs}).ipsi.(bands{iBand}){ii, jj}(max_id),'%.2f')], 'color', [.7 .7 .7], 'fontsize', cfg.lgd_size, 'FontName', 'helvetica');
                                
                            end
                            
                            ylim([y_lim(1) y_lim(end)]);xlim([x_lim(1) x_lim(end)])
                            if strcmp(measures{iMs}, 'AMP_AC')
                                v_zero = vline(0, 'k');
                            else
                                % put boxes for the gamma bands
                                rectangle('position', [cfg.filter1.f(1), y_lim(1), (cfg.filter1.f(2)-cfg.filter1.f(1)), gamma_box_h],  'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3])
                                rectangle('position', [cfg.filter2.f(1), y_lim(1), (cfg.filter2.f(2)-cfg.filter2.f(1)), gamma_box_h],  'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
                                chi=get(gca, 'Children');
                                set(gca, 'Children',flipud(chi))
                            end
                            set(findall(gca, 'Type', 'Line'),'LineWidth',2)
                            if strcmp(measures{iMs}, 'AMP_AC')
                                if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
                                h_ip.LineWidth = 3;
                                end
                                h_con.LineWidth = 3;
                                v_zero.LineWidth = 1;
                                uistack(v_zero, 'bottom')
                            end
                            set(gca, 'ytick', y_lim, 'xtick',x_lim)
                            %                                 ylabel([]); xlabel([]); %remove the labels for clerity in the figure.
                            
                            
                            cfg.set_fig.ft_size = 36;
                            SetFigure(cfg.set_fig, gcf)
                            
                            mkdir(PARAMS.inter_dir, 'Phase_plots')
                            % save the version with the ipsi for the
                            % suppliment
                            saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '_' labels{ii,jj} '_' bands{iBand} 'S2.fig']);
                            saveas_eps([iSub{1} '_' measures{iMs} '_' labels{ii,jj} '_' bands{iBand} 'S2'], [PARAMS.inter_dir 'Phase_plots/'])
                            print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '_' labels{ii,jj} '_' bands{iBand} 'S2.png'],'-r300')
                            
                            
                            if strcmp(labels{ii,jj}, 'NAc_OFC') || strcmp(labels{ii,jj}, 'OFC_NAc') || strcmp(labels{ii,jj}, 'CG_OFC') || strcmp(labels{ii,jj}, 'OFC_CG')

                            % same thing but without the ipsi line
                            if ~isempty(all_mean.(measures{iMs}).ipsi.(bands{iBand}){ii, jj})
                                set(h2.patch, 'visible', 'off'); set(h2.edge, 'visible', 'off'); set(h2.mainLine, 'visible', 'off')
                                if strcmp(measures{iMs}, 'AMP_AC')
                                    h_ip.Visible = 'off';
                                    lgd = legend(h_con,['Contra lag: ' num2str(all_mean.AMP_LAG.contra.(bands{iBand}){ii, jj}(con_idx)*1000,2) 'ms']);
                                    lgd.FontSize = 24; lgd.Box = 'off';
                                    clear h_ip
                                end
                            end
                            % remove ipsi text
                            
                            mean_txt_ip.Visible = 'off';
                            
                            
                            set(gca, 'ytick', y_lim, 'xtick',x_lim)
                            cfg.set_fig.ft_size = 36;
                            SetFigure(cfg.set_fig, gcf)
                            ylim([y_lim(1) y_lim(end)]);xlim([x_lim(1) x_lim(end)])
                            %                                     ylabel([]); xlabel([]); %remove the labels for clerity in the figure.
                            
                            %                                 rectangle('position', [cfg.filter1.f(1), y_lim(1), (cfg.filter1.f(2)-cfg.filter1.f(1)), gamma_box_h],  'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3])
                            %                                 rectangle('position', [cfg.filter2.f(1), y_lim(1), (cfg.filter2.f(2)-cfg.filter2.f(1)), gamma_box_h],  'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
                            chi=get(gca, 'Children');
                            set(gca, 'Children',flipud(chi))
                            
                            saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '_' labels{ii,jj} '_' bands{iBand} 'F3.fig']);
                            saveas_eps([iSub{1} '_' measures{iMs} '_' labels{ii,jj} '_' bands{iBand} 'F3'], [PARAMS.inter_dir 'Phase_plots/'])
                            print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '_' labels{ii,jj} '_' bands{iBand} 'F3.png'],'-r300')
                            end
                            close all
                        end
                        
%                     end
                end
            end
        end
    end
end

%% plot the single values for each pair in matrix format
mean_measures = fieldnames(mean_out);
for iMs = 1:length(mean_measures)
    
    Phase_c = linspecer(4);
    low_mat = tril(ones(size(labels)),-1);
    figure(1113);
    for iPhase =1:length(PARAMS.Phases)
        
        subplot(1,4,iPhase)
        
        if strcmp(mean_measures{iMs}, 'Coh_out')
            caxis([0 1]);
            nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}))
            [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}), 2,  10); % adds numerics to imagesc
        elseif strcmp(mean_measures{iMs}, 'Phase_diff_out')
            caxis([0 180]);
            nan_imagesc_ec(abs(rad2deg(mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}))))
            [~, hStrings] =add_num_imagesc(gca,rad2deg(mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase})), 0,  10); % adds numerics to imagesc
        elseif strcmp(mean_measures{iMs}, 'AMP_AC_max_out')
            caxis([0 1]);
            nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}))
            [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}), 2,  10); % adds numerics to imagesc
        elseif strcmp(mean_measures{iMs}, 'AMP_LAG_max_out')
            caxis([-0.02 0.02]);
            nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}))
            [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}), 4,  10); % adds numerics to imagesc
        elseif strcmp(mean_measures{iMs}, 'Evt_count_out')
            caxis([0 4000])
            nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}))
            [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{iPhase}), 0,  10); % adds numerics to imagesc
        end
        
        %     hcb = colorbar();
        colormap('parula')
        
        %     set(hcb, 'Ytick', 0:0.1:.5)
        set(gca,'xtick', 1:length(single_labels), 'ytick', 1:length(single_labels),  'xticklabel', single_labels, 'yticklabel', single_labels)
        set(gca,'XTickLabelRotation',65)
        rectangle('position', [0.5, 0.5, 7, 7], 'EdgeColor', Phase_c(iPhase,:), 'linewidth', 4)
        
    end
    %
    Square_subplots
    tightfig
    cfg_count_fig.resize = 1;
    cfg_count_fig.pos = [3 401 1440 400];
    SetFigure(cfg_count_fig, gcf)
    
    saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' mean_measures{iMs} '.fig']);
    %         saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' mean_measures{iMs}], 'epsc')
    saveas_eps([iSub{1} '_' mean_measures{iMs}], [PARAMS.inter_dir 'Phase_plots/'])
    print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' mean_measures{iMs} '.png'], '-r300')
    
end
%% use a difference matrix
mean_measures = fieldnames(mean_out);
ft_size = 18;
for iMs = 1:length(mean_measures)
    figure(1114);
    if strcmp(mean_measures{iMs}, 'Coh_out')
        caxis([0 1]);
        nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})- mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}))
        [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})- mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}), 2,  ft_size); % adds numerics to imagesc
    elseif strcmp(mean_measures{iMs}, 'Phase_diff_out')
        caxis([0 180]);
        nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})- mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}))
        [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})- mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}), 2,  ft_size); % adds numerics to imagesc
    elseif strcmp(mean_measures{iMs}, 'AMP_AC_max_out')
        caxis([0 1]);
        nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})-mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}))
        [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})-mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}), 2,  ft_size); % adds numerics to imagesc
    elseif strcmp(mean_measures{iMs}, 'AMP_LAG_max_out')
        caxis([-0.02 0.02]);
        nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})-mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}))
        [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})-mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}), 4,  ft_size); % adds numerics to imagesc
    elseif strcmp(mean_measures{iMs}, 'Evt_count_out')
        caxis([0 4000])
        nan_imagesc_ec(mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})- mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}))
        [~, hStrings] =add_num_imagesc(gca,mean_out.(mean_measures{iMs}).(PARAMS.Phases{3})- mean_out.(mean_measures{iMs}).(PARAMS.Phases{2}), 2,  ft_size); % adds numerics to imagesc
        
    end
    
    %     hcb = colorbar();
    colormap('parula')
    
    %     set(hcb, 'Ytick', 0:0.1:.5)
    set(gca,'xtick', 1:length(single_labels), 'ytick', 1:length(single_labels),  'xticklabel', single_labels, 'yticklabel', single_labels)
    set(gca,'XTickLabelRotation',65)
    %                 tightfig
    %         cfg_count_fig.resize = 1;
    %         cfg_count_fig.pos = [3 401 1440 400];
    SetFigure([], gcf)
    
    saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' mean_measures{iMs} '_diff.fig']);
    %         saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' mean_measures{iMs} '_diff'], 'epsc')
    saveas_eps([iSub{1} '_' mean_measures{iMs} '_diff'], [PARAMS.inter_dir 'Phase_plots/'])
    print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' mean_measures{iMs} '_diff.png'], '-r300')
    close all
end


%%
stats_log = fopen([PARAMS.inter_dir iSub{1} '_Phase_stats.txt'], 'w');
ver = version;
ver = str2double(ver(end-5:end-2));
for iPhase =1:length(PARAMS.Phases)
    %         phase_mean_mat.(PARAMS.Phases{iPhase}) = NaN(length(labels), length(labels));
    str_vals.low = cell(7,7); str_vals.high = cell(7,7);
    str_vals_piri.low = cell(7,7); str_vals_piri.high = cell(7,7);
    phase_mean = NaN(7,7); vector_r = NaN(7,7);
    
    pol_legend = {}; piri_legend = {};
    n = 1000;
    figure(1115)
    all_h = [];
    %         for kk = 1:4
    %            hold on
    %            subplot(2,2,kk)
    %            h_temp = polar([0 0], [0 1]);
    %            set(h_temp, 'visible', 'off') ;
    %
    %         end
    %%
    for ii = 1:size(labels,1)
        for jj = 1:size(labels,2)
            iC = find(ismember(c_labels, labels{ii,jj}));
            if ~isnan(mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj)) && low_mat(ii,jj)
                iC = find(ismember(c_labels, labels{ii,jj}));
                
                phase_mean(ii,jj) = circ_mean(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).low{ii,jj});
                vector_r(ii,jj) = circ_r(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).low{ii,jj});
                %                     spacing_str = '________';
                %                     spacing_str =  spacing_str(1:end-length(labels{ii,jj}))
                
                %                     str_vals.low{ii,jj} = sprintf('low %s: %.0f ,MVL .2f',labels{ii,jj}, rad2deg(phase_mean), vector_r); % with values.
                %                     str_vals.low{ii, jj} = labels{ii, jj}; % just the labels
                if isempty(strfind(labels{ii,jj}, 'Piri'))
                    subplot(2,2,1)
                    pol_legend = [pol_legend, labels{ii,jj}];
                    str_vals.low{ii,jj} = sprintf('%s: & %.0f & %.2f',labels{ii,jj}, rad2deg(phase_mean(ii,jj)), vector_r(ii,jj)); % with values.
                else
                    subplot(2,2,3)
                    piri_legend = [piri_legend, labels{ii,jj}];
                    str_vals_piri.low{ii,jj} = sprintf('%s: & %.0f & %.2f',labels{ii,jj}, rad2deg(phase_mean(ii,jj)), vector_r(ii,jj)); % with values.
                    
                end
                % save the values for a summary table later.
                %                     if isempty(phase_mean)
                %                         phase_mean_mat.(PARAMS.Phases{iPhase})(ii, jj) = NaN;
                %                         phase_vec_mat.(PARAMS.Phases{iPhase})(ii, jj) = NaN;
                %                     else
                %                         phase_mean_mat.(PARAMS.Phases{iPhase})(ii, jj) = phase_mean;
                %                         phase_vec_mat.(PARAMS.Phases{iPhase})(ii, jj) = vector_r;
                %                     end
                
                %                 [t,r]= rose(mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj),n);
                
                %                 h = polar(t,r);
                
                %polar plot the median and the vector length
                if ver >= 2017
                    h = polarplot([0 phase_mean(ii,jj)],[0 vector_r(ii,jj)]);
                    rlim([0 1])
                    set(gca, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top','RTick' ,[0.5 1], 'FontSize', 22, 'FontName', 'Helvetica')
                    cfg_polar.axes_labels = 0; % turns off any scaling of the x and y label which don't work in polar plots
                else
                    h_temp = polar([0 0], [0 1]);
                    hold on
                    h = polar([0 phase_mean(ii,jj)], [0 vector_r(ii,jj)]);
                    set(h_temp, 'visible', 'off') ;
                    cfg_polar.axes_labels = 1; % turns off any scaling of the x and y label which don't work in polar plots
                end
                
                S = strsplit(labels{ii,jj}, '_');
                if strcmp(S{1}, 'NAc') && strcmp(S{2}, 'OFC') || strcmp(S{1}, 'PiriN') && strcmp(S{2}, 'PiriO')
                    set(h, 'color', 'k', 'linewidth', 7)
                else
                    set(h, 'color', c_ord(iC,:), 'linewidth', 4)
                end
                hold all
                %                     view([90, 90])
                
            elseif ~isnan(mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj)) && high_mat(ii,jj) && ~isempty(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj})
                iC = find(ismember(c_labels, labels{jj,ii}));
                
                phase_mean(ii,jj) = circ_mean(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj});
                vector_r(ii,jj) = circ_r(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj});
                
                if isempty(strfind(labels{ii,jj}, 'Piri')) % contains no piriform channel
                    subplot(2,2,2)
                    str_vals.high{ii,jj} = sprintf('%s: & %.0f  & %.2f',labels{jj,ii}, rad2deg(phase_mean(ii,jj)), vector_r(jj,ii)); % with values.
                    
                else % if it is a piriform - x comparison
                    subplot(2,2,4)
                    str_vals_piri.high{ii,jj} = sprintf('%s: & %.0f  & %.2f',labels{jj,ii}, rad2deg(phase_mean(ii,jj)), vector_r(jj,ii)); % with values.
                end
                %                 [t,r]= rose(mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj),n);
                %
                %                 h = polar(t,r);
                %                 phase_mean = circ_mean(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj});
                %                 vector_r = circ_r(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj});
                %                 h = compass(phase_mean,vector_r);
                
                
                
                % save the values for a summary table later.
                %                     if isempty(phase_mean)
                %                         phase_mean_mat.(PARAMS.Phases{iPhase})(ii, jj) = NaN;
                %                         phase_vec_mat.(PARAMS.Phases{iPhase})(ii, jj) = NaN;
                %                     else
                %                         phase_mean_mat.(PARAMS.Phases{iPhase})(ii, jj) = phase_mean;
                %                         phase_vec_mat.(PARAMS.Phases{iPhase})(ii, jj) = vector_r;
                %                     end
                %polar plot the median and the vector length
                if ver >= 2017
                    h = polarplot( [0 phase_mean(ii,jj)],[0 vector_r(ii,jj)]);
                    rlim([0 1])
                    set(gca, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top','RTick' ,[0.5 1], 'FontSize', 22, 'FontName', 'Helvetica')
                    cfg_polar.axes_labels = 0; % turns off any scaling of the x and y label which don't work in polar plots
                else
                    h_temp = polar([0 0], [0 1]);
                    hold on
                    h = polar([0 phase_mean(ii,jj)], [0 vector_r(ii,jj)]);
                    set(h_temp, 'visible', 'off') ;
                    cfg_polar.axes_labels = 1; % turns off any scaling of the x and y label which don't work in polar plots
                end
                
                
                S = strsplit(labels{ii,jj}, '_');
                if strcmp(S{2}, 'NAc') && strcmp(S{1}, 'OFC') || strcmp(S{1}, 'PiriN') && strcmp(S{2}, 'PiriO')
                    set(h, 'color', 'k', 'linewidth', 7)
                else
                    set(h, 'color', c_ord(iC,:), 'linewidth', 4)
                end
                hold all
                %                     view([270, 90])
                %                     view([0, 90])
            end
        end
    end
    %     text(.5,1.15,'low','units','normalized','fontsize', 24,'horizontalalignment','center','verticalalignment','top');
    %         text(2,1.15,'high','units','normalized','fontsize', 24,'horizontalalignment','center','verticalalignment','top');
    
    th = findall(gcf,'Type','text');
    for i = 1:length(th)
        set(th(i),'FontSize',20)
    end
    subplot(2,2,1)
    hL = legend(pol_legend);
    set(hL,'Position', [0.5 0.7 0.05 0.1], 'Units', 'normalized')
    
    subplot(2,2,3)
    hL = legend(piri_legend);
    set(hL,'Position', [0.5 0.25 0.05 0.1], 'Units', 'normalized')
    %         tightfig
    cfg_polar.resize = 1;
    cfg_polar.pos = [205    50  1200   800];
    SetFigure(cfg_polar, gcf);
    
    saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' (PARAMS.Phases{iPhase}) '_polar_phase.fig']);
    %         saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' (PARAMS.Phases{iPhase}) '_polar_phase'], 'epsc')
    saveas_eps([iSub{1} '_' (PARAMS.Phases{iPhase}) '_polar_phase'], [PARAMS.inter_dir 'Phase_plots/'])
    print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' (PARAMS.Phases{iPhase})  '_polar_phase.png'], '-r300')
    close all
    
    %% print the statstistics
    fprintf(stats_log,['\n\n' PARAMS.Phases{iPhase} '\n\n'])
    
    bands = {'low', 'high'};
    for iBands = 1:length(bands)
        fprintf(stats_log, '####    Circ Stats for %s gamma    ####\n', bands{iBands})
        for ii = 1:size(str_vals.(bands{iBands}),1)
            for jj = 1:size(str_vals.(bands{iBands}),2)
                if ~isempty(str_vals.(bands{iBands}){ii, jj})
                    fprintf(stats_log,str_vals.(bands{iBands}){ii, jj})
                    fprintf(stats_log,'\n')
                end
                
            end
        end
    end
    % same for the piriform pairs
    fprintf(stats_log,'-------- Trans piriform---------')
    bands = {'low', 'high'};
    for iBands = 1:length(bands)
        fprintf(stats_log,'####    Circ Stats for %s gamma    ####\n', bands{iBands})
        for ii = 1:size(str_vals_piri.(bands{iBands}),1)
            for jj = 1:size(str_vals_piri.(bands{iBands}),2)
                if ~isempty(str_vals_piri.(bands{iBands}){ii, jj})
                    fprintf(stats_log,str_vals_piri.(bands{iBands}){ii, jj})
                    fprintf(stats_log,'\n')
                end
                
            end
        end
    end
    
    %% stats for latex table:
    for ii = 1:size(phase_mean,1)
        for jj = 1:size(phase_mean,2)
            if ~isempty(phase_mean(ii,jj)) && ~isnan(phase_mean(ii,jj)) &&  low_mat(ii,jj) &&isempty(strfind(labels{ii,jj}, 'Piri'))
                if strcmp(labels{ii,jj}, 'PiriN_PiriO') || strcmp(labels{ii,jj}, 'NAc_OFC')
                    fprintf(stats_log,'%s}: & %.0f & %.2f & %.0f & %.2f', ['\textbf{' strrep(labels{ii,jj}, '_', '-')], rad2deg(phase_mean(ii,jj)), vector_r(ii,jj), rad2deg(phase_mean(jj,ii)), vector_r(jj,ii))
                else
                    fprintf(stats_log,'%s: & %.0f & %.2f & %.0f & %.2f', strrep(labels{ii,jj}, '_', '-'), round(rad2deg(phase_mean(ii,jj))), vector_r(ii,jj), round(rad2deg(phase_mean(jj,ii))), vector_r(jj,ii))
                end
                fprintf(stats_log,'\\\\ \n ')
                fprintf(stats_log,'%s', '\hline')
                fprintf(stats_log,'\n')
            end
            
        end
    end
    fprintf(stats_log,'\n\n % Trans-Piriform\n\n')
    % same thing but with piriform
    for ii = 1:size(phase_mean,1)
        for jj = 1:size(phase_mean,2)
            if ~isempty(phase_mean(ii,jj)) && ~isnan(phase_mean(ii,jj)) &&  low_mat(ii,jj) && ~isempty(strfind(labels{ii,jj}, 'Piri'))
                if strcmp(labels{ii,jj}, 'PiriN_PiriO') || strcmp(labels{ii,jj}, 'NAc_OFC')
                    fprintf(stats_log,'%s}: & %.0f & %.2f & %.0f & %.2f', ['\textbf{' strrep(labels{ii,jj}, '_', '-')], rad2deg(phase_mean(ii,jj)), vector_r(ii,jj), rad2deg(phase_mean(jj,ii)), vector_r(jj,ii))
                    
                else
                    fprintf(stats_log,'%s: & %.0f & %.2f & %.0f & %.2f', strrep(labels{ii,jj}, '_', '-'), round(rad2deg(phase_mean(ii,jj))), vector_r(ii,jj), round(rad2deg(phase_mean(jj,ii))), vector_r(jj,ii))
                end
                fprintf(stats_log,'\\\\ \n ')
                fprintf(stats_log,'%s', '\hline')
                fprintf(stats_log,'\n')
            end
            
        end
    end
end

% %% make a table of the pairs and corresponding values
% %% check the rayleigh score of each and report
% R = [];
% for iPhase =1:length(PARAMS.Phases)
%     for iBand = 1:length(bands);
%         for ii = 1:size(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).(bands{iBand}),1)
%             for jj = 1:size(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).(bands{iBand}),2)
%                 if ~isempty(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
%                     R.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = circ_rtest(...
%                         all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj});
%                 else
%                     R.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} = NaN;
%                 end
%             end
%         end
%     end
% end
% figure
% subplot(1,2,1)
% h = nan_imagesc_ec(cell2mat(R.contra.low));
% add_num_imagesc(h, cell2mat(R.contra.low), 4)
% caxis([0 0.001])
% subplot(1,2,2)
% h = nan_imagesc_ec(cell2mat(R.contra.high));
% add_num_imagesc(h, cell2mat(R.contra.high), 4)
%
% caxis([0 0.001])


%% same plot as above but for same side of PC vs opposite

% for iPhase =1:length(PARAMS.Phases)
%     %         phase_mean_mat.(PARAMS.Phases{iPhase}) = NaN(length(labels), length(labels));
%     str_vals.low = cell(7,7); str_vals.high = cell(7,7);
%     str_vals_piri.low = cell(7,7); str_vals_piri.high = cell(7,7);
%     phase_mean = NaN(7,7); vector_r = NaN(7,7);
%
%     pol_legend = {}; piri_legend = {};
%     n = 1000;
%     figure(1116)
%     all_h = [];
%
%     %
%     for ii = 1:size(labels,1)
%         for jj = 1:size(labels,2)
%             iC = find(ismember(c_labels, labels{ii,jj}));
%             if ~isnan(mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj)) && low_mat(ii,jj)
%                 %                 iC = find(ismember(c_labels, labels{ii,jj}));
%
%                 phase_mean(ii,jj) = circ_mean(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).low{ii,jj});
%                 vector_r(ii,jj) = circ_r(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).low{ii,jj});
%
%                 if isempty(strfind(labels{ii,jj}, 'Piri'))
%                     colour = PARAMS.fig_pink;
%                     %                     if ~isempty(strcmp(labels{ii,jj}, 'PiriN_PiriO'))
%                     subplot(2,2,1)
%                     pol_legend = [pol_legend, labels{ii,jj}];
%                     str_vals.low{ii,jj} = sprintf('%s: & %.0f & %.2f',labels{ii,jj}, rad2deg(phase_mean(ii,jj)), vector_r(ii,jj)); % with values.
%                     %                     end
%                 elseif strcmp(labels{ii,jj}, 'PiriN_PiriO') || strcmp(labels{ii,jj}, 'PiriO_PiriN')
%                     %                     colour = PARAMS.fig_green;
%                     colour = [0.6 0.6 0.6];
%                     subplot(2,2,1)
%                     piri_legend = [piri_legend, labels{ii,jj}];
%                     str_vals_piri.low{ii,jj} = sprintf('%s: & %.0f & %.2f',labels{ii,jj}, rad2deg(phase_mean(ii,jj)), vector_r(ii,jj)); % with values.
%                 else
%                     colour = PARAMS.fig_green;
%                     subplot(2,2,1)
%                     piri_legend = [piri_legend, labels{ii,jj}];
%                     str_vals_piri.low{ii,jj} = sprintf('%s: & %.0f & %.2f',labels{ii,jj}, rad2deg(phase_mean(ii,jj)), vector_r(ii,jj)); % with values.
%
%                 end
%
%                 if ver >= 2017
%                     h = polarplot([0 phase_mean(ii,jj)],[0 vector_r(ii,jj)]);
%                     rlim([0 1])
%                     set(gca, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top','RTick' ,[0.5 1], 'FontSize', 22, 'FontName', 'Helvetica')
%                     cfg_polar.axes_labels = 0; % turns off any scaling of the x and y label which don't work in polar plots
%                 else
%                     h_temp = polar([0 0], [0 1]);
%                     hold on
%                     h = polar([0 phase_mean(ii,jj)], [0 vector_r(ii,jj)]);
%                     set(h_temp, 'visible', 'off') ;
%                     cfg_polar.axes_labels = 1; % turns off any scaling of the x and y label which don't work in polar plots
%                 end
%
%                 S = strsplit(labels{ii,jj}, '_');
%                 %                 if strcmp(S{1}, 'NAc') && strcmp(S{2}, 'OFC') || strcmp(S{1}, 'PiriN') && strcmp(S{2}, 'PiriO')
%                 %                     set(h, 'color', 'k', 'linewidth', 7)
%                 %                 else
%                 set(h, 'color', colour, 'linewidth', 4)
%                 %                 end
%                 hold all
%                 %                     view([90, 90])
%
%             elseif ~isnan(mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj)) && high_mat(ii,jj) && ~isempty(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj})
%                 iC = find(ismember(c_labels, labels{jj,ii}));
%
%                 phase_mean(ii,jj) = circ_mean(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj});
%                 vector_r(ii,jj) = circ_r(all_mean.Phase_diff_all.(PARAMS.Phases{iPhase}).high{ii,jj});
%
%                 if isempty(strfind(labels{ii,jj}, 'Piri')) % contains no piriform channel
%                     %                     if ~isempty(strcmp(labels{ii,jj}, 'PiriN_PiriO'))
%                     colour = PARAMS.fig_blue;
%                     subplot(2,2,3)
%                     str_vals.high{ii,jj} = sprintf('%s: & %.0f  & %.2f',labels{jj,ii}, rad2deg(phase_mean(ii,jj)), vector_r(jj,ii)); % with values.
%                     %                     end
%                 elseif strcmp(labels{ii,jj}, 'PiriN_PiriO') || strcmp(labels{ii,jj}, 'PiriO_PiriN')
%                     colour = [0.6 0.6 0.6];
%                     subplot(2,2,3)
%                     str_vals_piri.high{ii,jj} = sprintf('%s: & %.0f  & %.2f',labels{jj,ii}, rad2deg(phase_mean(ii,jj)), vector_r(jj,ii)); % with values.
%                 else % if it is a piriform - x comparison
%                     colour = PARAMS.fig_green;
%                     subplot(2,2,3)
%                     str_vals_piri.high{ii,jj} = sprintf('%s: & %.0f  & %.2f',labels{jj,ii}, rad2deg(phase_mean(ii,jj)), vector_r(jj,ii)); % with values.
%                 end
%
%                 if ver >= 2017
%                     h = polarplot( [0 phase_mean(ii,jj)],[0 vector_r(ii,jj)]);
%                     rlim([0 1])
%                     set(gca, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top','RTick' ,[0.5 1], 'FontSize', 22, 'FontName', 'Helvetica')
%                     cfg_polar.axes_labels = 0; % turns off any scaling of the x and y label which don't work in polar plots
%                 else
%                     h_temp = polar([0 0], [0 1]);
%                     hold on
%                     h = polar([0 phase_mean(ii,jj)], [0 vector_r(ii,jj)]);
%                     set(h_temp, 'visible', 'off') ;
%                     cfg_polar.axes_labels = 1; % turns off any scaling of the x and y label which don't work in polar plots
%                 end
%
%
%                 %                 S = strsplit(labels{ii,jj}, '_');
%                 %                 if strcmp(S{2}, 'NAc') && strcmp(S{1}, 'OFC') || strcmp(S{1}, 'PiriN') && strcmp(S{2}, 'PiriO')
%                 %                     set(h, 'color', 'k', 'linewidth', 7)
%                 %                 else
%                 set(h, 'color', colour, 'linewidth', 4)
%                 %                 end
%                 hold all
%                 %                     view([270, 90])
%                 %                     view([0, 90])
%             end
%         end
%     end
%     %     text(.5,1.15,'low','units','normalized','fontsize', 24,'horizontalalignment','center','verticalalignment','top');
%     %         text(2,1.15,'high','units','normalized','fontsize', 24,'horizontalalignment','center','verticalalignment','top');
%
%     th = findall(gcf,'Type','text');
%     for i = 1:length(th)
%         set(th(i),'FontSize',20)
%     end
%     subplot(2,2,1)
%     hL = legend(pol_legend);
%     set(hL,'Position', [0.5 0.7 0.05 0.1], 'Units', 'normalized')
%
%     subplot(2,2,3)
%     hL = legend(piri_legend);
%     set(hL,'Position', [0.5 0.25 0.05 0.1], 'Units', 'normalized')
%     %         tightfig
%     cfg_polar.resize = 1;
%     cfg_polar.pos = [205    50  1200   800];
%     SetFigure(cfg_polar, gcf);
%
%     saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' (PARAMS.Phases{iPhase}) '_polar_phase_across.fig']);
%     %         saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' (PARAMS.Phases{iPhase}) '_polar_phase'], 'epsc')
%     saveas_eps([iSub{1} '_' (PARAMS.Phases{iPhase}) '_polar_phase_across'], [PARAMS.inter_dir 'Phase_plots/'])
%     print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' (PARAMS.Phases{iPhase})  '_polar_phase_across.png'], '-r300')
%     close all
% end
%% try to make the compass plots for the phase differences
% conver the tril or triu to a linear
%     for iPhase =1:length(PARAMS.Phases)
%         for iBand = 1:2
%             Polar_phase.(PARAMS.Phases{iPhase}).(bands{iBand}) = [];
%         end
%     end
%
%     pol_legend = {};
%     for ii = 1:size(labels,1)
%         for jj = 1:size(labels,2)
%             if low_mat(ii,jj)
%                 Polar_phase.(PARAMS.Phases{iPhase}).low =[Polar_phase.(PARAMS.Phases{iPhase}).low, mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj)];
%                 pol_legend = [pol_legend, labels{ii,jj}];
%             elseif ~low_mat(ii,jj) && ii == jj
%                 Polar_phase.(PARAMS.Phases{iPhase}).high = [Polar_phase.(PARAMS.Phases{iPhase}).high, mean_out.Phase_diff_out.(PARAMS.Phases{iPhase})(ii,jj)];
%
%             end
%         end
%     end


%%

%
%     iC = find(ismember(c_labels, labels{ii,jj}));
%     [t,r]= rose(Polar_phase.(PARAMS.Phases{iPhase}).low);
%     h = compass(t,r);
%     polarscatter(t,r,100*(ones(1,length(r))),1:length(r),'filled','MarkerFaceAlpha',.5)
%
%     legend(pol_legend)
%     for ii=1:length(h)
%         set(h(ii), 'color',c_ord(iC,:), 'linewidth', 2)
%     end
%     legend(pol_legend)
%
%
%     SetFigure([], gcf)
%


%%

%     %% skittles version
%     skittles = linspecer(length(labels)-1);
%     for iMs = 1:length(measures)
%         % get the line plots for certain veriables
%         if strcmp(measures{iMs}, 'Phase_lag_cxy') || strcmp(measures{iMs}, 'PS_slope') || strcmp(measures{iMs}, 'COH_cxx')
%             Phase_c = linspecer(4);
%             low_mat = tril(ones(size(labels)),-1);
%             for iPhase =1:length(PARAMS.Phases)
%                 figure(iPhase);
%
%                 % plot the coherence
%                 labs = {};
%                 handles_p = {};
%                 for ii = 1:size(labels,1)
%                     subplot(1,7,ii)
%                     for jj = 1:size(labels,2)
%                         if ~isempty(all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).low{ii, jj})
%                             iC = find(ismember(c_labels, labels{ii,jj}));
%                             hold on
%                             % for the coherence
%                             if strcmp(measures{iMs}, 'COH_cxx')
%                                 h(iC) =shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).low{ii, jj}, all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).low{ii, jj}/sqrt(length(all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).low{ii, jj})));
%                                 h(iC).mainLine.Color = c_ord(iC,:);
%                                 h(iC).patch.FaceColor = c_ord(iC,:);
%                                 h(iC).patch.FaceAlpha = .5;
%                                 handles_p = [handles_p, h(iC).patch];
%                             end
%                             % for the phase lag
%                             if strcmp(measures{iMs}, 'Phase_lag_cxy')
%                                 sem =  rad2deg(all_std.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj});%/sqrt(length(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}));
%                                 hE(iC) = errorbar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), sem, 'o', 'color', c_ord(iC,:));
%                                 handles_p = [handles_p, hE(iC)];
%                             end
%                             % for the phase slope
%                             if strcmp(measures{iMs}, 'PS_slope')
%                                 h(iC) =shadedErrorBar(all_mean.PS_F.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).low{ii, jj}, all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).low{ii, jj}/sqrt(length(all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).low{ii, jj})));
%                                 h(iC).mainLine.Color = c_ord(iC,:);
%                                 h(iC).patch.FaceColor = c_ord(iC,:);
%                                 h(iC).patch.FaceAlpha = .5;
%                                 handles_p = [handles_p, h(iC).patch];
%                             end
%                             %                             hV = vline([cfg.filter1.f,cfg.filter2.f], {'--b', '--b', '--g', '--g'});
%                             %                             for ihV = 1:length(hV)
%                             %                                 hV(ihV).LineWidth= 2;
%                             %                             end
%
%                             labs = [labs, labels{ii, jj}];
%                             iC = iC+1;
%                         end
%                     end
%                 end
%                 if strcmp(measures{iMs}, 'PS_slope')
%                     ylim([-5 5])
%                     ylabel('Phase_slope (deg/Hz)');
%                     xlim([30 100])
%                 elseif strcmp(measures{iMs}, 'Phase_lag_cxy')
%                     ylim([-180 180])
%                     ylabel('Phase lag (deg)');
%                     xlim([0 100])
%                 else
%                     ylim([0 1])
%                     ylabel('Coherence');
%                     xlim([0 100])
%                 end
%                 if iPhase ==4
%                     hL = legend(handles_p, labs);
%                     set(hL,'Position', [0.45 0.5 0.05 0.1], 'Units', 'normalized')
%                 end
%                 text(.5,1.02,PARAMS.Phases{iPhase},'units','normalized','fontsize', 32,'horizontalalignment','center','verticalalignment','top', 'color', Phase_c(iPhase,:));
%                 xlabel('Frequency (Hz)');
%                 % put boxes for the gamma bands
%                 rectangle('position', [cfg.filter1.f(1), 0, (cfg.filter1.f(2)-cfg.filter1.f(1)), .05], 'FaceColor',[0.5 0.5 0.5], 'edgecolor', [0.5 0.5 0.5])
%                 rectangle('position', [cfg.filter2.f(1), 0, (cfg.filter2.f(2)-cfg.filter2.f(1)), .05], 'FaceColor',[0.7 0.7 0.7], 'edgecolor', [0.7 0.7 0.7])
%                 rectangle('position', [0, 0, 100, 1], 'EdgeColor', Phase_c(iPhase,:), 'linewidth', 2)
%             end
%             Square_subplots
%             tightfig
%             cfg_count_fig.resize = 1;
%             cfg_count_fig.pos = [100 40 1190 765];
%             SetFigure(cfg_count_fig, gcf)
%             saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '.fig']);
%             saveas(gcf, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs}], 'epsc')
%             print(gcf, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_' measures{iMs} '.png'])
%             close all
%         end
%     end
%

%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% plot ethe coherence
%     Phase_c = linspecer(4);
%     low_mat = tril(ones(size(labels)),-1);
%     coh_f =figure(1111);
%     % for iBand = 1:2
%     for iPhase =1:length(PARAMS.Phases)
%         % plot the coherence
%         subplot(2,4,iPhase)
%         cxx_labs = {};
%         handles_p = {};
%         for ii = 1:size(labels,1)
%             for jj = 1:size(labels,2)
%                 if ~isempty(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj}) && low_mat(ii, jj)
%                     iC = find(ismember(c_labels, labels{ii,jj}));
%                     hold on
%                     %                 plot(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj})
%                     h(iC) =shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_std.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj}/sqrt(length(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj})));
%                     %                 shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_std.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}/sqrt(length(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj})),'b')
%                     h(iC).mainLine.Color = c_ord(iC,:);
%                     h(iC).patch.FaceColor = c_ord(iC,:);
%                     h(iC).patch.FaceAlpha = .5;
%                     % vline(cfg_filter.f)
%                     %                 vline([cfg_filter.f,cfg_filter2.f], {'--b', '--b', '--g', '--g'})
%                     % p_95 = prctile(shuf_coh,95, 1);
%                     % plot(shift_fxx, p_95, '--r');
%                     % plot(shift_fxx, p_95, '--r')
%                     handles_p = [handles_p, h(iC).patch];
%                     cxx_labs = [cxx_labs, labels{ii, jj}];
%                     iC = iC+1;
%                 end
%             end
%         end
%         xlim([0 100])
%         ylim([0.1 1])
%         legend(handles_p, cxx_labs)
%         text(.5,1,PARAMS.Phases{iPhase},'units','normalized','horizontalalignment','center','verticalalignment','top', 'color', Phase_c(iPhase,:));
%         xlabel('Frequency (Hz)');
%         ylabel('Coherence');
%
%         %% add the image SC to the bottom
%         %     for iPhase =1:length(PARAMS.Phases)
%
%         subplot(2,4,4+iPhase)
%         nan_imagesc_ec(Coh_out.(PARAMS.Phases{iPhase}))
%         [~, hStrings] =add_num_imagesc(gca,Coh_out.(PARAMS.Phases{iPhase}) , 2,  10); % adds numerics to imagesc
%         %     hcb = colorbar();
%         colormap('parula')
%         % caxis([0 .5]);
%         %     set(hcb, 'Ytick', 0:0.1:.5)
%         set(gca,'xtick', 1:length(single_labels), 'ytick', 1:length(single_labels),  'xticklabel', single_labels, 'yticklabel', single_labels)
%         xtickangle(65)
%     end
%     SetFigure([], gcf)
%     tightfig
%     saveas(coh_f, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_coh.fig']);
%     saveas(coh_f, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_coh'], 'epsc')
%     print(coh_f, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_coh.png'])
%     close all
%     %% plot the phase offset
%     c_ord = linspecer(24);
%     Phase_c = linspecer(4);
%     low_mat = tril(ones(size(labels)),-1);
%     Phase_off_f = figure(1112);
%     % for iBand = 1:2
%     for iPhase =1:length(PARAMS.Phases)
%         % plot the coherence
%         subplot(3,4,iPhase)
%         cxx_labs = {};
%         iC = 1; handles_p = {};
%         for ii = 1:size(labels,1)
%             for jj = 1:size(labels,2)
%                 if ~isempty(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}) &&  low_mat(ii, jj);
%                     hold on
%                     sem =  rad2deg(all_std.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj});%/sqrt(length(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}));
%                     errorbar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), sem, 'o', 'color', c_ord(iC,:))
%                     %                 plot(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), 'o', 'color', c_ord(iC,:))
%                     %                 plot(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj})
%                     %                 h(iC) =shadedErrorBar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), rad2deg(all_std.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj})/sqrt(length(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj})));
%                     %                 shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_std.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}/sqrt(length(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj})),'b')
%                     %                 h(iC).mainLine.Color = c_ord(iC,:);
%                     %                 h(iC).patch.FaceColor = c_ord(iC,:);
%                     %                 h(iC).patch.FaceAlpha = .5;
%                     % vline(cfg_filter.f)
%                     %                 vline([cfg_filter.f,cfg_filter2.f], {'--b', '--b', '--g', '--g'})
%                     % p_95 = prctile(shuf_coh,95, 1);
%                     % plot(shift_fxx, p_95, '--r');
%                     % plot(shift_fxx, p_95, '--r')
%                     %                 handles_p = [handles_p, h(iC).patch];
%                     cxx_labs = [cxx_labs, labels{ii, jj}];
%                     iC = iC+1;
%                 end
%             end
%         end
%         xlim([0 100])
%         %     ylim([0.1 1])
%         legend(handles_p, cxx_labs)
%         text(.5,1,PARAMS.Phases{iPhase},'units','normalized','horizontalalignment','center','verticalalignment','top', 'color', Phase_c(iPhase,:));
%         xlabel('Frequency (Hz)');
%         ylabel('Phase offset');
%         %% make a polar plot
%         %     subplot(3,4,4+iPhase)
%         %     n = 100;
%         %     [t,r] = rose(Phase_diff_out.(PARAMS.Phases{iPhase}){ii, jj},n)
%
%         %% add the image SC to the bottom
%         %     for iPhase =1:length(PARAMS.Phases)
%
%         subplot(3,4,8+iPhase)
%         nan_imagesc_ec(rad2deg(Phase_diff_out.(PARAMS.Phases{iPhase})))
%         [~, hStrings] =add_num_imagesc(gca,rad2deg(Phase_diff_out.(PARAMS.Phases{iPhase})), 2,  10); % adds numerics to imagesc
%         %     hcb = colorbar();
%         colormap('parula')
%         % caxis([0 .5]);
%         %     set(hcb, 'Ytick', 0:0.1:.5)
%         set(gca,'xtick', 1:length(single_labels), 'ytick', 1:length(single_labels),  'xticklabel', single_labels, 'yticklabel', single_labels)
%         xtickangle(65)
%     end
%
%     SetFigure([], gcf)
%     tightfig
%     saveas(Phase_off_f, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_Phase_off.fig']);
%     saveas(Phase_off_f, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_Phase_off'], 'epsc')
%     print(Phase_off_f, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_Phase_off.png'])
%
%     close all
%     %% same for the number of events
%     c_ord = linspecer(24);
%     Phase_c = linspecer(4);
%     low_mat = tril(ones(size(labels)),-1);
%     Evt_count_f = figure(1113);
%     for iPhase =1:length(PARAMS.Phases)
%
%         subplot(1,4,+iPhase)
%         nan_imagesc_ec(Evt_count_out.(PARAMS.Phases{iPhase}))
%         [~, hStrings] =add_num_imagesc(gca,Evt_count_out.(PARAMS.Phases{iPhase}), 0,  10); % adds numerics to imagesc
%         %     hcb = colorbar();
%         colormap('parula')
%         % caxis([0 .5]);
%         %     set(hcb, 'Ytick', 0:0.1:.5)
%         set(gca,'xtick', 1:length(single_labels), 'ytick', 1:length(single_labels),  'xticklabel', single_labels, 'yticklabel', single_labels)
%         xtickangle(65)
%     end
%%
%     Square_subplots
%     tightfig
%     cfg_count_fig.resize = 1;
%     cfg_count_fig.pos = [3 451 1438 333];
%     SetFigure(cfg_count_fig, gcf)
%
%     saveas(Evt_count_f, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_count.fig']);
%     saveas(Evt_count_f, [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_count'], 'epsc')
%     print(Evt_count_f, '-dpng', [PARAMS.inter_dir 'Phase_plots/' iSub{1} '_count.png'])
%     %%
close all
% clear everything except the subject
clearvars -except Subjects iSub PARAMS
end
