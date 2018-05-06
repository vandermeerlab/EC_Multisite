%% for debugging
run('/Users/jericcarmichael/Documents/GitHub/EC_Multisite/MS_initialize.m')

load([PARAMS.inter_dir '/Phase_outputs/R108_phase_out.mat']);
mat_in = Phase_mat;
clear Phase_mat;
cfg_in = [];
%%
% function MS_plot_phase(cfg_in, mat_in)
%% MS_plot_phase: plots the outputs from the three phase measures in the
%       MS project: Phase correlation, amplitude cross-corelation, phase
%       slope index.  mat_in is an input taken from the
%       MS_get_phase_metrics.
%
%   Inputs:
%       - cfg_in [strcut]: contains configuration parameters.  Inludes:
%                  - cfg_in.type : what type of phase measure to plot. can
%                  be 'phase_coh", "amp_xcorr", 'psi'.
%
%       - mat_in [struct]:
%
%
%% default configurations
cfg_def = [];
cfg_def.filter = [];
cfg_def.filter1.f = [45 65];
cfg_def.filter2.f = [70 90];

cfg = ProcessConfig2(cfg_def, cfg_in);



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
    for iMs = 1:length(measures)
        for iBand = 1:2
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
                    if isempty(temp{1,1,1})
                        all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  [];
                        all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  [];
                    else
                        %                         if size(temp{1,1,1},1) > size(temp{1,1,1},2)
                        %                             temp(cellfun(@(temp) any(isnan(temp)),temp)) = [];
                        %                             temp = cell2mat(squeeze(temp)');
                        %                         else
                        temp(cellfun(@(temp) any(isnan(temp)),temp)) = [];
                        temp = cell2mat(squeeze(temp));
                        %                         end
                        if strcmp(measures{iMs}, 'Phase_lag_cxy')  || strcmp(measures{iMs}, 'Phase_lag_mean') || strcmp(measures{iMs}, 'PS_slope')
                            all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj}=  circ_mean(temp);
                            all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  circ_std(temp);
                            
                        elseif strcmp(measures{iMs}, 'Phase_lag_F')  || strcmp(measures{iMs}, 'COH_fxx') || strcmp(measures{iMs}, 'PS_D')
                            all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanmean(temp, 1);
                            
                        else
                            all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanmedian(temp,1);
                            all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanstd(temp,[],1);
                        end
                    end
                    clear temp
                end
            end
        end
    end
end
%% get the average within the frequency bands
for iPhase = 1:length(PARAMS.Phases)
    for iBand = 1:length(bands)
        for ii = size(labels,1):-1:1
            for jj = size(labels,2):-1:1
                if isempty(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    Coh_mean.(bands{iBand})(ii, jj) = NaN;
                    % if it is in the lower triangular put the low
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    this_cxx = all_mean.COH_cxx.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    this_f = all_mean.COH_fxx.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj};
                    Coh_mean.(bands{iBand})(ii, jj) = nanmean(this_cxx(nearest_idx(cfg.(['filter' num2str(iBand)]).f(1), this_f):nearest_idx(cfg.(['filter' num2str(iBand)]).f(2), this_f)));
                end
                
                % same for the phase lag
                if isempty(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    Phase_diff.(bands{iBand})(ii, jj) = NaN;
                    % if it is in the lower triangular put the low
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    Phase_diff.(bands{iBand})(ii, jj) = circ_mean(all_mean.Phase_lag_mean.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj});
                end
                
                % for the phase slope
                if isempty(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj})
                    Phase_diff.(bands{iBand})(ii, jj) = NaN;
                    % if it is in the lower triangular put the low
                elseif ismember(Idx_mat(ii, jj), tril(Idx_mat,-1))
                    Phase_diff.(bands{iBand})(ii, jj) = circ_mean(all_mean.Phase_lag_mean.(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj});
                end
            end
        end
    end
    Coh_out.(PARAMS.Phases{iPhase}) = tril(Coh_mean.low) + tril(Coh_mean.high,-1)';
    Phase_diff_out.(PARAMS.Phases{iPhase}) = tril(Phase_diff.low) + tril(Phase_diff.high,-1)';

end
%% plot ethe coherence
c_ord = linspecer(44);
Phase_c = linspecer(4);
low_mat = tril(ones(size(labels)),-1);
coh_f =figure(1111);
% for iBand = 1:2
for iPhase =1:length(PARAMS.Phases)
    % plot the coherence
    subplot(2,4,iPhase)
    cxx_labs = {};
    iC = 1; handles_p = {}; 
    for ii = 1:size(labels,1)
        for jj = 1:size(labels,2)
            if ~isempty(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj});
                hold on
                %                 plot(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj})
                h(iC) =shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_std.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj}/sqrt(length(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj})));
                %                 shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_std.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}/sqrt(length(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj})),'b')
                h(iC).mainLine.Color = c_ord(iC,:);
                h(iC).patch.FaceColor = c_ord(iC,:);
                h(iC).patch.FaceAlpha = .5;
                % vline(cfg_filter.f)
                %                 vline([cfg_filter.f,cfg_filter2.f], {'--b', '--b', '--g', '--g'})
                % p_95 = prctile(shuf_coh,95, 1);
                % plot(shift_fxx, p_95, '--r');
                % plot(shift_fxx, p_95, '--r')
                handles_p = [handles_p, h(iC).patch];
                cxx_labs = [cxx_labs, labels{ii, jj}];
                iC = iC+1;
            end
        end
    end
    xlim([0 100])
    ylim([0.1 1])
    legend(handles_p, cxx_labs)
    text(.5,1,PARAMS.Phases{iPhase},'units','normalized','horizontalalignment','center','verticalalignment','top', 'color', Phase_c(iPhase,:));
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    
    %% add the image SC to the bottom
%     for iPhase =1:length(PARAMS.Phases)

    subplot(2,4,4+iPhase)
    nan_imagesc_ec(Coh_out.(PARAMS.Phases{iPhase}))
    [~, hStrings] =add_num_imagesc(gca,Coh_out.(PARAMS.Phases{iPhase}) , 2,  10); % adds numerics to imagesc
%     hcb = colorbar();
    colormap('parula')
    % caxis([0 .5]);
%     set(hcb, 'Ytick', 0:0.1:.5)
    set(gca,'xtick', 1:length(single_labels), 'ytick', 1:length(single_labels),  'xticklabel', single_labels, 'yticklabel', single_labels)
end
tightfig
% saveas(coh_f, ['

%% plot the phase offset
c_ord = linspecer(24);
Phase_c = linspecer(4);
low_mat = tril(ones(size(labels)),-1);
Phase_off_f = figure(1112);
% for iBand = 1:2
for iPhase =1:length(PARAMS.Phases)
    % plot the coherence
    subplot(3,4,iPhase)
    cxx_labs = {};
    iC = 1; handles_p = {}; 
    for ii = 1:size(labels,1)
        for jj = 1:size(labels,2)
            if ~isempty(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}) &&  low_mat(ii, jj);
                hold on
                sem =  rad2deg(all_std.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj});%/sqrt(length(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}));
                errorbar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), sem, 'o', 'color', c_ord(iC,:))
%                 plot(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), 'o', 'color', c_ord(iC,:))
                %                 plot(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).low{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj})
%                 h(iC) =shadedErrorBar(all_mean.Phase_lag_F.(PARAMS.Phases{iPhase}).low{ii, jj}, rad2deg(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj}), rad2deg(all_std.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj})/sqrt(length(all_mean.Phase_lag_cxy.(PARAMS.Phases{iPhase}).low{ii, jj})));
                %                 shadedErrorBar(all_mean.COH_fxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}, all_std.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj}/sqrt(length(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).high{ii, jj})),'b')
%                 h(iC).mainLine.Color = c_ord(iC,:);
%                 h(iC).patch.FaceColor = c_ord(iC,:);
%                 h(iC).patch.FaceAlpha = .5;
                % vline(cfg_filter.f)
                %                 vline([cfg_filter.f,cfg_filter2.f], {'--b', '--b', '--g', '--g'})
                % p_95 = prctile(shuf_coh,95, 1);
                % plot(shift_fxx, p_95, '--r');
                % plot(shift_fxx, p_95, '--r')
%                 handles_p = [handles_p, h(iC).patch];
                cxx_labs = [cxx_labs, labels{ii, jj}];
                iC = iC+1;
            end
        end
    end
    xlim([0 100])
%     ylim([0.1 1])
    legend(handles_p, cxx_labs)
    text(.5,1,PARAMS.Phases{iPhase},'units','normalized','horizontalalignment','center','verticalalignment','top', 'color', Phase_c(iPhase,:));
    xlabel('Frequency (Hz)');
    ylabel('Phase offset');
    %% make a polar plot
%     subplot(3,4,4+iPhase)
%     n = 100;
%     [t,r] = rose(Phase_diff_out.(PARAMS.Phases{iPhase}){ii, jj},n)
    
    %% add the image SC to the bottom
%     for iPhase =1:length(PARAMS.Phases)

    subplot(3,4,8+iPhase)
    nan_imagesc_ec(rad2deg(Phase_diff_out.(PARAMS.Phases{iPhase})))
    [~, hStrings] =add_num_imagesc(gca,rad2deg(Phase_diff_out.(PARAMS.Phases{iPhase})), 2,  10); % adds numerics to imagesc
%     hcb = colorbar();
    colormap('parula')
    % caxis([0 .5]);
%     set(hcb, 'Ytick', 0:0.1:.5)
    set(gca,'xtick', 1:length(single_labels), 'ytick', 1:length(single_labels),  'xticklabel', single_labels, 'yticklabel', single_labels)
end