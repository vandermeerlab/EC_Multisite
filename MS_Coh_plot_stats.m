function [stats_out] = MS_Coh_plot_stats(Coh_mat)
%% MS_Coh_plot_stats: takes the matrix containing the coherence vaules from
%  MS_event_pairs and averages the matrices across sessions, plots the heat
%  maps, and then runs very basic stats.
%
%
%  inputs:
%     - Coh_mat: [struct] contains the matrices for all the mean coherence
%     values for each pair of electrodes in each subject.
%
%
%  outputs:
%     - stats_out: [struct] contains the statisics for the coherence values
%
%
global PARAMS
%% Covert the matrices to a 3d matrix for each
phases = PARAMS.Phases;
bands = {'low', 'high'};
types = {'evt', 'sess_amp', 'sess'};
%initialize the all_mats
for iPhase = 1:length(phases)
    for iBand = 1:length(bands)
        for iType = 1:length(types)
        All_mat.(phases{iPhase}).(bands{iBand}).(types{iType}) = [];
        end
    end
end
% collect all the data
subs = fieldnames(Coh_mat);
for iSub = 1:length(subs)
    for iPhase = 1:length(phases)
        if strcmp(phases{iPhase}, 'labels')
            continue
        else
            for iBand = 1:2
                for iType = 1:length(types)
                    % collect the evt coherence matrices
                    temp.(phases{iPhase}).(bands{iBand}).(types{iType})  = cat(3, Coh_mat.(subs{iSub}).(sess_list{1}).(phases{iPhase}).(types{iType}).(bands{iBand}),...
                        Coh_mat.(subs{iSub}).(sess_list{2}).(phases{iPhase}).(types{iType}).(bands{iBand}),...
                        Coh_mat.(subs{iSub}).(sess_list{3}).(phases{iPhase}).(types{iType}).(bands{iBand}),...
                        Coh_mat.(subs{iSub}).(sess_list{4}).(phases{iPhase}).(types{iType}).(bands{iBand}));
                    % keep the individual subjects as well if they are needed.
                    Subject_out.(subs{iSub}).(types{iType}) = temp;
                    All_mat.(phases{iPhase}).(bands{iBand}).(types{iType}) = cat(3,All_mat.(phases{iPhase}).(bands{iBand}).(types{iType}), temp.(phases{iPhase}).(bands{iBand}).(types{iType}));
                end
            end
        end
    end
end


%% get some stats



%% make some plot
labels = {};
for ii = 1:length(Coh_mat.R102.R102_2016_09_23.labels)
    t_label = strsplit(Coh_mat.R102.R102_2016_09_23.labels{ii,1}, '_');
    labels{ii} = t_label{1};
end
    
for iType = 1:length(types)
    figure(iType)
    
    for iPhase = 1:4
    plot_mat = tril(nanmean(All_mat.(phases{iPhase}).low.(types{iType}),3),-1) +(tril(nanmean(All_mat.(phases{iPhase}).high.(types{iType}),3),-1)');
    s=size(plot_mat,1);
    plot_mat(1:s+1:s*s) = NaN;
    
    subplot(1,4,iPhase)
    h = nan_imagesc_ec(plot_mat);
    add_num_imagesc(h, plot_mat)
    caxis([0 0.8])
    Square_subplots
    set(gca, 'xticklabel', labels,'ytick', 1:length(labels), 'yticklabel',(labels), 'xaxisLocation','top');
    title([phases{iPhase} '_' types{iType}])
    end
    
    
end

