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
%% Covert the matrices to a 3d matrix for each
all_mat = [];
subs = fieldnames(Coh_mat)
for iSub = 1:length(subs)
    sess_list = fieldnames(Coh_mat.(subs{iSub}))
    for iSess = 1:length(sess_list)
        phases = PARAMS.Phases;
        for iPhase = 1:length(phases)
            if strcmp(phases{iPhase}, 'labels')
                continue
            else
                
                
                
                
                
            end
        end
    end
end

