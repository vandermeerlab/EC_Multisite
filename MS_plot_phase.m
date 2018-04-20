%% for debugging
MS_initialize

mat_in = load([PARAMS.inter_dir 'R107_phase_out.mat']);
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
%
%
%% pre allocate all the sufields of the different analyses;

sess_list = fieldnames(mat_in);
bands = {'low', 'high'};

for iPhase =1:length(PARAMS.Phases)
    measures = fieldnames(mat_in.(sess_list{1}).(PARAMS.Phases{iPhase}));
    for iMs = 1:length(measures)
        for iBand = 1:2
            all_out.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}) = {};
        end
    end
end
%% temp
sess_list = fieldnames(mat_in);

for iSess = [1, 3]%1:length(sess_list)
    for iPhase =1:length(PARAMS.Phases)
        for iMs = 1:length(measures)
            for iBand = 1:2
                all_out.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}) = cat(3,all_out.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}), mat_in.(sess_list{iSess}).(PARAMS.Phases{iPhase}).(measures{iMs}).(bands{iBand}));
            end
        end
    end
end


%%
labels = {'PL', 'OFC', 'NAc', 'CG'};
for ii =1:size(PS,1)
    for jj =1:size(PS,2)
        mean_mat(ii, jj) = mean(cell2mat([PS{ii, jj,:}]));
        S_mat(ii, jj) = length(cell2mat([PS{ii, jj,:}]));
    end
end

%% plot
nan_imagesc_ec(mean_mat)
legend
[h, hStrings] =add_num_imagesc(gca,mean_mat*100 , 2,  18); % adds numerics to imagesc
hcb = colorbar();
colormap('parula')
% caxis([0 .5]);
set(hcb, 'Ytick', 0:0.1:.5)
set(gca,'xtick', 1:4, 'ytick', 1:4,  'xticklabel', labels, 'yticklabel', labels)
SetFigure([], gcf)
