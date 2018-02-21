function MS_plot_phase(cfg_in, mat_in)
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

%% temp
sess_list = fieldnames(mat_out.R104);
PS = {}
for iSess = 1:length(fieldnames(mat_out.R104))
    PS = cat(3,PS, mat_out.R104.(sess_list{iSess}).contra.low.PS)
end
    
    
    
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
