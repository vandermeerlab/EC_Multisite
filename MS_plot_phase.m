%% for debugging
run('/Users/jericcarmichael/Documents/GitHub/EC_Multisite/MS_initialize.m')

load([PARAMS.inter_dir 'R107_phase_out.mat']);
mat_in = Phase_mat;
clear Phase_mat;
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
labels = mat_in.(sess_list{1}).labels;
bands = {'low', 'high'};
all_mean = []; all_std = [];
measures = fieldnames(mat_in.(sess_list{1}).(PARAMS.Phases{1}));
%% move everything to the all_out structure
sess_list = fieldnames(mat_in);

for iPhase =2%:length(PARAMS.Phases)
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
                        if size(temp{1,1,1},1) > size(temp{1,1,1},2)
                            temp(cellfun(@(temp) any(isnan(temp)),temp)) = [];
                            temp = cell2mat(squeeze(temp)');
                        else
                            temp(cellfun(@(temp) any(isnan(temp)),temp)) = [];
                            temp = cell2mat(squeeze(temp));
                        end
                        if strcmp(measures{iMs}, 'Phase_lag_cxy')  || strcmp(measures{iMs}, 'Phase_lag_mean') || strcmp(measures{iMs}, 'PS_slope')
                            all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj}=  circ_median(temp');
                            all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  circ_std(temp');
                            
                        elseif strcmp(measures{iMs}, 'Phase_lag_F')  || strcmp(measures{iMs}, 'COH_fxx') || strcmp(measures{iMs}, 'PS_D')
                            all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanmean(temp, 2);
                            
                        else
                            all_mean.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanmedian(temp,2);
                            all_std.(measures{iMs}).(PARAMS.Phases{iPhase}).(bands{iBand}){ii,jj} =  nanstd(temp,[],2);
                        end
                    end
                    clear temp
                end
            end
        end
    end
end


%% plot everything
c_ord = linspecer(7);
Phase_c = linspecer(4);
low_mat = tril(ones(size(labels)),-1);
figure(1111)
% for iBand = 1:2
for iPhase =2%:length(PARAMS.Phases)
    % plot the coherence
    subplot(2,4,iPhase)
    cxx_labs = {};
    iC = 1; handles_p = {};
    for ii = 1:size(labels,1)
        for jj = 1:size(labels,2)
            if ~isempty(all_mean.COH_cxx.(PARAMS.Phases{iPhase}).low{ii, jj}) && low_mat(ii, jj);
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
end
%%
nan_imagesc_ec(mean_mat)
legend
[h, hStrings] =add_num_imagesc(gca,mean_mat*100 , 2,  18); % adds numerics to imagesc
hcb = colorbar();
colormap('parula')
% caxis([0 .5]);
set(hcb, 'Ytick', 0:0.1:.5)
set(gca,'xtick', 1:4, 'ytick', 1:4,  'xticklabel', labels, 'yticklabel', labels)
SetFigure([], gcf)
%     end
% end
%