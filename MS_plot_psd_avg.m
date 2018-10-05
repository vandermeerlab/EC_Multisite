function MS_plot_psd_avg(cfg_in, Naris_in)
%% MS_plot__psd_avg: plots average power spectral densities for the data files
% in the "Naris" structure (output from MS_collect_psd)
%
% inputs:
%    -cfg_in: [struct] contains configuration paramters
%         - cfg_in.method = 'mean' or 'median'
%
%
%    -Naris: [struct] contains power and frequency values for each channel
%    for each subject/session/phase
%
%    this script currently uses a global parameter set to determine where
%    to save the output figures

%% set up defaults

cfg_def = [];
cfg_def.type = 'both'; % whether to output the 'standard' or "white" filtered PSD
cfg_def.plot_type = 'raw';
cfg_def.pot_trk = '';
cfg_def.linewidth = 4;
cfg_def.color.blue = double([158,202,225])/255;
cfg_def.color.green = double([168,221,181])/255;
cfg_def.filter = [45 65; 70 90];
cfg = ProcessConfig2(cfg_def, cfg_in);
global PARAMS
c_ord = [linspecer(length(PARAMS.Phases)); [.6 .6 .6]];

%% Collect all the power ratios across the different sites per subject/session
if isempty(cfg.pot_trk)
    rec_type = {'pot', 'trk'};
else
    rec_type = cfg.pot_trk;
end

types = {'Pxx', 'White_Pxx', 'F', 'White_F'};
sites = {'PL'    'IL'    'OFC'    'Piri_O'    'NAc'    'Piri_N'    'CG'};

for iRec= 1:length(rec_type)
    for iType = 1:length(types)
        %         for iSite = 1:length(sites)
        for iPhase = 1:4
            all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})= [];
        end
        all_PSD.(rec_type{iRec}).control.(types{iType})= [];
        %         end
    end
    
    all_low.(rec_type{iRec}).Pxx = [];
    all_high.(rec_type{iRec}).Pxx = [];
    all_low.(rec_type{iRec}).White_Pxx = [];
    all_high.(rec_type{iRec}).White_Pxx = [];
end

% collect PSDs
subjects = fieldnames(Naris_in);
for iRec= 1:length(rec_type)
    for iSub = 1:length(subjects)
        sess_list = fieldnames(Naris_in.(subjects{iSub}));
        this_low.Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        this_high.Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        this_low.White_Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        this_high.White_Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        for iSess = 1:length(sess_list);
            for iSite = 1:length(sites)
                if sum(ismember(fieldnames(Naris_in.(subjects{iSub}).(sess_list{iSess}).pre), [sites{iSite} '_' rec_type{iRec}])) >0
                    for iPhase = 1:length(PARAMS.Phases)
                        for iType = 1:length(types)
                            current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).([sites{iSite} '_' rec_type{iRec}]).psd.(types{iType});
                            %                             all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}).(sites{iSite}) = cat(3,all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}).(sites{iSite}),current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess));
                        end
                        this_F = current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).F(iSite,:,iSess);
                        this_low.Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                        this_high.Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                        this_low.White_Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                        this_high.White_Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                    end
                    % set up control with average of pre and post
                    
                    for iType = 1:length(types)
                        current_PSD.(rec_type{iRec}).control.(types{iType})(iSite,:,iSess) = mean([Naris_in.(subjects{iSub}).(sess_list{iSess}).pre.([sites{iSite} '_' rec_type{iRec}]).psd.(types{iType}),Naris_in.(subjects{iSub}).(sess_list{iSess}).post.([sites{iSite} '_' rec_type{iRec}]).psd.(types{iType})],2);
                        %                         all_PSD.(rec_type{iRec}).control.(types{iType}).(sites{iSite}) = cat(3,all_PSD.(rec_type{iRec}).control.(types{iType}).(sites{iSite}),current_PSD.(rec_type{iRec}).control.(types{iType})(iSite,:,iSess));
                        
                    end
                    this_F = current_PSD.(rec_type{iRec}).control.F(iSite,:,iSess);
                    this_low.Pxx(iSite, 5, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).control.Pxx(iSite,nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                    this_high.Pxx(iSite, 5, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).control.Pxx(iSite,nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                    this_low.White_Pxx(iSite, 5, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).control.White_Pxx(iSite,nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                    this_high.White_Pxx(iSite, 5, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).control.White_Pxx(iSite,nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                else
                    for iType = 1:length(types)
                        for iPhase = 1:length(PARAMS.Phases)
                            current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess) = zeros(1,2049);
                        end
                            current_PSD.(rec_type{iRec}).control.(types{iType})(iSite,:,iSess) = zeros(1,2049);
                    end
                end
            end
        end
        for iType = 1:length(types)
            for iPhase = 1:length(PARAMS.Phases)
                all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}) = cat(3, all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}),current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}));
            end
            all_PSD.(rec_type{iRec}).control.(types{iType}) = cat(3, all_PSD.(rec_type{iRec}).control.(types{iType}),current_PSD.(rec_type{iRec}).control.(types{iType}));
        end
        
        
        
        all_low.(rec_type{iRec}).Pxx = cat(3,all_low.(rec_type{iRec}).Pxx, this_low.Pxx);
        all_high.(rec_type{iRec}).Pxx = cat(3,all_high.(rec_type{iRec}).Pxx, this_high.Pxx);
        all_low.(rec_type{iRec}).White_Pxx = cat(3,all_low.(rec_type{iRec}).White_Pxx, this_low.White_Pxx);
        all_high.(rec_type{iRec}).White_Pxx = cat(3,all_high.(rec_type{iRec}).White_Pxx, this_high.White_Pxx);
        
        clear current_PSD;
        
    end
end
all_low.labels = sites;
all_high.labels = sites;
%% deal with missing sites that will come out as zeros
for iRec= 1:length(rec_type)
    for iPhase = 1:length(PARAMS.Phases)
        for iType = 1:length(types)
            for iSite = 1:length(sites)
                for iSess = 1:size(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}),3)
                    if all(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess) ==0)
                        all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess) = NaN(1,length(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess)));
                    end
                end
            end
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot the average PSD across all sessions   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for this_Type = {'White', 'Pxx'};
    for iRec= 1:length(rec_type)
        for iSite = [1,2,3,4,5,6,7]
            h_site.(['n' num2str(iSite)]) = figure(10 +(iSite));
            y_val = [-150 200];
            rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
            rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
            for iPhase = 1:length(PARAMS.Phases)
                hold on
                if strcmp(this_Type, 'White')
                    this_mean_F = nanmean(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_F(iSite,:,:),3);
                    this_mean = nanmean(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,:,:)),3);
                    this_sem = nanstd(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,:,:)),0,3)./sqrt(size(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:),3));
                else
                    this_mean_F = nanmean(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).F(iSite,:,:),3);
                    this_mean = nanmean(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:)),3);
                    this_sem = nanstd(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:)),0,3)./sqrt(size(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:),3));
                end
                h1 = shadedErrorBar(this_mean_F, this_mean, this_sem);
                h1.mainLine.Color = c_ord(iPhase, :);
                h1.mainLine.LineWidth = cfg.linewidth;
                h1.patch.FaceColor = c_ord(iPhase, :);
                h1.patch.EdgeColor = c_ord(iPhase, :);
                h1.patch.FaceAlpha = .2;
                h1.patch.EdgeAlpha = .2;
            end
            xlim([0 100])
            if strcmp(this_Type, 'White')
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-140 -110])
                text(10, -138, sites{iSite}, 'fontsize', 48)
            else
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-120 -70])
                text(10, -118, sites{iSite}, 'fontsize', 48)
            end
            if iSite == 1
                set(gca, 'ytick', [-140 -110]); % put these back since it is the first figure. 
            end
            cfg_f.ft_size = 24;
            SetFigure(cfg_f, h_site.(['n' num2str(iSite)]))
            
            if ~isunix
                saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD\all_'  (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1}], 'png')
                saveas_eps([ 'all_' (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1}],[PARAMS.inter_dir 'PSD\'])
            else
                saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD/all_'  (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1}], 'png')
                saveas_eps([ 'all_'  (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1}],[PARAMS.inter_dir 'PSD/'])
            end
            close(h_site.(['n' num2str(iSite)]));
            %% traditional version
            h_site.(['n' num2str(iSite)]) = figure(100 +(iSite));
            rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.3]);
            rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3]);
            for iPhase = 1:length(PARAMS.Phases)
                hold on
                if strcmp(this_Type, 'White')
                    this_mean_F = nanmean(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_F(iSite,:,:),3);
                    this_mean = nanmean(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,:,:)),3);
                    this_sem = nanstd(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,:,:)),0,3);%./sqrt(size(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:),3));
                else
                    this_mean_F = nanmean(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).F(iSite,:,:),3);
                    this_mean = nanmean(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:)),3);
                    this_sem = nanstd(10*log10(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:)),0,3);%./sqrt(size(all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,:,:),3));
                end
                plot(this_mean_F,this_mean,'color', c_ord(iPhase,:), 'linewidth', cfg.linewidth);
            end
            
            
            H = get(gca, 'Children');
            set(gca, 'Children', [H(3), H(4), H(5), H(6), H(2), H(1)])
            clear H
            xlim([0 100])
            if strcmp(this_Type, 'White')
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-140 -110])
                text(10, -138, sites{iSite}, 'fontsize', 48)
            else
                set(gca, 'xtick', [0 100], 'ytick', [], 'ylim', [-120 -70])
                text(10, -118, sites{iSite}, 'fontsize', 48)
            end
            
            % move site text
            %                     legend(PARAMS.Phases)
            cfg_f.ft_size = 24;
            SetFigure(cfg_f, h_site.(['n' num2str(iSite)]))
            
            if ~isunix
                saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD\all_basic'  (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1} ], 'png')
                saveas_eps([ 'all_basic' (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1}],[PARAMS.inter_dir 'PSD\'])
            else
                saveas(h_site.(['n' num2str(iSite)]), [PARAMS.inter_dir 'PSD/all_basic'  (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1}], 'png')
                saveas_eps([ 'all_basic'  (rec_type{iRec}) '_' sites{iSite} '_' this_Type{1}],[PARAMS.inter_dir 'PSD/'])
            end
            close(h_site.(['n' num2str(iSite)]));
            
        end
    end
end
%% make a legend
figure(101)
hold on
for iPhase = 1:4
    plot(this_mean_F,this_mean,'color', c_ord(iPhase,:), 'linewidth', 4);
end
legend(PARAMS.Phases, 'location', 'south', 'orientation', 'horizontal');
legend boxoff
axis off
cfg_plt1.pos = [600 50 560*1.4 560*1.8];
cfg_plt1.ft_size = 18;
SetFigure(cfg_plt1, gcf)


if isunix
    %     saveas(gcf, [PARAMS.inter_dir '/AOC_fit/legend'], 'epsc')
    saveas_eps('legend_psd',[PARAMS.inter_dir '/AOC_fit/'])
else
    %     saveas(gcf, [PARAMS.inter_dir '\AOC_fit\legend'], 'epsc')
    saveas_eps('legend_psd',[PARAMS.inter_dir '\AOC_fit\'])
end
%% Stats for bar graph
for iRec= 1:length(rec_type)
    for iPhase = 1:length(PARAMS.Phases)+1 % cycle through pre, ipsi, contra, post and divide by the control to get a normalized value.
        all_norm_low.(rec_type{iRec}).Pxx(:,iPhase,:)  = all_low.(rec_type{iRec}).Pxx(:,iPhase,:)./all_low.(rec_type{iRec}).Pxx(:,5,:);
        all_norm_high.(rec_type{iRec}).Pxx(:,iPhase,:)  = all_high.(rec_type{iRec}).Pxx(:,iPhase,:)./all_high.(rec_type{iRec}).Pxx(:,5,:);
        all_norm_low.(rec_type{iRec}).White_Pxx(:,iPhase,:)  = all_low.(rec_type{iRec}).White_Pxx(:,iPhase,:)./all_low.(rec_type{iRec}).White_Pxx(:,5,:);
        all_norm_high.(rec_type{iRec}).White_Pxx(:,iPhase,:)  = all_high.(rec_type{iRec}).White_Pxx(:,iPhase,:)./all_high.(rec_type{iRec}).White_Pxx(:,5,:);
    end
    %     temp
end

%%
stats_file = fopen(['Raw_stats.txt'], 'w');
cfg_stats = [];
        cfg_stats.title = 'Raw_PSD_low_comp';
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri_O'  'NAc'  'Piri_N'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= [1 2 3 5 7];
        cfg_stats.ft_size= 20;
%         cfg_stats.save_dir= [PARAMS.inter_dir 'AOC_fit'];
        cfg_stats.stats_dir = stats_file;
        MS_stats(cfg_stats,permute(all_low.(rec_type{iRec}).White_Pxx,[2 1 3]));
        
%%
% make a bar plot
m_low.White_Pxx =  nanmean(all_norm_low.pot.White_Pxx(:,[2,3,5],:), 3);
m_high.White_Pxx =  nanmean(all_norm_high.pot.White_Pxx(:,[2,3,5],:), 3);
m_low.White_Pxx = circshift(m_low.White_Pxx,1,2);
m_high.White_Pxx = circshift(m_high.White_Pxx,1,2);


bar_names = sites;
c_ord = linspecer(7);
error_bar_low = (nanstd(all_norm_low.pot.White_Pxx(:,[2,3,5],:),1,3)./sqrt(length(all_norm_low.pot.White_Pxx(:,[2,3,5],:))))';
error_bar_high = (nanstd(all_norm_high.pot.White_Pxx(:,[2,3,5],:),1,3)./sqrt(length(all_norm_high.pot.White_Pxx(:,[2,3,5],:))))';
error_bar_low = circshift(error_bar_low,1,2);
error_bar_high = circshift(error_bar_high,1,2);

h_low = errorbar_groups(m_low.White_Pxx',error_bar_low, 'bar_colors', c_ord, 'bar_names', bar_names', 'FigID', 100);
h_high = errorbar_groups(m_high.White_Pxx',error_bar_high, 'bar_colors', c_ord,  'bar_names', bar_names', 'FigID', 200);

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Same thing but average in the gamma bands  %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % low gamma bar plot
% figure(200)
% figHandles = get(0,'Children');
% if sum(figHandles == 200) > 0
%     close(200)
% end
% cfg.conditions = {'control', 'ipsi', 'contra'};
% 
% F = figure(200);
% set(gcf, 'PaperPositionMode', 'auto', 'color', 'w')
% set(F, 'Position', [200, 200, 900 700])
% if isfield(cfg, 'whitefilter')
%     cfg.phases = cfg.conditions;
%     note = ('w/diff');
%     labels = {'Control', 'Ipsi', 'Contra'};
% else
%     cfg.phases = {'Pre', 'Ipsi', 'Contra', 'Post'};
%     labels = cfg.phases;
%     note = '';
% end
% hold on
% for iphase = 1:length(cfg.phases)
%     for id= 1:length(ids)
%         if id ==1 % "*iphase-(1/10) is to give a slight horizontal offset to the data for easier viewing.
%             h1 = plot(ones(4,3)*iphase-(1/10),comp_data.(cfg.phases{iphase})(:,id),  cfg.marker_ord{id}, 'MarkerEdgeColor', cfg.c_order_pow(id,:), 'MarkerSize', 14, 'LineWidth', 3) ;
%         elseif id==2
%             h2 = plot(ones(4,3)*iphase,comp_data.(cfg.phases{iphase})(:,id),  cfg.marker_ord{id}, 'MarkerEdgeColor', cfg.c_order_pow(id,:), 'MarkerSize', 14, 'LineWidth', 3);
%         elseif id ==3
%             h3 = plot(ones(4,3)*iphase+(1/10),comp_data.(cfg.phases{iphase})(:,id),  cfg.marker_ord{id}, 'MarkerEdgeColor', cfg.c_order_pow(id,:), 'MarkerSize', 14, 'LineWidth', 3);
%         end
%         disp(num2str(comp_data.(cfg.phases{iphase})(:,id)))
%     end
%     %         legend([h1(1), h2(1), h3(1)], 'R2', 'R4', 'R5', 'orientation', 'horizontal', 'fontsize', 16)
%     
%     legend([h1(1), h2(1), h3(1)],{'R5', 'R6', 'R7'}, 'location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 16)
%     line_of_best_fit(1,iphase) = mean(nanmean(comp_data.(cfg.phases{iphase}),1));
% end
% 
% plot(1:length(cfg.phases), line_of_best_fit, 'k', 'LineWidth', 1)
% % set the figure properties
% box on
% set(gcf, 'color', [1 1 1])
% ylabh = get(gca,'yLabel');
% % ylab_str('interpreter','latex','string','\fontsize{20}{0}\selectfont$Power$\fontsize{16}{0}\selectfont$(Normalized)$');
% % ylabel('\fontsize{20}{0}\selectfont$Power$\n\fontsize{16}{0}\selectfont$(Normalized)$','Interpreter','LaTex')
% ylabel('Power (normalized)')%, 'position', get(ylabh,'Position') - [.2 0 0] )
% text(2, -3, 'Contidtion', 'FontSize', 20)
% SetFigure([], gcf)





end
