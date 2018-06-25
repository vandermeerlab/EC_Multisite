function MS_Figure_S2(cfg_in, all_Naris)



%% set up parameters
global PARAMS

cfg_def = [];
cfg_def.measure = 'coh';
cfg_def.color.blue = double([158,202,225])/255;
cfg_def.color.green = double([168,221,181])/255;
cfg_def.linewidth = 4;
cfg_def.ylim = [0 1];
cfg_def.filter = [45 65; 70 90];
cfg_def.plot_type = 'no_piri';
cfg_def.legend = 'off';

cfg = ProcessConfig2(cfg_def, cfg_in);


%% get the average phase across each subject
for iPair = 1:length(PARAMS.all_pairs)
    for iPhase = 1:length(PARAMS.Phases)
        all_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
        all_fxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
        all_low_xx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
        all_high_xx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
    end
end

Subjects = fieldnames(all_Naris);
for iSub = 1:length(Subjects)
    sess_list = fieldnames(all_Naris.(Subjects{iSub}));
    for iSess = 1:length(sess_list)
        sub_pairs = fieldnames(all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx);
        for iPair = 1:length(sub_pairs)
            for iPhase = 1:length(PARAMS.Phases)
                if strcmp(cfg.measure, 'coh')
                    all_cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}));
                    
                    all_fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}));
                    
                    % get the mean within the gamma bands
                    F = all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase});
                    
                    m_low = all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(1,1), F):nearest_idx(cfg.filter(1,2), F));
                    m_high = all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(2,1), F):nearest_idx(cfg.filter(2,2), F));
                    
                    all_low_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_low_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_low));
                    all_high_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_high_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_high));
                    
                    %  get the mean per subject as well.
                    
                elseif strcmp(cfg.measure, 'amp')
                    % these need to be transposed
                    all_cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})');
                    
                    all_fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.f.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})');
                    
                    % get the mean within the gamma bands
                    F = all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.f.(sub_pairs{iPair}).(PARAMS.Phases{iPhase});
                    
                    m_low = all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(1,1), F):nearest_idx(cfg.filter(1,2), F));
                    m_high = all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(2,1), F):nearest_idx(cfg.filter(2,2), F));
                    
                    all_low_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_low_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_low));
                    all_high_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_high_xx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_high));
                    
                end
            end
        end
    end
end
%% get the mean values per site pair
for iPair = 1:length(PARAMS.all_pairs)
    for iPhase = 1:length(PARAMS.Phases)
        
        mean_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = mean(all_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}),2);
        std_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = std(all_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}),0,2);
        sem_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = std_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase})./sqrt(size(all_cxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}),2));
        mean_fxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}) = mean(all_fxx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}),2);
        
    end
end



%% plot ipsi vs contra for each pair with SEM
for iPair = 1:length(PARAMS.all_pairs)
    
    figure(iPair)
    hold on
    rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
    rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
    c_idx = find(strcmp(PARAMS.all_pairs,PARAMS.all_pairs{iPair}));
    
    h1 = shadedErrorBar(mean_fxx.(PARAMS.all_pairs{iPair}).contra, mean_cxx.(PARAMS.all_pairs{iPair}).contra, sem_cxx.(PARAMS.all_pairs{iPair}).contra);
    h1.mainLine.Color = PARAMS.pair_c_ord(c_idx, :);
    h1.mainLine.LineWidth = cfg.linewidth;
    h1.patch.FaceColor = PARAMS.pair_c_ord(c_idx, :);
    h1.patch.EdgeColor = PARAMS.pair_c_ord(c_idx, :);
    h1.patch.FaceAlpha = .5;
    
    h2 = shadedErrorBar(mean_fxx.(PARAMS.all_pairs{iPair}).ipsi, mean_cxx.(PARAMS.all_pairs{iPair}).ipsi, sem_cxx.(PARAMS.all_pairs{iPair}).ipsi);
    h2.mainLine.LineWidth = cfg.linewidth;
    h2.mainLine.Color = [.6 .6 .6];
    h2.patch.FaceColor =[.6 .6 .6];
    h2.patch.EdgeColor = [.6 .6 .6];
    h2.patch.FaceAlpha = .5;
    xlim([0 100])
    
    set(findall(gca, 'Type', 'Line'),'LineWidth',2)
    set(gca, 'ytick', [0 1], 'xtick',[0 100])
    ylabel([]); xlabel([]);
    cfg.set_fig.ft_size = 36;
    ylim(cfg.ylim)
    SetFigure(cfg.set_fig, gcf)
    
    mkdir([PARAMS.inter_dir 'sess'], 'all')
    if strcmp(cfg.measure, 'coh')
        saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_coh'], 'png')
        saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_coh'], 'fig')
        %         saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_coh'], 'epsc')
        saveas_eps(['all_' PARAMS.all_pairs{iPair} '_coh'], [PARAMS.inter_dir 'sess/all/'])
    elseif strcmp(cfg.measure, 'amp')
        saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_amp'], 'png')
        saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_amp'], 'fig')
        saveas_eps(['all_' PARAMS.all_pairs{iPair} '_amp'], [PARAMS.inter_dir 'sess/all/'])
        %         saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_amp'], 'epsc')
    end
    close all
    
    if strcmp(PARAMS.all_pairs{iPair}, 'OFC_NAc')
        figure(333)
        hold on
        rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
        rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
        c_idx = find(strcmp(PARAMS.all_pairs,PARAMS.all_pairs{iPair}));
        
        h1 = shadedErrorBar(mean_fxx.(PARAMS.all_pairs{iPair}).post, mean_cxx.(PARAMS.all_pairs{iPair}).post, sem_cxx.(PARAMS.all_pairs{iPair}).post);
        h1.mainLine.Color = PARAMS.pair_c_ord(c_idx, :);
        h1.mainLine.LineWidth = cfg.linewidth+4;
        h1.patch.FaceColor = PARAMS.pair_c_ord(c_idx, :);
        h1.patch.EdgeColor = PARAMS.pair_c_ord(c_idx, :);
        h1.patch.FaceAlpha = .5;
        h1.edge(1).Color = PARAMS.pair_c_ord(c_idx, :);
        h1.edge(2).Color = PARAMS.pair_c_ord(c_idx, :);
        
        set(findall(gca, 'Type', 'Line'),'LineWidth',2)
        h1.mainLine.LineWidth = cfg.linewidth+4;
        
        set(gca, 'ytick', [0 1], 'xtick',[0 100])
        xlim([0 100])
        ylabel([]); xlabel([]);
        cfg.set_fig.ft_size = 36;
        ylim(cfg.ylim)
        SetFigure(cfg.set_fig, gcf)
        %save the figure
        if strcmp(cfg.measure, 'coh')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_coh_OFC_NAc'], 'png')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_coh_OFC_NAc'], 'fig')
            %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
            saveas_eps(['all_sess_coh_OFC_NAc']',[PARAMS.inter_dir 'sess/'])
        elseif strcmp(cfg.measure, 'amp')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_amp_OFC_NAc'], 'png')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_amp_OFC_NAc'], 'fig')
            %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
            saveas_eps(['all_sess_amp_OFC_NAc']',[PARAMS.inter_dir 'sess/'])
        end
        
        set(findall(gca, 'Type', 'Line'),'LineWidth',6)
        ylabel([]); xlabel([]);
        set(gca, 'ytick', [0 1], 'xtick',[0 100])
        cfg.set_fig.ft_size = 36;
        SetFigure(cfg.set_fig, gcf)
        legend off
        %         saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc_small'], 'epsc')
        if strcmp(cfg.measure, 'coh')
            saveas_eps(['all_sess_coh_OFC_NAc_small']',[PARAMS.inter_dir 'sess/'])
        elseif strcmp(cfg.measure, 'amp')
            saveas_eps(['all_sess_amp_OFC_NAc_small']',[PARAMS.inter_dir 'sess/'])
        end
        close all
    elseif strcmp(PARAMS.all_pairs{iPair}, 'OFC_CG')
        figure(333)
        hold on
        rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
        rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
        c_idx = find(strcmp(PARAMS.all_pairs,PARAMS.all_pairs{iPair}));
        
        h1 = shadedErrorBar(mean_fxx.(PARAMS.all_pairs{iPair}).post, mean_cxx.(PARAMS.all_pairs{iPair}).post, sem_cxx.(PARAMS.all_pairs{iPair}).post);
        h1.mainLine.Color = PARAMS.pair_c_ord(c_idx, :);
        h1.mainLine.LineWidth = cfg.linewidth+4;
        h1.patch.FaceColor = PARAMS.pair_c_ord(c_idx, :);
        h1.patch.EdgeColor = PARAMS.pair_c_ord(c_idx, :);
        h1.patch.FaceAlpha = .5;
        h1.edge(1).Color = PARAMS.pair_c_ord(c_idx, :);
        h1.edge(2).Color = PARAMS.pair_c_ord(c_idx, :);
        
        set(findall(gca, 'Type', 'Line'),'LineWidth',2)
        h1.mainLine.LineWidth = cfg.linewidth+4;
        
        set(gca, 'ytick', [0 1], 'xtick',[0 100])
        xlim([0 100])
        ylabel([]); xlabel([]);
        cfg.set_fig.ft_size = 36;
        ylim(cfg.ylim)
        SetFigure(cfg.set_fig, gcf)
        %save the figure
        if strcmp(cfg.measure, 'coh')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_coh_OFC_CG'], 'png')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_coh_OFC_CG'], 'fig')
            %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
            saveas_eps(['all_sess_coh_OFC_CG']',[PARAMS.inter_dir 'sess/'])
        elseif strcmp(cfg.measure, 'amp')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_amp_OFC_CG'], 'png')
            saveas(gcf, [PARAMS.inter_dir 'sess/' 'all_sess_mean_amp_OFC_CG'], 'fig')
            %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
            saveas_eps(['all_sess_amp_OFC_CG']',[PARAMS.inter_dir 'sess/'])
        end
        set(findall(gca, 'Type', 'Line'),'LineWidth',6)
        ylabel([]); xlabel([]);
        set(gca, 'ytick', [0 1], 'xtick',[0 100])
        cfg.set_fig.ft_size = 36;
        SetFigure(cfg.set_fig, gcf)
        legend off
        %         saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc_small'], 'epsc')
        if strcmp(cfg.measure, 'coh')
            saveas_eps(['all_sess_coh_OFC_CG_small']',[PARAMS.inter_dir 'sess/'])
        elseif strcmp(cfg.measure, 'amp')
            saveas_eps(['all_sess_amp_OFC_CG_small']',[PARAMS.inter_dir 'sess/'])
        end
        close all
    end
    
end
%% create the matrix figures for Fig 6

labels = PARAMS.all_sites;
if strcmp(cfg.plot_type, 'no_piri')
    isempty(strfind(labels, 'Piri'))
    labels(strcmp('PiriO', labels)) = [];
    labels(strcmp('PiriN', labels)) = [];
end



for ii = 1:size(labels,2)
    for jj = 1:size(labels,2)
        Idx_mat(ii, jj) = str2double([num2str(ii) num2str(jj)]);
    end
end

mat_out =[];
mat_labels = cell(length(labels));

for iPair = 1:length(PARAMS.all_pairs)
    
    if isempty(strfind(PARAMS.all_pairs{iPair}, 'Piri'))
        
        for iPhase = 1:length(PARAMS.Phases)
            
            S = strsplit(PARAMS.all_pairs{iPair}, '_');
            
            ii = strfind(labels, S{1});
            ii = find(not(cellfun('isempty', ii)));
            
            jj = strfind(labels, S{2});
            jj = find(not(cellfun('isempty', jj)));
            
            
            if strcmp(PARAMS.all_pairs{iPair}, 'IL_PL')
                old = [ii, jj];
                ii = old(2);
                jj = old(1);
            end
            mat_out.(PARAMS.Phases{iPhase})(ii, jj) = mean(all_high_xx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}));
            mat_out.(PARAMS.Phases{iPhase})(jj, ii) = mean(all_low_xx.(PARAMS.all_pairs{iPair}).(PARAMS.Phases{iPhase}));
            
        end
        mat_labels{ii, jj} = [S{1} '_' S{2} '_h'];
        mat_labels{jj, ii} = [S{1} '_' S{2} '_l'];
    end
end
%% plot the matrix of mean values
Phase_c = linspecer(4);
figure(1111)
for iPhase = 1:length(PARAMS.Phases)
    mat_out.(PARAMS.Phases{iPhase})(logical(eye(size(mat_out.(PARAMS.Phases{iPhase}))))) = NaN;
    
    
    subplot(1,4,iPhase)
    nan_imagesc_ec(mat_out.(PARAMS.Phases{iPhase}), 'nan_colour', 'w')
    
    
    [~, hStrings] =add_num_imagesc(gca,mat_out.(PARAMS.Phases{iPhase}), 2,  12); % adds numerics to imagesc
    
    colormap('parula')
    
    %     set(hcb, 'Ytick', 0:0.1:.5)
    set(gca,'xtick', 1:length(labels), 'ytick', 1:length(labels), 'xaxisLocation','top', 'xticklabel', labels, 'yticklabel', labels, 'XTickLabelRotation', 65)
    
    rectangle('position', [0.5, 0.5, length(labels), length(labels)], 'EdgeColor', Phase_c(iPhase,:), 'linewidth', 4)
    caxis([0 .6]);
end
%
Square_subplots
tightfig
cfg_count_fig.resize = 1;
cfg_count_fig.pos = [0 401 1440 380];
cfg_count_fig.ft_size = 14;
SetFigure(cfg_count_fig, gcf)



if strcmp(cfg.measure, 'coh')
    saveas(gcf, [PARAMS.inter_dir 'sess/all/all_coh_mat'], 'png')
    saveas(gcf, [PARAMS.inter_dir 'sess/all/all_coh_mat'], 'fig')
    saveas_eps('all_coh_mat',[PARAMS.inter_dir 'sess/all/'])
    
    %     saveas(gcf, [PARAMS.inter_dir 'sess/all/all_coh_mat'], 'epsc')
    %     pushdir([PARAMS.inter_dir 'sess/all/'])
    %     eval(sprintf('print -depsc2 -tiff  -r300 -painters %s','all_coh_mat'));
    %     popdir
elseif strcmp(cfg.measure, 'amp')
    saveas(gcf, [PARAMS.inter_dir 'sess/all/all_amp_mat'], 'png')
    saveas(gcf, [PARAMS.inter_dir 'sess/all/all_amp_mat'], 'fig')
    %     saveas(gcf, [PARAMS.inter_dir 'sess/all/all_amp_mat'], 'epsc')
    
    saveas_eps('all_amp_mat',[PARAMS.inter_dir 'sess/all/'])
    %     pushdir([PARAMS.inter_dir 'sess/all/'])
    %     eval(sprintf('print -depsc2 -tiff  -r300 -painters %s','all_amp_mat'));
    %     popdir
    
end
close all


%% same thing for within subject averages

for iSub = 1:length(Subjects)
    sess_list = fieldnames(all_Naris.(Subjects{iSub}));
    %     for iSess = 1:length(sess_list)
    sub_pairs = fieldnames(all_Naris.(Subjects{iSub}).(sess_list{1}).coh.cxx);
    % pre allocate
    for iPair = 1:length(sub_pairs)
        for iPhase = 1:length(PARAMS.Phases)
            all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
            all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
            all_sub_low_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
            all_sub_high_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = [];
        end
    end
    
    for iPair = 1:length(sub_pairs)
        for iPhase = 1:length(PARAMS.Phases)
            for iSess = 1:length(sess_list)
                if strcmp(cfg.measure, 'coh')
                    all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}));
                    
                    all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase}));
                    
                    % get the mean within the gamma bands
                    F = all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.fxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase});
                    
                    m_low = all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(1,1), F):nearest_idx(cfg.filter(1,2), F));
                    m_high = all_Naris.(Subjects{iSub}).(sess_list{iSess}).coh.cxx.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(2,1), F):nearest_idx(cfg.filter(2,2), F));
                    
                    all_sub_low_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_low_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_low));
                    all_sub_high_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_high_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_high));
                    
                    %  get the mean per subject as well.
                    
                elseif strcmp(cfg.measure, 'amp')
                    % these need to be transposed
                    all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})');
                    
                    all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),   all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.f.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})');
                    
                    % get the mean within the gamma bands
                    F = all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.f.(sub_pairs{iPair}).(PARAMS.Phases{iPhase});
                    
                    m_low = all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(1,1), F):nearest_idx(cfg.filter(1,2), F));
                    m_high = all_Naris.(Subjects{iSub}).(sess_list{iSess}).amp.ac.(sub_pairs{iPair}).(PARAMS.Phases{iPhase})(nearest_idx(cfg.filter(2,1), F):nearest_idx(cfg.filter(2,2), F));
                    
                    all_sub_low_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_low_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_low));
                    all_sub_high_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = cat(2,all_sub_high_xx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}), mean(m_high));
                end
            end
        end
    end
end
%% get the mean values per site pair for each subject
for iSub = 1:length(Subjects)
    sub_pairs = fieldnames(all_sub_cxx.(Subjects{iSub}));
    for iPair = 1:length(sub_pairs)
        for iPhase = 1:length(PARAMS.Phases)
            
            all_sub_mean_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = mean(all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),2);
            all_sub_std_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = std(all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),0,2);
            all_sub_sem_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = all_sub_std_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase})./sqrt(size(all_sub_cxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),2));
            all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}) = mean(all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).(PARAMS.Phases{iPhase}),2);
            
        end
    end
end
%% plot the mean for each subject.
for iSub =1:length(Subjects)
    sub_pairs = fieldnames(all_sub_mean_cxx.(Subjects{iSub}));
    
    for iPair = 1:length(sub_pairs)
        
        figure(iSub)
        hold on
        rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
        rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
        c_idx = find(strcmp(PARAMS.all_pairs,sub_pairs{iPair}));
        
        h1 = shadedErrorBar(all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).post, all_sub_mean_cxx.(Subjects{iSub}).(sub_pairs{iPair}).post, all_sub_sem_cxx.(Subjects{iSub}).(sub_pairs{iPair}).post);
        h1.mainLine.Color = PARAMS.pair_c_ord(c_idx, :);
        h1.mainLine.LineWidth = cfg.linewidth;
        h1.patch.FaceColor = PARAMS.pair_c_ord(c_idx, :);
        h1.patch.EdgeColor = PARAMS.pair_c_ord(c_idx, :);
        h1.patch.FaceAlpha = .5;
        h1.edge(1).Color = PARAMS.pair_c_ord(c_idx, :);
        h1.edge(2).Color = PARAMS.pair_c_ord(c_idx, :);
        
        %     h2 = shadedErrorBar(mean_fxx.(PARAMS.all_pairs{iPair}).ipsi, mean_cxx.(PARAMS.all_pairs{iPair}).ipsi, sem_cxx.(PARAMS.all_pairs{iPair}).ipsi);
        %     h2.mainLine.LineWidth = cfg.linewidth;
        %     h2.mainLine.Color = [.6 .6 .6];
        %     h2.patch.FaceColor =[.6 .6 .6];
        %     h2.patch.EdgeColor = [.6 .6 .6];
        %     h2.patch.FaceAlpha = .5;
        xlim([0 100])
        
        set(findall(gca, 'Type', 'Line'),'LineWidth',2)
        set(gca, 'ytick', [0 1], 'xtick',[0 100])
        ylabel([]); xlabel([]);
        cfg.set_fig.ft_size = 36;
        ylim(cfg.ylim)
        SetFigure(cfg.set_fig, gcf)
        
        mkdir([PARAMS.inter_dir 'sess'], 'coh')
        mkdir([PARAMS.inter_dir 'sess'], 'amp')
        
        if strcmp(cfg.measure, 'coh')
            saveas(gcf, [PARAMS.inter_dir 'sess/coh/' (Subjects{iSub}) '_' sub_pairs{iPair} '_coh'], 'png')
            saveas(gcf, [PARAMS.inter_dir 'sess/coh/' (Subjects{iSub}) '_' sub_pairs{iPair} '_coh'], 'fig')
            %         saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_coh'], 'epsc')
            saveas_eps([(Subjects{iSub}) '_' PARAMS.all_pairs{iPair} '_coh'], [PARAMS.inter_dir 'sess/coh/'])
        elseif strcmp(cfg.measure, 'amp')
            saveas(gcf, [PARAMS.inter_dir 'sess/amp/' (Subjects{iSub}) '_' sub_pairs{iPair} '_amp'], 'png')
            saveas(gcf, [PARAMS.inter_dir 'sess/amp/' (Subjects{iSub}) '_' sub_pairs{iPair} '_amp'], 'fig')
            saveas_eps([(Subjects{iSub}) '_' PARAMS.all_pairs{iPair} '_amp'], [PARAMS.inter_dir 'sess/amp/'])
            %         saveas(gcf, [PARAMS.inter_dir 'sess/all/all_' PARAMS.all_pairs{iPair} '_amp'], 'epsc')
        end
        close all
        
        if strcmp(sub_pairs{iPair}, 'OFC_NAc')
            figure(333)
            hold on
            rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
            rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
            c_idx = find(strcmp(PARAMS.all_pairs,sub_pairs{iPair}));
            
            h1 = shadedErrorBar(all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).post, all_sub_mean_cxx.(Subjects{iSub}).(sub_pairs{iPair}).post, all_sub_sem_cxx.(Subjects{iSub}).(sub_pairs{iPair}).post);
            h1.mainLine.Color = PARAMS.pair_c_ord(c_idx, :);
            h1.patch.FaceColor = PARAMS.pair_c_ord(c_idx, :);
            h1.patch.EdgeColor = PARAMS.pair_c_ord(c_idx, :);
            h1.patch.FaceAlpha = .5;
            h1.edge(1).Color = PARAMS.pair_c_ord(c_idx, :);
            h1.edge(2).Color = PARAMS.pair_c_ord(c_idx, :);
            
            set(findall(gca, 'Type', 'Line'),'LineWidth',2)
            h1.mainLine.LineWidth = cfg.linewidth+4;
            
            set(gca, 'ytick', [0 1], 'xtick',[0 100])
            xlim([0 100])
            ylabel([]); xlabel([]);
            cfg.set_fig.ft_size = 36;
            ylim(cfg.ylim)
            SetFigure(cfg.set_fig, gcf)
            %save the figure
            if strcmp(cfg.measure, 'coh')
                saveas(gcf, [PARAMS.inter_dir 'sess/coh/' (Subjects{iSub}) 'mean_coh_OFC_NAc'], 'png')
                saveas(gcf, [PARAMS.inter_dir 'sess/coh/' (Subjects{iSub}) 'mean_coh_OFC_NAc'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
                saveas_eps([(Subjects{iSub}) 'coh_OFC_NAc'],[PARAMS.inter_dir 'sess/coh/'])
            elseif strcmp(cfg.measure, 'amp')
                saveas(gcf, [PARAMS.inter_dir 'sess/amp/' (Subjects{iSub}) 'mean_amp_OFC_NAc'], 'png')
                saveas(gcf, [PARAMS.inter_dir 'sess/amp/' (Subjects{iSub}) 'mean_amp_OFC_NAc'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
                saveas_eps([(Subjects{iSub}) 'amp_OFC_NAc'],[PARAMS.inter_dir 'sess/amp/'])
            end
            
            set(findall(gca, 'Type', 'Line'),'LineWidth',6)
            ylabel([]); xlabel([]);
            set(gca, 'ytick', [0 1], 'xtick',[0 100])
            cfg.set_fig.ft_size = 36;
            SetFigure(cfg.set_fig, gcf)
            legend off
            %         saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc_small'], 'epsc')
            if strcmp(cfg.measure, 'coh')
                saveas_eps([(Subjects{iSub}) 'coh_OFC_NAc_small'],[PARAMS.inter_dir 'sess/coh/'])
            elseif strcmp(cfg.measure, 'amp')
                saveas_eps([(Subjects{iSub}) 'amp_OFC_NAc_small'],[PARAMS.inter_dir 'sess/amp/'])
            end
            close all
        end
        if strcmp(sub_pairs{iPair}, 'OFC_CG')
            figure(333)
            hold on
            rectangle('position', [45, 0.001, 20, 1], 'facecolor', [cfg.color.blue 0.3], 'edgecolor', [cfg.color.blue 0.5])
            rectangle('position', [70, 0.001, 20, 1], 'facecolor', [cfg.color.green 0.3], 'edgecolor', [cfg.color.green 0.3])
            c_idx = find(strcmp(PARAMS.all_pairs,sub_pairs{iPair}));
            
            h1 = shadedErrorBar(all_sub_fxx.(Subjects{iSub}).(sub_pairs{iPair}).post, all_sub_mean_cxx.(Subjects{iSub}).(sub_pairs{iPair}).post, all_sub_sem_cxx.(Subjects{iSub}).(sub_pairs{iPair}).post);
            h1.mainLine.Color = PARAMS.pair_c_ord(c_idx, :);
            h1.patch.FaceColor = PARAMS.pair_c_ord(c_idx, :);
            h1.patch.EdgeColor = PARAMS.pair_c_ord(c_idx, :);
            h1.patch.FaceAlpha = .5;
            h1.edge(1).Color = PARAMS.pair_c_ord(c_idx, :);
            h1.edge(2).Color = PARAMS.pair_c_ord(c_idx, :);
            
            set(findall(gca, 'Type', 'Line'),'LineWidth',2)
            h1.mainLine.LineWidth = cfg.linewidth+4;
            
            set(gca, 'ytick', [0 1], 'xtick',[0 100])
            xlim([0 100])
            ylabel([]); xlabel([]);
            cfg.set_fig.ft_size = 36;
            ylim(cfg.ylim)
            SetFigure(cfg.set_fig, gcf)
            %save the figure
            if strcmp(cfg.measure, 'coh')
                saveas(gcf, [PARAMS.inter_dir 'sess/coh/' (Subjects{iSub}) 'mean_coh_OFC_CG'], 'png')
                saveas(gcf, [PARAMS.inter_dir 'sess/coh/' (Subjects{iSub}) 'mean_coh_OFC_CG'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
                saveas_eps([(Subjects{iSub}) 'coh_OFC_CG'],[PARAMS.inter_dir 'sess/coh/'])
            elseif strcmp(cfg.measure, 'amp')
                saveas(gcf, [PARAMS.inter_dir 'sess/amp/' (Subjects{iSub}) 'mean_amp_OFC_CG'], 'png')
                saveas(gcf, [PARAMS.inter_dir 'sess/amp/' (Subjects{iSub}) 'mean_amp_OFC_CG'], 'fig')
                %                 saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc'], 'epsc')
                saveas_eps([(Subjects{iSub}) 'amp_OFC_CG'],[PARAMS.inter_dir 'sess/amp/'])
            end
            
            set(findall(gca, 'Type', 'Line'),'LineWidth',6)
            ylabel([]); xlabel([]);
            set(gca, 'ytick', [0 1], 'xtick',[0 100])
            cfg.set_fig.ft_size = 36;
            SetFigure(cfg.set_fig, gcf)
            legend off
            %         saveas(gcf, [PARAMS.inter_dir '/sess/' Subjects{iSub} '_sess_amp_OFC_NAc_small'], 'epsc')
            if strcmp(cfg.measure, 'coh')
                saveas_eps([(Subjects{iSub}) 'coh_OFC_CG_small'],[PARAMS.inter_dir 'sess/coh/'])
            elseif strcmp(cfg.measure, 'amp')
                saveas_eps([(Subjects{iSub}) 'amp_OFC_CG_small'],[PARAMS.inter_dir 'sess/amp/'])
            end
            close all
        end
        
    end
end


end



