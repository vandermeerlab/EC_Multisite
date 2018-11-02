function [out] = MS_stats(cfg_in, data_in)
%% MS_stats: takes in data in the form of a 3d array [phases, sites, sessions].
%
%
%
%
%   Input:
%       - cfg_in:  [struct] contains configuration paramters
%       - data_in: [3d array] in the form of [sites, phases, sessions]
%             - example:      [PL_pre, IL_pre, OFC_pre, ...]
%                             [PL_ipsi, IL_ipsi, OFC_ipsi, ...]
%                                   ...
%
%
%
%  Outputs:
%       - out: [struct]






%% defaults
% global PARAMS
cfg_def = [];
cfg_def.title = [];
cfg_def.method = 'median';
cfg_def.row_names = {'PL', 'IL', 'OFC', 'Piri_O', 'NAc', 'Piri_N', 'CG'};
cfg_def.col_names = {'pre', 'ipsi', 'contra','post'};
cfg_def.s_idx = 1:length(cfg_def.row_names); % corresponds to the sites to plot.
cfg_def.save_dir = cd; % just put it here unless otherwise specified with dir here.
cfg_def.stats_dir = []; % can be used to append to a text file
cfg = ProcessConfig2(cfg_def, cfg_in);

%% collect the mean and SD and transpose for grouping later
switch cfg.method
    case 'median'
        avg_vals = nanmedian(data_in, 3)';
    case 'mean'
        avg_vals = nanmean(data_in, 3)';
end

std_vals = nanstd(data_in, [], 3)';
SEM_vals = (nanstd(data_in,[],3)./sqrt(size(data_in,3)))';

fprintf(cfg.stats_dir,['**************** ' cfg.title ' ****************\n']);
fprintf(cfg.stats_dir,['**************** ' date ' ****************\n']);

fprintf(cfg.stats_dir, ['Using ' cfg.method '\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Get the stats %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% KS test for normality
ks = [];
for iR = 1:length(cfg.row_names)
    %     labels = {'ipsi', 'contra', 'control'};
    h = kstest(data_in(:, iR,:));
    if h
        disp('***************************************************************')
        disp(['KS test FAIL for low ' cfg.row_names{iR}])
        disp('***************************************************************')
        ks = [ks ; 1];
    end
    ks = [ks; 0];
end

%% test for differences
for iR = 1:length(cfg.row_names)
    %ipsi
    this_ipsi = squeeze(data_in(2, iR,:));
    this_ipsi(isnan(this_ipsi)) = [];
    % contra
    this_con= squeeze(data_in(3, iR,:));
    this_con(isnan(this_con)) = [];
    % control
    this_ctr = squeeze(data_in(5, iR,:));
    this_ctr(isnan(this_ctr)) = [];
    
    if sum(ks)>=1 || sum(ksh)>=1
        disp('Using Sign-Rank')
        [p_ip_con(iR), h_ip_con(iR)] = ranksum(this_ipsi, this_con);
        [p_ip_ctr(iR), h_ip_ctr(iR)] = ranksum(this_ipsi, this_ctr);
        [p_con_ctr(iR), h_con_ctr(iR)] = ranksum(this_con, this_ctr);
        
    else
        disp('Using T-Test')
        [h_ip_con(iR), p_ip_con(iR), ~, l_stats_ip_con(iR)] = ttest2(this_ipsi, this_con);
        [h_ip_ctr(iR), p_ip_ctr(iR), ~, l_stats_ip_ctr(iR)] = ttest2(this_ipsi, this_ctr);
        [h_con_ctr(iR), p_con_ctr(iR), ~,l_stats_con_ctr(iR)] = ttest2(this_con, this_ctr);
        
    end
end

%% display stats
if sum(ks) >=1
    fprintf(cfg.stats_dir,'\nWilcoxin Sign Rank test\n');
    fprintf(cfg.stats_dir,[cfg.title ': ' cfg.method '\n']);
    fprintf(cfg.stats_dir,'                                      ');
    for iR = 1:length(cfg.row_names)
        fprintf(cfg.stats_dir,[cfg.row_names{iR} '        '] );
        fprintf(cfg.stats_dir,repmat('\b', 1, length(cfg.row_names{iR})-2));
    end
    fprintf(cfg.stats_dir,['\nIpsilateral   vs. Contralateral:    P:' num2str(p_ip_con, '%10.4f') '\n' ]);
    fprintf(cfg.stats_dir,['Ipsilateral   vs. Control:          P:' num2str(p_ip_ctr, '%10.4f') '\n' ]);
    fprintf(cfg.stats_dir,['Contralateral vs. Control:          P:' num2str(p_con_ctr, '%10.4f') '\n' ]);
    
else
    fprintf('\nPaired T-Test\n')
    for iR = 1:length(cfg.row_names)
        fprintf(cfg.stats_dir,['\nLow Gamma   ' cfg.row_names{iR} '\n'])
        fprintf(cfg.stats_dir,[cfg.row_names{iR} ' Ipsilateral   vs. Contralateral:   df(' num2str(l_stats_ip_con(iR).df) ')   t:' num2str(l_stats_ip_con(iR).tstat, '%4.4f') '  P:' num2str(p_ip_con(iR), '%4.4f') '\n' ]);
        fprintf(cfg.stats_dir,[cfg.row_names{iR} ' Ipsilateral   vs. Control:         df(' num2str(l_stats_ip_ctr(iR).df) ')   t:' num2str(l_stats_ip_ctr(iR).tstat, '%4.4f') '  P:' num2str(p_ip_ctr(iR), '%4.4f') '\n' ]);
        fprintf(cfg.stats_dir,[cfg.row_names{iR} ' Contralateral vs. Control:         df(' num2str(l_stats_con_ctr(iR).df) ')   t:' num2str(l_stats_con_ctr(iR).tstat, '%4.4f') '  P:' num2str(p_con_ctr(iR), '%4.4f') '\n' ]);
    end
end


%% plot the output
% shift control to the first column
bar_c_ord = linspecer(3);
% select the row names to remove
bar_names= cfg.row_names;
bar_names(~ismember(1:length(cfg.row_names),cfg.s_idx)) = [];

bar_temp = circshift(avg_vals,1,2);
SEM_bar_temp= circshift(SEM_vals,1,2);


h_low = errorbar_groups(bar_temp(cfg.s_idx,[1,3,4])',SEM_bar_temp(cfg.s_idx,[1,3,4])', 'bar_colors', bar_c_ord, 'bar_names', bar_names, 'FigID', 100);
title([cfg.title])

%% add sig markers
%set up spacing
spacing_ip = 2:3:(3*length(bar_names));
spacing_con = spacing_ip + 0.9;
spacing_ctr = spacing_ip -0.9;
range =max(max(avg_vals))- min(min(avg_vals));
% use only the p and h values for the corresponding 
h_ip_con(~ismember(1:length(cfg.row_names),cfg.s_idx)) = [];
h_ip_ctr(~ismember(1:length(cfg.row_names),cfg.s_idx)) = [];
h_con_ctr(~ismember(1:length(cfg.row_names),cfg.s_idx)) = [];

p_ip_con(~ismember(1:length(cfg.row_names),cfg.s_idx)) = [];
p_ip_ctr(~ismember(1:length(cfg.row_names),cfg.s_idx)) = [];
p_con_ctr(~ismember(1:length(cfg.row_names),cfg.s_idx)) = [];

%
for iR  = 1:length(bar_names)
%     disp(bar_names{iR})
    hold on
    l_width = 1;
    ip_con = [];
    ip_con(1,:) = linspace(max(max(avg_vals)),max(max(avg_vals)), 50);
    ip_con(2,:) = ones(1,length(ip_con))*iR;
    
    % add the bars
    if h_ip_ctr(iR) == 1; % ipsi Vs control
        ip_con(1,:) = ip_con(1,:);
        plot(spacing_ip(iR), ip_con(1,:), '-k','linewidth', l_width);
        plot(spacing_ctr(iR), ip_con(1,:), '-k','linewidth', l_width)
        plot(linspace(spacing_ctr(iR),spacing_ip(iR), 50), ip_con(1,end)*ones(1,length(linspace(spacing_ctr(iR),spacing_ip(iR), 50))), '-k', 'linewidth', l_width)
        if p_ip_ctr(iR) <0.001;
            text(mean([spacing_ip(iR),spacing_ctr(iR)])-.5, ip_con(1,end)+.02, '***', 'FontSize', cfg.ft_size);
        elseif p_ip_ctr(iR) >0.00101 && p_ip_ctr(iR) <0.005;
            text(mean([spacing_ip(iR),spacing_ctr(iR)])-.3, ip_con(1,end)+.02, '**', 'FontSize', cfg.ft_size);
        elseif p_ip_ctr(iR) >0.0051 && p_ip_ctr(iR) <0.05;
            text(mean([spacing_ip(iR),spacing_ctr(iR)])-.1, ip_con(1,end)+.02, '*', 'FontSize', cfg.ft_size);
        end
    end
    
    if h_ip_con(iR) == 1; % ipsi Vs contra
        ip_con(1,:) = ip_con(1,:) +range*.05;
        plot(spacing_ip(iR), ip_con(1,:), '-k','linewidth', l_width);
        plot(spacing_con(iR), ip_con(1,:), '-k','linewidth', l_width)
        plot(linspace(spacing_ip(iR),spacing_con(iR), 50), ip_con(1,end)*ones(1,length(linspace(spacing_ip(iR),spacing_con(iR), 50))), '-k', 'linewidth', l_width)
        if p_ip_con(iR) <=0.001;
            text(mean([spacing_ip(iR),spacing_con(iR)])-.5, ip_con(1,end)+.02, '***', 'FontSize', cfg.ft_size);
        elseif p_ip_con(iR) >0.00101 && p_ip_con(iR) <0.005;
            text(mean([spacing_ip(iR),spacing_con(iR)])-.3, ip_con(1,end)+.02, '**', 'FontSize', cfg.ft_size);
        elseif p_ip_con(iR) >0.0051 && p_ip_con(iR) <0.05;
            text(mean([spacing_ip(iR),spacing_con(iR)])-.1, ip_con(1,end)+.02, '*', 'FontSize', cfg.ft_size);
        end
    end
    
    if h_con_ctr(iR) == 1; % contra Vs control
        ip_con(1,:) = ip_con(1,:) +range*.1;
        plot(spacing_ctr(iR), ip_con(1,:), '-k','linewidth', l_width);
        plot(spacing_con(iR), ip_con(1,:), '-k','linewidth', l_width)
        plot(linspace(spacing_ctr(iR),spacing_con(iR), 50), ip_con(1,end)*ones(1,length(linspace(spacing_ctr(iR),spacing_con(iR), 50))), '-k', 'linewidth', l_width)
        if p_con_ctr(iR) <=0.001;
            text(mean([spacing_ctr(iR),spacing_con(iR)])-.5, ip_con(1,end)+.02, '***', 'FontSize', cfg.ft_size);
        elseif p_con_ctr(iR) >0.00101 && p_con_ctr(iR) <0.005;
            text(mean([spacing_ctr(iR),spacing_con(iR)])-.3, ip_con(1,end)+.02, '**', 'FontSize', cfg.ft_size);
        elseif p_con_ctr(iR) >0.0051 && p_con_ctr(iR) <0.05;
            text(mean([spacing_ctr(iR),spacing_con(iR)])-.1, ip_con(1,end)+.02, '*', 'FontSize', cfg.ft_size);
        end
    end
    
end


%% set the default figure layout 
SetFigure([], gcf);


%% save the figure
      mkdir(cfg.save_dir)
      save_name = strrep(cfg.title, ' ', '_');
      if isunix
          fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '/' save_name '\n\n']);
          saveas(gcf, [ cfg.save_dir '/' save_name])
          saveas(gcf, [cfg.save_dir '/' save_name], 'png')
          saveas_eps(save_name,[cfg.save_dir '/'])
          %             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'epsc')
      else
          fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '\' save_name '\n\n']);
          saveas(gcf, [cfg.save_dir '\' save_name])
          saveas(gcf, [cfg.save_dir '\' save_name], 'png')
          saveas_eps(save_name,[cfg.save_dir '\'])
      end
      
      %close the text file
end



