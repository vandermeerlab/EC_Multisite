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
global PARAMS
cfg_def = [];
cfg_def.title = [];
cfg_def.method = 'median';
cfg_def.NaN_correct = 0; % can be used to correct for NaNs (work around for MS_get_naris_phase_distance
cfg_def.row_names = {'PL', 'IL', 'OFC', 'PiriO', 'NAc', 'PiriN', 'CG'};
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

% get the number of valid session for computing the 
for ii = size(data_in,1):-1:1;
    for jj = size(data_in,2):-1:1;
        this_site = data_in(ii,jj,:); 
        size_mat(ii, jj) = length(this_site(~isnan(this_site)));
        sqrt_mat(ii, jj) = sqrt(size_mat(ii, jj)); 
%         std_mat(ii,jj) = nanstd(this_site); 
%         SEM_vals(ii,jj) = std_mat(ii,jj)./sqrt_mat(ii,jj); % double check. 
    end
end
% SEM_vals = SEM_vals'; 
SEM_vals = (nanstd(data_in,[],3)./sqrt_mat)';

% SEM_vals = (nanstd(data_in,[],3)./sqrt(size(data_in,3)))';

fprintf(cfg.stats_dir,['**************** ' cfg.title ' ****************\n']);
fprintf(cfg.stats_dir,['**************** ' date ' ****************\n']);

fprintf(cfg.stats_dir, ['Using ' cfg.method '\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Get the stats %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(cfg.stats_dir, ['Analyzing with ' cfg.stats_method '\n']);
if strcmp(cfg.stats_method, 'sign_rank')
    
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
            [p_ip_con(iR), h_ip_con(iR)] = signrank(this_ipsi, this_con);
            [p_ip_ctr(iR), h_ip_ctr(iR)] = signrank(this_ipsi, this_ctr);
            [p_con_ctr(iR), h_con_ctr(iR)] = signrank(this_con, this_ctr);
            
        else
            disp('Using T-Test')
            [h_ip_con(iR), p_ip_con(iR), ci_ip_con(iR), l_stats_ip_con(iR)] = ttest2(this_ipsi, this_con);
            [h_ip_ctr(iR), p_ip_ctr(iR), ci_ip_con(iR), l_stats_ip_ctr(iR)] = ttest2(this_ipsi, this_ctr);
            [h_con_ctr(iR), p_con_ctr(iR), ci_ip_con(iR),l_stats_con_ctr(iR)] = ttest2(this_con, this_ctr);
            
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
    
elseif strcmp(cfg.stats_method, 'lme')
    
    %% try with an LME
    % prepare data for LME format
    % build sess, subject, condition 3D arrays
    sess_id = ones(size(data_in));
    for iSub = 1:4:size(sess_id,3)
        for ii = 0:3
            sess_id(:,:,iSub+ii) = ones(size(sess_id,1), size(sess_id,2)).*ii+1;
        end
    end
    
    subject_id = [];
    for iSub = 1:4:size(sess_id,3)
        subject_id = cat(3,subject_id, repmat(iSub, size(data_in,1), size(data_in,2), 4));
    end
    
    % condition
    condition_id = cell(size(data_in));
    for ii =1 :size(data_in,2)
        for jj = 1:size(data_in,3)
            condition_id(:,ii,jj) = {'pre', 'ispi', 'contra', 'post','control'};
        end
    end
    %% loop over targets.  One target per LME for ipsi_v_contra
            fprintf(cfg.stats_dir, '\n<strong>%s</strong>', 'Ipsi-Contra'); 

    for iSite = 1:length(cfg.row_names)
        sess_1d = reshape(squeeze(sess_id(1,iSite,:)), 1, size(sess_id,3));
        subject_1d = reshape(squeeze(subject_id(1,iSite,:)), 1, size(subject_id,3));
        ipsi_1d = reshape(squeeze(data_in(2,iSite,:)), 1, size(data_in,3));
        contra_1d = reshape(squeeze(data_in(3,iSite,:)), 1,  size(data_in,3));
        ipsi_id_1d = reshape(squeeze(condition_id(2,iSite,:)), 1,  size(data_in,3));
        contra_id_1d = reshape(squeeze(condition_id(3,iSite,:)), 1,  size(data_in,3));
        % cat all the data together into a 1d array for each variable
        sess_lme = [sess_1d, sess_1d];
        subject_lme = [subject_1d, subject_1d];
        
        power_lme = [ipsi_1d, contra_1d];
        condition_lme = [ipsi_id_1d, contra_id_1d];
%         power_lme = [contra_1d,ipsi_1d];
%         condition_lme = [contra_id_1d,ipsi_id_1d];
        
        %% build the LME for this site
        D_power.tbl = table(subject_lme', sess_lme', power_lme', condition_lme','VariableNames',{'SubjectID','SessID', 'Power', 'Condition'});
        D_power.tbl.SubjectID = nominal(D_power.tbl.SubjectID);
        D_power.tbl.SessID = nominal(D_power.tbl.SessID);
        D_power.tbl.Condition = nominal(D_power.tbl.Condition);
        
        % remove NaNs
        if cfg.NaN_correct == 1
       D_power.tbl(~any(ismissing(D_power.tbl),2),:)
        end
        
        % disp([cfg.row_names{iSite} '___________________'])
        D_power.lme = fitlme(D_power.tbl,'Power~1+Condition+(1|SubjectID)+(1|SessID)');
        % D_power.lme
        % anova(D_power.lme,'DFMethod','satterthwaite')
        
        % collect values
        Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Est =  D_power.lme.Coefficients.Estimate(2);
        Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).SE =  D_power.lme.Coefficients.SE(2);
        Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val =  D_power.lme.Coefficients.pValue(2);
        Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Lower =  D_power.lme.Coefficients.Lower(2);
        Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Upper =  D_power.lme.Coefficients.Upper(2);
        Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).df =  D_power.lme.Coefficients.DF(2);

        % T-stats
        Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).tstat =  D_power.lme.Coefficients.tStat(2);
        
        % hold the P value for plotting later
        p_ip_con(iSite) = Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val;
        
        
        fprintf(cfg.stats_dir,['\n' PARAMS.all_sites{iSite} '& Ipsi-Contra  ']);
        if Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val >= 0.05
            p_val_string = ['p$=$ ' num2str(Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val,'%10.2e')];
            h_ip_con(iSite) = 0; % use for assigning markers later.
            
        elseif (0.049999 > Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val >= 0.01);
            p_val_string = '*';
            h_ip_con(iSite) = 1;
            
        elseif (0.009999 > Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val >= 0.001);
            p_val_string = '**';
            h_ip_con(iSite) = 1;
            
        elseif 0.0009999 > Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val
            p_val_string = '***';
            h_ip_con(iSite) = 1;
        end
        
        
            fprintf(cfg.stats_dir,'& %4.2f & (%2.0f) %4.2f %s & %4.2f & %4.2f',...
            Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Est,...
            Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).df,...
            Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).tstat,...
            p_val_string,...
            Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Lower,...
            Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Upper); 
        
        
    end
    
    %% loop over targets.  One target per LME for ipsi_v_control
    
    fprintf(cfg.stats_dir, '\n<strong>%s</strong>', 'Ipsi-Control'); 
    for iSite = 1:length(cfg.row_names)
        sess_1d = reshape(squeeze(sess_id(1,iSite,:)), 1, size(sess_id,3));
        subject_1d = reshape(squeeze(subject_id(1,iSite,:)), 1, size(subject_id,3));
        ipsi_1d = reshape(squeeze(data_in(2,iSite,:)), 1, size(data_in,3));
        ctrl_1d = reshape(squeeze(data_in(5,iSite,:)), 1,  size(data_in,3));
        ipsi_id_1d = reshape(squeeze(condition_id(2,iSite,:)), 1,  size(data_in,3));
        ctrl_id_1d = reshape(squeeze(condition_id(5,iSite,:)), 1,  size(data_in,3));
        % cat all the data together into a 1d array for each variable
        sess_lme = [sess_1d, sess_1d];
        subject_lme = [subject_1d, subject_1d];
        power_lme = [ipsi_1d, ctrl_1d];
        condition_lme = [ipsi_id_1d, ctrl_id_1d];
        

        %% build the LME for this site
        D_power.tbl = table(subject_lme', sess_lme', power_lme', condition_lme','VariableNames',{'SubjectID','SessID', 'Power', 'Condition'});
        D_power.tbl.SubjectID = nominal(D_power.tbl.SubjectID);
        D_power.tbl.SessID = nominal(D_power.tbl.SessID);
        D_power.tbl.Condition = nominal(D_power.tbl.Condition);
        
        % disp([cfg.row_names{iSite} '___________________'])
        D_power.lme = fitlme(D_power.tbl,'Power~1+Condition+(1|SubjectID)+(1|SessID)');
        % D_power.lme
        % anova(D_power.lme,'DFMethod','satterthwaite')
        
        % collect values
        Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).Est =  D_power.lme.Coefficients.Estimate(2);
        Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).SE =  D_power.lme.Coefficients.SE(2);
        Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val =  D_power.lme.Coefficients.pValue(2);
        Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).Lower =  D_power.lme.Coefficients.Lower(2);
        Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).Upper =  D_power.lme.Coefficients.Upper(2);
        Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).df =  D_power.lme.Coefficients.DF(2);

        %tstats
        Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).tstat =  D_power.lme.Coefficients.tStat(2);

        
        % hold the P value for plotting later
        p_ip_ctr(iSite) = Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val;
        
        
        fprintf(cfg.stats_dir,['\n' PARAMS.all_sites{iSite} '& Ipsi-Control  ']);
        if Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val >= 0.05
            p_val_string = ['p$=$ ' num2str(Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val,'%10.2e')];
            
            h_ip_ctr(iSite) = 0; % use for assigning markers later.
            
        elseif (0.0499999 > Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val >= 0.01);
            p_val_string = '*';
            h_ip_ctr(iSite) = 1;
            
        elseif (0.0099999> Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val >= 0.001);
            p_val_string = '**';
            h_ip_ctr(iSite) = 1;
            
        elseif 0.0009999 > Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).P_val
            p_val_string = '***';
            h_ip_ctr(iSite) = 1;
            
        end
        
            fprintf(cfg.stats_dir,'& %4.2f & (%2.0f) %4.2f %s & %4.2f & %4.2f',...
            Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).Est,...
            Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).df,...
            Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).tstat,...
            p_val_string,...
            Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).Lower,...
            Stats_out.ipsi_ctrl.(PARAMS.all_sites{iSite}).Upper); % p $<$ 0.001
        
    end
    %% loop over targets.  One target per LME for Control_v_contra
        fprintf(cfg.stats_dir, '\n<strong>%s</strong>', 'Contra-Control'); 

    for iSite = 1:length(cfg.row_names)
        sess_1d = reshape(squeeze(sess_id(1,iSite,:)), 1, size(sess_id,3));
        subject_1d = reshape(squeeze(subject_id(1,iSite,:)), 1, size(subject_id,3));
        contra_1d = reshape(squeeze(data_in(3,iSite,:)), 1, size(data_in,3));
        ctrl_1d = reshape(squeeze(data_in(5,iSite,:)), 1,  size(data_in,3));
        contra_id_1d = reshape(squeeze(condition_id(3,iSite,:)), 1,  size(data_in,3));
        ctrl_id_1d = reshape(squeeze(condition_id(5,iSite,:)), 1,  size(data_in,3));
        % cat all the data together into a 1d array for each variable
        sess_lme = [sess_1d, sess_1d];
        subject_lme = [subject_1d, subject_1d];
        power_lme = [contra_1d, ctrl_1d];
        condition_lme = [contra_id_1d, ctrl_id_1d];
        
        % double checking.  Error bars overlap for NAc con-ctr but LME says
        % sig diff. this is because we have random effects and a mean +/-SEM plot
        % will not take that into acount
%         this_SEM_contra = nanstd(contra_1d)./sqrt(length(contra_1d(~isnan(contra_1d))));
%         this_SEM_ctrl = nanstd(ctrl_1d)./sqrt(length(ctrl_1d(~isnan(ctrl_1d))));

        %% build the LME for this site
        D_power.tbl = table(subject_lme', sess_lme', power_lme', condition_lme','VariableNames',{'SubjectID','SessID', 'Power', 'Condition'});
        D_power.tbl.SubjectID = nominal(D_power.tbl.SubjectID);
        D_power.tbl.SessID = nominal(D_power.tbl.SessID);
        D_power.tbl.Condition = nominal(D_power.tbl.Condition);
        
        % disp([cfg.row_names{iSite} '___________________'])
        D_power.lme = fitlme(D_power.tbl,'Power~1+Condition+(1|SubjectID)+(1|SessID)');
        % D_power.lme
        % anova(D_power.lme,'DFMethod','satterthwaite')
        
        % collect values
        Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).Est =  D_power.lme.Coefficients.Estimate(2);
        Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).SE =  D_power.lme.Coefficients.SE(2);
        Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val =  D_power.lme.Coefficients.pValue(2);
        Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).Lower =  D_power.lme.Coefficients.Lower(2);
        Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).Upper =  D_power.lme.Coefficients.Upper(2);
        Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).df =  D_power.lme.Coefficients.DF(2);

        Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).tstat =  D_power.lme.Coefficients.tStat(2);

        % hold the P value for plotting later
        p_con_ctr(iSite) = Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val;
        
        
        fprintf(cfg.stats_dir,['\n' PARAMS.all_sites{iSite} ' & Contra-Control  ']);
        if Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val >= 0.05
            p_val_string = ['p$=$ ' num2str(Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val,'%10.2e')];
            h_con_ctr(iSite) = 0; % use for assigning markers later.
            
        elseif (0.0499999 > Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val >= 0.01);
            p_val_string = '*';
            h_con_ctr(iSite) = 1;
            
        elseif (0.0099999 > Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val >= 0.001);
            p_val_string = '**';
            h_con_ctr(iSite) = 1;
            
        elseif 0.00099999 > Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).P_val
            p_val_string = '***';
            h_con_ctr(iSite) = 1;
            
        end
        
        fprintf(cfg.stats_dir,'& %4.2f & (%2.0f) %4.2f %s & %4.2f & %4.2f',...
            Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).Est,...
            Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).df,...
            Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).tstat,...
            p_val_string,...
            Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).Lower,...
            Stats_out.contra_ctrl.(PARAMS.all_sites{iSite}).Upper); % p $<$ 0.001
        
    end
else
    error('Statistical method not specified in cfg.stats_method')
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
if exist(cfg.save_dir, 'dir'); mkdir(cfg.save_dir); end % check is the stats dir exists, if not make the folder. 
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

out.(cfg.method) =avg_vals; 
out.SEM = SEM_vals; 
out.LME = Stats_out; 


end



