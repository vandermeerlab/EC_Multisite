function out = MS_get_phase_distance_LMEs(cfg_in, all_Naris)
%% MS_get_naris_dist: cycles through expkeys or specified table to get the
% distance to the nearest piriform layer in the coronal plane.  This is
% then used to determine the amount of gamma suppresion from the contra to
% the ipsi condition.
%
%
%
%   INPUTS
%      - cfg_in [struct] : contains all configuration paramters

%      - all_Naris [struct] output from Master_Multisite_postprocess
%
%
%   * currently uses internal list of distance instead of those from
%      ExpKeys.  To be added at some point.
%% internal list of distance per site, per subject.  Should be replaced with something computed from the ExpKeys in the future

% pl distances (not used, swapped to by subjects)
Pl_dist = repmat([NaN,    4.841  NaN   NaN  4.327  3.976  NaN],4,1)';
IL_dist = repmat([3.124   NaN    NaN   NaN  NaN    3.329    3.561],4,1)';
OFC_dist = repmat([NaN    1.887  1.265 1    0.825  1.166  0.894],4,1)';
NAc_dist = repmat([0.447  1.897  0.6   1    1.414  0.6    1.4],4,1)';
CG_dist = repmat([NaN     6.251  NaN   NaN  5.855  5.492  6.030],4,1)';


% if strcmp(cfg.traget, 'OB')
R102 = repmat([NaN 1.6 NaN 3.24 NaN],4,1)';
R104 = repmat([ 1.6 NaN 2.16 3.24 4.20],4,1)';
R122 = repmat([NaN NaN 2.16 3.24 NaN],4,1)';
R123 = repmat([NaN NaN NaN 4.20  NaN],4,1)';
R107 = repmat([1.6 NaN 2.16 5.2800 3.24],4,1)';
R108 = repmat([1.6 1.6 2.64 3.24 3.24],4,1)';
R112 = repmat([1.6 1.6 2.16 NaN 3.24],4,1)';
distance_ob = cat(3,R102, R104, R107, R108, R112, R122, R123);

% elseif strcmp(cfg.traget, 'PC')
R102 = repmat([NaN 3.124 NaN 0.447 NaN],4,1)';
R104 = repmat([ 4.841 NaN 1.887 1.897 6.251],4,1)';
R122 = repmat([NaN NaN 1.265 0.6 NaN],4,1)';
R123 = repmat([NaN NaN 1 1.077  NaN],4,1)';
R107 = repmat([4.327 NaN 0.825 1.414 5.855],4,1)';
R108 = repmat([3.967 3.329 1.166 0.6 5.492],4,1)';
R112 = repmat([4.1 3.561 0.894 NaN 6.030],4,1)';
% end

dist_labels = {'PL', 'IL','OFC', 'NAc', 'CG'};
distance_pc = cat(3,R102, R104, R107, R108, R112, R122, R123);
%  distance = cat(3, Pl_dist, IL_dist, OFC_dist, NAc_dist, CG_dist);


%% setup configuration
global PARAMS

cfg_def = [];
cfg_def.pot_trk = {'pot'};
cfg_def.type = 'both'; % whether to output the 'standard' or "white" filtered PSD
cfg_def.plot_type = 'raw';
cfg_def.pot_trk = '';
cfg_def.linewidth = 4;
cfg_def.color.blue = double([158,202,225])/255;
cfg_def.color.green = double([168,221,181])/255;
cfg_def.filter = [45 65; 70 90];
cfg = ProcessConfig(cfg_def, cfg_in);

%% collect all sessions/subjects
if isempty(cfg.pot_trk)
    rec_type = {'pot', 'trk'};
else
    rec_type = cfg.pot_trk;
end

% sites = {'PL'    'IL'    'OFC'    'Piri_O'    'NAc'    'Piri_N'    'CG'};

%% collect Coh and AMP values
subjects = fieldnames(all_Naris);
% set up matrices
for iPhase = 1:length(PARAMS.Phases)
    out_matrix.(PARAMS.Phases{iPhase}) = [];
end

used_pairs = [];
comp_low_coh_array = []; % for low gamma coherence
comp_high_coh_array = [];% for high gamma coherence
comp_low_amp_array = [];% for low gamma amplitude coherence
comp_high_amp_array = [];% for high gamma amplitude coherence

out_matrix.contrast = [];
labels = {'Distance', 'Subjects', 'Session', 'Site', 'Amp_low', 'Amp_high', 'Coh_low', 'Coh_high'}; % labels for the colums of the matrices

count = 0;
for iSub = 1:length(subjects) % loop subjects
    sess_list = fieldnames(all_Naris.(subjects{iSub}));
    for iSess = 1:length(sess_list) % loop sessions
        these_pairs = fieldnames(all_Naris.(subjects{iSub}).(sess_list{iSess}).amp.ac);
        count = count+1;
        for iPair = 1:length(these_pairs) % loop for the pairs of electrodes in this subject/session
            % find the pair idx
            temp_idx = strfind(PARAMS.all_pairs,these_pairs{iPair});
            pair_idx = find(not(cellfun('isempty',temp_idx)));
            this_pair = PARAMS.all_pairs{pair_idx};
            
            if strfind(this_pair, 'Piri')
                comp_low_coh_array(1:4,pair_idx,count) = 0;
                comp_high_coh_array(1:4,pair_idx,count) = 0;
                comp_low_amp_array(1:4,pair_idx,count) = 0;
                comp_high_amp_array(1:4,pair_idx,count) = 0;
                
                continue
            else
                
                for iPhase = 1:length(PARAMS.Phases)
                    % get the mean amp
                    this_amp_F = all_Naris.(subjects{iSub}).(sess_list{iSess}).amp.f.(this_pair).(PARAMS.Phases{iPhase});
                    this_amp = all_Naris.(subjects{iSub}).(sess_list{iSess}).amp.ac.(this_pair).(PARAMS.Phases{iPhase});
                    % low amp
                    this_low_amp = nanmean(this_amp(nearest_idx(cfg.filter(1,1), this_amp_F):nearest_idx(cfg.filter(1,2), this_amp_F))); % get the mean in the low gamma band
                    % high amp
                    this_high_amp = nanmean(this_amp(nearest_idx(cfg.filter(2,1), this_amp_F):nearest_idx(cfg.filter(2,2), this_amp_F))); % get the mean in the low gamma band
                    
                    % get the mean coh
                    this_coh_F = all_Naris.(subjects{iSub}).(sess_list{iSess}).coh.fxx.(this_pair).(PARAMS.Phases{iPhase});
                    this_coh = all_Naris.(subjects{iSub}).(sess_list{iSess}).coh.cxx.(this_pair).(PARAMS.Phases{iPhase});
                    % low amp
                    this_low_coh = nanmean(this_coh(nearest_idx(cfg.filter(1,1), this_coh_F):nearest_idx(cfg.filter(1,2), this_coh_F))); % get the mean in the low gamma band
                    % high amp
                    this_high_coh = nanmean(this_coh(nearest_idx(cfg.filter(2,1), this_coh_F):nearest_idx(cfg.filter(2,2), this_coh_F))); % get the mean in the low gamma band
                    
                    %get the distance from the Piri for the furthest site
                    sites = strsplit(this_pair, '_');
                    S1_dist_idx = find(not(cellfun('isempty',strfind(dist_labels, sites{1}))));
                    S2_dist_idx = find(not(cellfun('isempty',strfind(dist_labels, sites{2}))));
                    
                    
%                          % distance to PC using mean distance
%                         dist_val = mean([distance_pc(S1_dist_idx, 1,iSub),distance_pc(S1_dist_idx, 1,iSub)]);
%                                   % distance to PC using mean distance
%                         dist_val = sum([distance_pc(S1_dist_idx, 1,iSub),distance_pc(S1_dist_idx, 1,iSub)]);
%                   
%                     
%                     % distance to PC using max val 
                    if distance_pc(S1_dist_idx, 1,iSub) > distance_pc(S2_dist_idx, 1,iSub)
                        dist_val = distance_pc(S1_dist_idx, 1,iSub);
                    elseif distance_pc(S1_dist_idx, 1,iSub) < distance_pc(S2_dist_idx, 1,iSub)
                        dist_val = distance_pc(S2_dist_idx, 1,iSub);
                    else
                        dist_val = distance_pc(S1_dist_idx, 1,iSub);
                    end
                    
%                         % distance to PC using min val 
%                     if distance_pc(S1_dist_idx, 1,iSub) > distance_pc(S2_dist_idx, 1,iSub)
%                         dist_val = distance_pc(S2_dist_idx, 1,iSub);
%                     elseif distance_pc(S1_dist_idx, 1,iSub) < distance_pc(S2_dist_idx, 1,iSub)
%                         dist_val = distance_pc(S1_dist_idx, 1,iSub);
%                     else
%                         dist_val = distance_pc(S2_dist_idx, 1,iSub);
%                     end
                    
                    % hold the segment values for each case
                    temp.(PARAMS.Phases{iPhase}).low_coh = this_low_coh;
                    temp.(PARAMS.Phases{iPhase}).high_coh = this_high_coh;
                    temp.(PARAMS.Phases{iPhase}).low_amp = this_low_amp;
                    temp.(PARAMS.Phases{iPhase}).high_amp = this_high_amp;
                    
                    comp_low_coh_array(iPhase,pair_idx,count) =  temp.(PARAMS.Phases{iPhase}).low_coh;
                    comp_high_coh_array(iPhase,pair_idx,count) =  temp.(PARAMS.Phases{iPhase}).high_coh;
                    comp_low_amp_array(iPhase,pair_idx,count) =  temp.(PARAMS.Phases{iPhase}).low_amp;
                    comp_high_amp_array(iPhase,pair_idx,count) =  temp.(PARAMS.Phases{iPhase}).high_coh;
                    
                    
                    out_matrix.(PARAMS.Phases{iPhase}) = cat(1, out_matrix.(PARAMS.Phases{iPhase}),[dist_val, iSub, iSess, pair_idx, this_low_amp, this_high_amp, this_low_coh, this_high_coh]);
                end %  phases
                
                % get the contrast between the contra and ipsi segments.
                out_matrix.contrast = cat(1, out_matrix.contrast,[dist_val, iSub, iSess, pair_idx,...
                    (temp.contra.low_amp - temp.ipsi.low_amp)./(temp.contra.low_amp + temp.ipsi.low_amp),...
                    (temp.contra.high_amp - temp.ipsi.high_amp)./(temp.contra.high_amp + temp.ipsi.high_amp),...
                    (temp.contra.low_coh - temp.ipsi.low_coh)./(temp.contra.low_coh + temp.ipsi.low_coh),...
                    (temp.contra.high_coh - temp.ipsi.high_coh)./(temp.contra.high_coh + temp.ipsi.high_coh)]);
                
%                                 out_matrix.contrast = cat(1, out_matrix.contrast,[dist_val, iSub, iSess, pair_idx,...
%                     (temp.contra.low_amp - temp.ipsi.low_amp),...
%                     (temp.contra.high_amp - temp.ipsi.high_amp),...
%                     (temp.contra.low_coh - temp.ipsi.low_coh),...
%                     (temp.contra.high_coh - temp.ipsi.high_coh)]);
                
                % hold values for comparisons later on
                comp_low_coh_array(5,pair_idx,count) =  nanmean([temp.pre.low_coh, temp.post.low_coh]);
                comp_high_coh_array(5,pair_idx,count) =  nanmean([temp.pre.low_coh, temp.post.high_coh]);
                comp_low_amp_array(5,pair_idx,count) =  nanmean([temp.pre.low_coh, temp.post.low_amp]);
                comp_high_amp_array(5,pair_idx,count) =  nanmean([temp.pre.low_coh, temp.post.high_amp]);
                
                
                used_pairs = [used_pairs, pair_idx];
            end % skip if has Piriform
        end % pairs
    end % sessions
end % subjects


%% get indiviual comparisons between ipsi contra contr
comp_low_coh_array(comp_low_coh_array==0) = NaN;
comp_high_coh_array(comp_high_coh_array==0) = NaN;
comp_low_amp_array(comp_low_amp_array==0) = NaN;
comp_high_amp_array(comp_high_amp_array==0) = NaN;

%remove non-existant pairs of electrodes
unused_pairs = ~ismember(1:length(PARAMS.all_pairs), unique(used_pairs)); % remove pairs without any data;
comp_low_coh_array(:,unused_pairs, :) = [];
comp_high_coh_array(:,unused_pairs, :) = [];
comp_low_amp_array(:,unused_pairs, :) = [];
comp_high_amp_array(:,unused_pairs, :) = [];

stats_file = fopen([PARAMS.stats_dir 'COH_low_LME_stats_2020.txt'], 'w');


cfg_stats = [];
cfg_stats.NaN_correct = 1; % workaround for
cfg_stats.title = 'Coh low gamma';
cfg_stats.method= 'median';
cfg_stats.row_names= PARAMS.all_pairs(unique(used_pairs));
cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
cfg_stats.s_idx= 1:length(cfg_stats.row_names);
cfg_stats.ft_size= 20;
cfg_stats.stats_method = 'lme';
cfg_stats.save_dir= [PARAMS.inter_dir 'Phase_Stats_low_coh_2020'];
cfg_stats.stats_dir = stats_file;
MS_stats_Pairs(cfg_stats,comp_low_coh_array)

close all
% same thing for high coh
cfg_stats.title = 'Coh high gamma';
cfg_stats.save_dir= [PARAMS.inter_dir 'Phase_Stats_high_coh_2020'];
stats_file = fopen([PARAMS.stats_dir 'COH_high_LME_stats_2020.txt'], 'w');
cfg_stats.stats_dir = stats_file;
MS_stats_Pairs(cfg_stats,comp_high_coh_array)
close all

% same thing for low amp
cfg_stats.title = 'Amp low gamma';
cfg_stats.save_dir= [PARAMS.inter_dir 'Phase_Stats_low_amp_2020'];
stats_file = fopen([PARAMS.stats_dir 'AMP_low_LME_stats_2020.txt'], 'w');
cfg_stats.stats_dir = stats_file;
MS_stats_Pairs(cfg_stats,comp_low_amp_array)
close all

% same thing for high amp
cfg_stats.title = 'Amp high gamma';
cfg_stats.save_dir= [PARAMS.inter_dir 'Phase_Stats_high_amp_2020'];
stats_file = fopen([PARAMS.stats_dir 'AMP_high_LME_stats_2020.txt'], 'w');
cfg_stats.stats_dir = stats_file;
MS_stats_Pairs(cfg_stats,comp_high_amp_array)
close all

%% match distance to each subject/site/session for the difference between contra and ipsi power
Subjects = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7'};
sites = {'PL', 'IL', 'OFC', 'NAc', 'CG'};
pairs = unique(out_matrix.contra(:,4));
%% apply the distance to OB in a corresponding array.
c_ord = linspecer(length(pairs));
m_ord = {'o', '+', '*', 'x', 's', 'd', 'p'};
figure(100)
hold on
for iSub = 1:length(Subjects)
    for iSite = 1:length(pairs)
        p_idx = find(out_matrix.contra(:,4) == pairs(iSite) & out_matrix.contra(:,2) == iSub);
%         plot(distance_ob(iSite,:,iSub), out_matrix.contra(:,7),m_ord{iSub})
        plot(out_matrix.contra(p_idx,1), out_matrix.contra(p_idx,5),m_ord{iSub},'MarkerEdgeColor', c_ord(iSite,:), 'markersize', 10)
    end
end
xlabel('Distance from PC (mm)')
ylabel('Ipsi/Contra contrast index')
% annoying forced legend.  Avoids issue of markers and colors not working properly.
hold on
h = zeros(length(sites), 1);
for iSite = 1:length(pairs)
    h(iSite) = plot(NaN,NaN,'color', c_ord(iSite,:));
end
[~, hobj, ~, ~] = legend(h, PARAMS.all_pairs(pairs), 'location', 'northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',3);
% legend(sites, 'location', 'southeast')
SetFigure([], gcf)


%% make a talbe

clear D_coh
D_coh.tbl = table(out_matrix.contrast(:,1),out_matrix.contrast(:,2),out_matrix.contrast(:,3),out_matrix.contrast(:,4),out_matrix.contrast(:,5),...
    out_matrix.contrast(:,6),out_matrix.contrast(:,7),out_matrix.contrast(:,8),...
    'VariableNames',{'Distance_pc', 'RatID','SessID', 'Pair','Amp_low', 'Amp_high','Coh_low', 'Coh_high'});
D_coh.tbl.RatID = nominal(D_coh.tbl.RatID);
D_coh.tbl.SessID = nominal(D_coh.tbl.SessID);
D_coh.tbl.Pair = nominal(D_coh.tbl.Pair);
 %% basic comparisons between ipsi and contra for each pair.
% 
% D_coh.lme = fitlme(D_coh.tbl,'Coh_low~1+Condition+(1|SubjectID)+(1|SessID)');
% % D_coh.lme
% % anova(D_coh.lme,'DFMethod','satterthwaite')
% 
% % collect values
% Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Est =  D_coh.lme.Coefficients.Estimate(2);
% Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).SE =  D_coh.lme.Coefficients.SE(2);
% Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val =  D_coh.lme.Coefficients.pValue(2);
% Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Lower =  D_coh.lme.Coefficients.Lower(2);
% Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).Upper =  D_coh.lme.Coefficients.Upper(2);
% % T-stats
% Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).tstat =  D_coh.lme.Coefficients.tStat(2);
% 
% % hold the P value for plotting later
% p_ip_con(iSite) = Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val;
% 
% 
% fprintf(cfg.stats_dir,['\n' PARAMS.all_sites{iSite} '& Ipsi-Contra  ']);
% if Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val >= 0.05
%     fprintf(cfg.stats_dir,' & %4.2f &  p = %4.2f ',Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).tstat,...
%         Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val);
%     h_ip_con(iSite) = 0; % use for assigning markers later.
%     
% elseif (0.049 > Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val >= 0.01);
%     fprintf(cfg.stats_dir,' & %4.2f &  p $<$ 0.05',Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).tstat);
%     h_ip_con(iSite) = 1;
%     
% elseif (0.009 > Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val) && (Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val >= 0.001);
%     fprintf(cfg.stats_dir,' & %4.2f &  p $<$ 0.01',Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).tstat);
%     h_ip_con(iSite) = 1;
%     
% elseif 0.0009 > Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).P_val
%     fprintf(cfg.stats_dir,' & %4.2f &  p $<$ 0.001 ',Stats_out.ipsi_contra.(PARAMS.all_sites{iSite}).tstat);
%     h_ip_con(iSite) = 1;
%     
% end



%% run the distance LME

% make some models
% Amp low
D_coh.lme_baseline_amp_low = fitlme(D_coh.tbl,'Amp_low~1+(1|RatID)+(1|SessID)');
D_coh.lme_amp_low = fitlme(D_coh.tbl,'Amp_low~1+Distance_pc+(1|RatID)+(1|SessID)');
%Amp high
D_coh.lme_baseline_amp_high = fitlme(D_coh.tbl,'Amp_high~1+(1|RatID)+(1|SessID)');
D_coh.lme_amp_high = fitlme(D_coh.tbl,'Amp_high~1+Distance_pc+(1|RatID)+(1|SessID)');
% coh low
D_coh.lme_baseline_coh_low = fitlme(D_coh.tbl,'Coh_low~1+(1|RatID)+(1|SessID)');
D_coh.lme_coh_low = fitlme(D_coh.tbl,'Coh_low~1+Distance_pc+(1|RatID)+(1|SessID)');
% coh high
D_coh.lme_baseline_coh_high = fitlme(D_coh.tbl,'Coh_high~1+(1|RatID)+(1|SessID)');
D_coh.lme_coh_high = fitlme(D_coh.tbl,'Coh_high~1+Distance_pc+(1|RatID)+(1|SessID)');

% comparisons
%amp low
D_coh.comp_amp_low = compare(D_coh.lme_baseline_amp_low,D_coh.lme_amp_low);

D_coh.comp_amp_high = compare(D_coh.lme_baseline_amp_high,D_coh.lme_amp_high);

D_coh.comp_coh_low = compare(D_coh.lme_baseline_coh_low,D_coh.lme_coh_low);

D_coh.comp_coh_high = compare(D_coh.lme_baseline_coh_high,D_coh.lme_coh_high);


%% write the output
if exist(['LME_phase_Naris' datestr(date, 'YY_mm_dd') '2020.txt'], 'file');
    delete(['LME_phase_Naris' datestr(date, 'YY_mm_dd') '2020.txt'])
end
clc
diary('on')
diary(['LME_phase_Naris' datestr(date, 'YY_mm_dd') '2020.txt'])
disp(' LME amp_low' )
disp(D_coh.lme_amp_low)
disp(' LME Anova Out amp_low')
anova(D_coh.lme_amp_low)
disp('Compare amp_low to baseline')
compare(D_coh.lme_baseline_amp_low,D_coh.lme_amp_low)
disp('     ')
disp(' LME amp_high' )
disp(D_coh.lme_amp_high)
disp(' LME Anova Out amp_high')
anova(D_coh.lme_amp_high)
disp('Compare amp_high to baseline')
compare(D_coh.lme_baseline_amp_high,D_coh.lme_amp_high)
% coherence
disp('     ')
disp(' LME coh_low' )
disp(D_coh.lme_coh_low)
disp(' LME Anova Out coh_low')
anova(D_coh.lme_coh_low)
disp('Compare coh_low to baseline')
compare(D_coh.lme_baseline_coh_low,D_coh.lme_coh_low)
% coh high
disp(' LME coh_high' )
disp(D_coh.lme_coh_high)
disp(' LME Anova Out coh_high')
anova(D_coh.lme_coh_high)
disp('Compare coh_high to baseline')
compare(D_coh.lme_baseline_coh_high,D_coh.lme_coh_high)
diary('off')
movefile(['LME_phase_Naris' datestr(date, 'YY_mm_dd') '2020.txt'], PARAMS.stats_dir);
% %% try it as a logistic for 'prox' vs 'dist'.  Did not use.
% % this didn't work.
% clear L_power
% % add new value for distances greater than 2mm or less than
% log_1d = cell(size(dist_1d));
% prox_idx = dist_1d <=2;
% for ii = length(log_1d):-1:1
%     if prox_idx(ii) ==1
%         log_1d{ii} = 'prox';
%     else
%         log_1d{ii} = 'dist';
%     end
% end
end
% %% odd attempt at glm
% L_power.tbl = table(rat_1d, sess_1d, prox_idx, pow_1d,'VariableNames',{'RatID','SessID', 'Distance', 'Power'});
% L_power.tbl.RatID = nominal(L_power.tbl.RatID);
% L_power.tbl.SessID = nominal(L_power.tbl.SessID);
% L_power.tbl.Distance = logical(L_power.tbl.Distance);
%
% % m_spec = 'Power ~ 1+ Distance +(1|RatID) + (1|SessID)';
% % m_spec = 'Distance ~ 1+ Power +(1|RatID) + (1|SessID)';
% m_spec = 'Distance ~ Power ';
%
% glm_out = fitglme(L_power.tbl, m_spec, 'distribution', 'binomial')
% % fitglm(L_power.tbl, m_spec)
%
% plotResiduals(glm_out_2,'fitted')


