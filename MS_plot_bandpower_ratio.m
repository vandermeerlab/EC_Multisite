function MS_plot_bandpower_ratio(cfg_in, Naris_in)
%% MS_plot_bandpower_ratio: Uses the bandpower measure to plot the ratios between naris conditions

%% set up defaults

cfg_def = [];
cfg_def.type = 'both'; % whether to output the 'standard' or "white" filtered PSD
cfg_def.linewidth = 2;
cfg_def.plot_type = 'raw';
cfg_def.ylims = [-75 75];
cfg_def.ylims_norm = [0 3];
cfg_def.pot_trk = 'pot';
cfg = ProcessConfig2(cfg_def, cfg_in);
global PARAMS
c_ord = [linspecer(length(PARAMS.Phases)); [.6 .6 .6]];

%% Collect all the power ratios across the different sites per subject/session
stats_file = fopen([PARAMS.stats_dir 'Bandpower_stats_2020.txt'], 'w');
% stats_file_high = fopen([PARAMS.stats_dir 'POW_stats_high_2020.txt'], 'w');


% set up some 
Phases = [PARAMS.Phases, 'control'];
all_POW_low = [];
all_POW_high = [];
all_POW_con_low = [];
all_POW_con_high = [];
norm_all_POW_low = [];
norm_all_POW_high = [];
norm_all_POW_con_low = [];
norm_all_POW_con_high = [];

subjects = fieldnames(Naris_in);

% loop subjects
for iSub = 1:length(subjects)
    sess_list = fieldnames(Naris_in.(subjects{iSub}));
    sites = {'PL'    'IL'    'OFC'    'Piri_O'    'NAc'    'Piri_N'    'CG'};
    % create the empty array
    POW_low = NaN(length(Phases),length(sites),4);
    POW_high = NaN(length(Phases),length(sites),4);
    POW_con_low = NaN(length(Phases),length(sites),4);
    POW_con_high = NaN(length(Phases),length(sites),4);
    norm_POW_low = NaN(length(Phases),length(sites),4);
    norm_POW_high = NaN(length(Phases),length(sites),4);
    norm_POW_con_low = NaN(length(Phases),length(sites),4);
    norm_POW_con_high = NaN(length(Phases),length(sites),4);
    
    % loop sessions
    for iSess = 1:length(sess_list)
        % loop sites
        
        for iPhase = 1:length(Phases)
            these_sites = fieldnames(Naris_in.(subjects{iSub}).(sess_list{iSess}).(Phases{iPhase})); % get the sites list for this subject/session
            
            for iSite = 1:length(sites)
                site_idx = strcmp(these_sites, [sites{iSite} '_' cfg.pot_trk]); % find the indices of this site.
                if sum(site_idx) ==1
                    site_idx = find(site_idx ==1);
                    POW_low(iPhase,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).(Phases{iPhase}).([sites{iSite} '_' cfg.pot_trk]).bandpow.low; 
                    POW_high(iPhase,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).(Phases{iPhase}).([sites{iSite} '_' cfg.pot_trk]).bandpow.high; 
                    
                    POW_cont_low(iPhase,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).(Phases{iPhase}).([sites{iSite} '_' cfg.pot_trk]).bandpow.cont_low;
                    POW_con_high(iPhase,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).(Phases{iPhase}).([sites{iSite} '_' cfg.pot_trk]).bandpow.cont_high; 
                    
                    % normalize to the control condition
                    norm_POW_low(:,iSite,iSess) = POW_low(:,iSite,iSess)./POW_low(5,iSite,iSess);
                    norm_POW_high(:,iSite,iSess) = POW_high(:,iSite,iSess)./POW_high(5,iSite,iSess);
                    norm_POW_con_low(:,iSite,iSess) = POW_con_low(:,iSite,iSess)./POW_con_low(5,iSite,iSess);
                    norm_POW_con_high(:,iSite,iSess) = POW_con_high(:,iSite,iSess)./POW_con_high(5,iSite,iSess);
                    
                end
            end
        end
    end
    all_POW_low = cat(3,all_POW_low, POW_low);
    all_POW_high = cat(3,all_POW_high, POW_high);
    all_POW_con_low = cat(3,all_POW_con_low, POW_con_low);
    all_POW_con_high = cat(3,all_POW_con_high, POW_con_high);

    norm_all_POW_low = cat(3,norm_all_POW_low, norm_POW_low);
    norm_all_POW_high = cat(3,norm_all_POW_high, norm_POW_high);
    norm_all_POW_con_low = cat(3,norm_all_POW_con_low, norm_POW_con_low);
    norm_all_POW_con_high = cat(3,norm_all_POW_con_high, norm_POW_con_high);
    
    %collect for each subject
    Rats.(subjects{iSub}).all_POW_low = POW_low;
    Rats.(subjects{iSub}).all_POW_high = POW_high;
    Rats.(subjects{iSub}).all_POW_con_low = POW_con_low;
    Rats.(subjects{iSub}).all_POW_con_high = POW_con_high;
    
    Rats.(subjects{iSub}).norm_all_POW_low = norm_POW_low;
    Rats.(subjects{iSub}).norm_all_POW_high = norm_POW_high;
    Rats.(subjects{iSub}).norm_all_POW_con_low = norm_POW_con_low;
    Rats.(subjects{iSub}).norm_all_POW_con_high = norm_POW_con_high;
end



%% stats
Exp= {'Four', 'Piri'};
for iExp = 1:length(Exp)
    if strcmp(Exp{iExp}, 'Four')
        s_idx = [1 2 3 5 7];
    elseif strcmp(Exp{iExp}, 'Piri')
        s_idx = [3 4 5 6];
    end
    
    cfg_stats = [];
    cfg_stats.title = strcat({'POW low gamma'},{' '}, Exp{iExp});
    cfg_stats.title = cfg_stats.title{1};
    cfg_stats.method= 'median';
    cfg_stats.stats_method = cfg.stats_method;
    cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
    cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
    cfg_stats.s_idx= s_idx;
    cfg_stats.ft_size= 20;
    cfg_stats.save_dir= [PARAMS.inter_dir 'POW_2020'];
    cfg_stats.stats_dir = stats_file;
    stats_out.(Exp{iExp}).low = MS_stats(cfg_stats, all_POW_low);
    
    close all
    
    % high gamma median using MS sites
    cfg_stats = [];
    cfg_stats.title = strcat({'POW high gamma'},{' '},Exp{iExp});
    cfg_stats.title = cfg_stats.title{1};
    cfg_stats.method= 'median';
    cfg_stats.stats_method = cfg.stats_method;
    cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
    cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
    cfg_stats.s_idx= s_idx;
    cfg_stats.ft_size= 20;
    cfg_stats.save_dir= [PARAMS.inter_dir 'POW_2020'];
    cfg_stats.stats_dir = stats_file;
    stats_out.(Exp{iExp}).high = MS_stats(cfg_stats, all_POW_high);
    close all
    
    % same but normalized
    cfg_stats = [];
    cfg_stats.title = strcat({'POW low gamma Normalized'},{' '}, Exp{iExp});
    cfg_stats.title = cfg_stats.title{1};
    cfg_stats.method= 'median';
    cfg_stats.stats_method = cfg.stats_method;
    cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
    cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
    cfg_stats.s_idx= s_idx;
    cfg_stats.ft_size= 20;
    cfg_stats.save_dir= [PARAMS.inter_dir 'POW_2020'];
    cfg_stats.stats_dir = stats_file;
    stats_out.(Exp{iExp}).low = MS_stats(cfg_stats, norm_all_POW_low);
    
    close all
    
    % high gamma median using MS sites
    cfg_stats = [];
    cfg_stats.title = strcat({'POW high gamma Normalized'},{' '},Exp{iExp});
    cfg_stats.title = cfg_stats.title{1};
    cfg_stats.method= 'median';
    cfg_stats.stats_method = cfg.stats_method;
    cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
    cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
    cfg_stats.s_idx= s_idx;
    cfg_stats.ft_size= 20;
    cfg_stats.save_dir= [PARAMS.inter_dir 'POW_2020'];
    cfg_stats.stats_dir = stats_file;
    stats_out.(Exp{iExp}).high = MS_stats(cfg_stats, all_POW_high);
    close all
    
    
end



%% descriptive stats
fid = fopen([PARAMS.stats_dir 'POW_descriptive2020.txt'], 'w');

bands = {'low', 'high'};
phases = {'pre'  'ipsi'  'contra'  'post', 'control'};

fprintf(fid, ['**************** ' date ' ****************\n']);

fprintf(fid, ['\nAUC using White_Pxx and ' iComp{1} '\n']);
for iBand= 1:length(bands)
    this_pow = [];
    if strcmp(bands{iBand}, 'low')
        this_pow = all_POW_low.White_Pxx;
    elseif strcmp(bands{iBand}, 'high')
        this_pow = all_POW_high.White_Pxx;
    end
    
    
    fprintf(fid,['\n-------- ' bands{iBand} '---------\n']);
    for iSite = 1:length(sites)
        n_spaces = 6 - length(sites{iSite});
        fprintf(fid,[sites{iSite} ':%s' ], repmat(' ', 1,n_spaces));
        fprintf(fid,repmat('\b', 1, length(sites{iSite})));
        for iPhase = 1:size(this_pow,1)
            
            %             all_AUC.(bands{iBand})(iSite, iPhase) = nanmedian(this_pow(iPhase, iSite,:));
            %             all_AUC_std.(bands{iBand})(iSite, iPhase) = nanstd(this_pow(iPhase, iSite,:))./sqrt(size(this_pow(iPhase, iSite,:),3));
            %             fprintf(fid,[phases{iPhase} ' median= %.2f +/- %.2f  '], all_AUC.(bands{iBand})(iSite, iPhase), all_AUC_std.(bands{iBand})(iSite,iPhase));
            these_vals = this_pow(iPhase, iSite,:);
            SEM = nanstd(these_vals)/   sqrt(length(these_vals(~isnan(these_vals))));
            fprintf(fid,[phases{iPhase} ' median= %.2f +/- %.2f  '], nanmedian(these_vals), SEM);
            
        end
        fprintf(fid,'\n');
    end
end


fclose(fid);
%% make a legend for all the plots.
figure(9999)
c_ord = linspecer(3);
leg_val = {'control', 'ipsi', 'contra'};
b = bar(magic(3), 'BaseValue', 1);

for iPhase = 1:3
    set(b(iPhase), 'FaceColor', c_ord(iPhase,:), 'visible', 'off')
end
% set(gca, 'xticklabel', sites(1:3), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
legend(leg_val, 'location', 'south', 'orientation', 'horizontal');
legend boxoff
axis off
cfg_plt1.pos = [600 50 560*1.4 560*1.8];
cfg_plt1.ft_size = 18;
SetFigure(cfg_plt1, gcf)


if isunix
    %     saveas(gcf, [PARAMS.inter_dir '/POW_fit/legend'], 'epsc')
    saveas_eps('legend',[PARAMS.inter_dir '/POW_fit/'])
else
    %     saveas(gcf, [PARAMS.inter_dir '\POW_fit\legend'], 'epsc')
    saveas_eps('legend',[PARAMS.inter_dir '\POW_fit\'])
end
end

