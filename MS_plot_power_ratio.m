function MS_plot_power_ratio(cfg_in, Naris_in)
%% MS_plot_psd: plots multiple power spectral densities for the data files
% in the "Naris" structure (output from MS_collect_psd)
%
% inputs:
%    -cfg_in: [struct] contains configuration paramters
%    -Naris: [struct] contains power and frequency values for each channel
%    for each subject/session/phase
%
%    this script currently uses a global parameter set to determine where
%    to save the output figures

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
stats_file = fopen([PARAMS.stats_dir 'AOC_stats.txt'], 'w');

for iComp = {'Pink', 'Exp2'}
    
    AOC_low = []; AOC_high = [];
    AOC_con_low = []; AOC_con_high = [];
    types = {'Pxx', 'White_Pxx'};
    for iType = 1:length(types)
        all_AOC_low.(types{iType}) = [];
        all_AOC_high.(types{iType}) = [];
        all_AOC_con_low.(types{iType}) = [];
        all_AOC_con_high.(types{iType}) = [];
        norm_all_AOC_low.(types{iType}) = [];
        norm_all_AOC_high.(types{iType}) = [];
        norm_all_AOC_con_low.(types{iType}) = [];
        norm_all_AOC_con_high.(types{iType}) = [];
    end
    subjects = fieldnames(Naris_in);
    for iType = 1:length(types)
        it_log = {};
        for iSub = 1:length(subjects)
            sess_list = fieldnames(Naris_in.(subjects{iSub}));
            sites = {'PL'    'IL'    'OFC'    'Piri_O'    'NAc'    'Piri_N'    'CG'};
            % create the empty array
            AOC_low.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            AOC_high.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            AOC_con_low.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            AOC_con_high.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            norm_AOC_low.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            norm_AOC_high.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            norm_AOC_con_low.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            norm_AOC_con_high.(types{iType}) = NaN(length(PARAMS.Phases)+1, length(sites),4);
            for iSess = 1:length(sess_list);
                %
                for iSite = 1:length(sites)
                    site_idx = strcmp(Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_labels, [sites{iSite} '_' cfg.pot_trk]);
                    if sum(site_idx) ==1
                        site_idx = find(site_idx ==1);
                        AOC_low.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio.(types{iType}).(iComp{1}).low(:,site_idx);
                        AOC_high.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio.(types{iType}).(iComp{1}).high(:,site_idx);
                        con_types = fieldnames(Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}).(iComp{1}));
                        AOC_con_low.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}).(iComp{1}).(con_types{1})(:,site_idx);
                        AOC_con_high.(types{iType})(:,iSite, iSess) = Naris_in.(subjects{iSub}).(sess_list{iSess}).ratio_con.(types{iType}).(iComp{1}).(con_types{1})(:,site_idx);
                        %
                        norm_AOC_low.(types{iType})(:,iSite,iSess) = AOC_low.(types{iType})(:,iSite,iSess)./AOC_low.(types{iType})(5,iSite,iSess);
                        norm_AOC_high.(types{iType})(:,iSite,iSess) = AOC_high.(types{iType})(:,iSite,iSess)./AOC_high.(types{iType})(5,iSite,iSess);
                        norm_AOC_con_low.(types{iType})(:,iSite,iSess) = AOC_con_low.(types{iType})(:,iSite,iSess)./AOC_con_low.(types{iType})(5,iSite,iSess);
                        norm_AOC_con_high.(types{iType})(:,iSite,iSess) = AOC_con_high.(types{iType})(:,iSite,iSess)./AOC_con_high.(types{iType})(5,iSite,iSess);
                        
                    end
                end
                it_log{end+1,1} =  sess_list{iSess};
            end
            all_AOC_low.(types{iType}) = cat(3,all_AOC_low.(types{iType}), AOC_low.(types{iType}));
            all_AOC_high.(types{iType}) = cat(3,all_AOC_high.(types{iType}), AOC_high.(types{iType}));
            all_AOC_con_low.(types{iType}) = cat(3,all_AOC_con_low.(types{iType}), AOC_con_low.(types{iType}));
            all_AOC_con_high.(types{iType}) = cat(3,all_AOC_con_high.(types{iType}), AOC_con_high.(types{iType}));
            % normalize to control in each session
            norm_all_AOC_low.(types{iType}) = cat(3,norm_all_AOC_low.(types{iType}), norm_AOC_low.(types{iType}));
            norm_all_AOC_high.(types{iType}) = cat(3,norm_all_AOC_high.(types{iType}), norm_AOC_high.(types{iType}));
            norm_all_AOC_con_low.(types{iType}) = cat(3,norm_all_AOC_con_low.(types{iType}), norm_AOC_con_low.(types{iType}));
            norm_all_AOC_con_high.(types{iType}) = cat(3,norm_all_AOC_con_high.(types{iType}), norm_AOC_con_high.(types{iType}));
            
            %collect for each subject
            Rats.(subjects{iSub}).all_AOC_low.(types{iType}) = AOC_low.(types{iType});
            Rats.(subjects{iSub}).all_AOC_high.(types{iType}) = AOC_high.(types{iType});
            Rats.(subjects{iSub}).all_AOC_con_low.(types{iType}) = AOC_con_low.(types{iType});
            Rats.(subjects{iSub}).all_AOC_con_high.(types{iType}) = AOC_con_high.(types{iType});
            
            Rats.(subjects{iSub}).norm_all_AOC_low.(types{iType}) = norm_AOC_low.(types{iType});
            Rats.(subjects{iSub}).norm_all_AOC_high.(types{iType}) = norm_AOC_high.(types{iType});
            Rats.(subjects{iSub}).norm_all_AOC_con_low.(types{iType}) = norm_AOC_con_low.(types{iType});
            Rats.(subjects{iSub}).norm_all_AOC_con_high.(types{iType}) = norm_AOC_con_high.(types{iType});
        end
    end
    
    
    %% stats
    Exp= {'Four', 'Piri'};
    for iExp = 1:length(Exp)
        if strcmp(Exp{iExp}, 'Four')
            s_idx = [1 2 3 5 7];
        else strcmp(Exp{iExp}, 'Piri')
            s_idx = [3 4 5 6];
        end
        types = { 'White_Pxx'}; % used to use 'Pxx' as well but it was not great
        for iType = 1:length(types);   % low gamma median using MS sites
            cfg_stats = [];
            cfg_stats.title = strcat({'AOC low gamma'},{' '}, Exp{iExp},{' '},types{iType},{' '},iComp);
            cfg_stats.title = cfg_stats.title{1};
            cfg_stats.method= 'median';
            cfg_stats.stats_method = cfg.stats_method; 
            cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
            cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
            cfg_stats.s_idx= s_idx;
            cfg_stats.ft_size= 20;
            cfg_stats.save_dir= [PARAMS.inter_dir 'AOC_fit'];
            cfg_stats.stats_dir = stats_file;
            MS_stats(cfg_stats, all_AOC_low.(types{iType}));
            
            close all
            
            % high gamma median using MS sites
            cfg_stats = [];
            cfg_stats.title = strcat({'AOC high gamma'},{' '},Exp{iExp},{' '},types{iType},{' '},iComp);
            cfg_stats.title = cfg_stats.title{1};
            cfg_stats.method= 'median';
            cfg_stats.stats_method = cfg.stats_method;
            cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
            cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
            cfg_stats.s_idx= s_idx;
            cfg_stats.ft_size= 20;
            cfg_stats.save_dir= [PARAMS.inter_dir 'AOC_fit'];
            cfg_stats.stats_dir = stats_file;
            
            MS_stats(cfg_stats, all_AOC_high.(types{iType}));
            close all
             end
    end
end
            
            
%% descriptive stats
fid = fopen([PARAMS.stats_dir 'AUC_descriptive.txt'], 'w');

bands = {'low', 'high'};
phases = {'pre'  'ipsi'  'contra'  'post', 'control'};

fprintf(fid, ['**************** ' date ' ****************\n']);

fprintf(fid, ['\nAUC using White_Pxx and ' iComp{1} '\n']);
for iBand= 1:length(bands)
    this_pow = [];
if strcmp(bands{iBand}, 'low')
    this_pow = all_AOC_low.White_Pxx; 
elseif strcmp(bands{iBand}, 'high')
    this_pow = all_AOC_high.White_Pxx; 
end
    fprintf(fid,['\n-------- ' bands{iBand} '---------\n']);
    for iSite = 1:length(sites)
        n_spaces = 6 - length(sites{iSite});
        fprintf(fid,[sites{iSite} ':%s' ], repmat(' ', 1,n_spaces));
        fprintf(fid,repmat('\b', 1, length(sites{iSite})));
        for iPhase = 1:size(this_pow,1)
            all_AUC.(bands{iBand})(iSite, iPhase) = nanmedian(this_pow(iPhase, iSite,:));
            all_AUC_std.(bands{iBand})(iSite, iPhase) = nanstd(this_pow(iPhase, iSite,:))./sqrt(size(this_pow(iPhase, iSite,:),1));
            fprintf(fid,[phases{iPhase} ' median= %.2f +/- %.2f  '], all_AUC.(bands{iBand})(iSite, iPhase), all_AUC_std.(bands{iBand})(iSite,iPhase)); 
        end
        fprintf(fid,'\n')
    end
end



fclose(fid)
            % same but using normalized
            
%             cfg_stats = [];
%             cfg_stats.title = strcat({'Norm AOC low gamma'},{' '}, Exp{iExp},{' '},types{iType},{' '},iComp);
%             cfg_stats.title = cfg_stats.title{1};
%             cfg_stats.method= 'median';
%             cfg_stats.stats_method = cfg.stats_method;
%             cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
%             cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
%             cfg_stats.s_idx= s_idx;
%             cfg_stats.ft_size= 20;
%             cfg_stats.save_dir= [PARAMS.inter_dir 'AOC_fit'];
%             cfg_stats.stats_dir = stats_file;
%             
%             MS_stats(cfg_stats, norm_all_AOC_low.(types{iType}));
%             
%             close all
%             
%             % high gamma median using MS sites
%             cfg_stats = [];
%             cfg_stats.title = strcat({'Norm AOC high gamma'},{' '}, Exp{iExp},{' '},types{iType},{' '},iComp);
%             cfg_stats.title = cfg_stats.title{1};
%             cfg_stats.method= 'median';
%             cfg_stats.stats_method = cfg.stats_method;
%             cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
%             cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post'};
%             cfg_stats.s_idx= s_idx;
%             cfg_stats.ft_size= 20;
%             cfg_stats.save_dir= [PARAMS.inter_dir 'AOC_fit'];
%             cfg_stats.stats_dir = stats_file;
%             
%             MS_stats(cfg_stats, norm_all_AOC_high.(types{iType}));
%             close all
            

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
    %     saveas(gcf, [PARAMS.inter_dir '/AOC_fit/legend'], 'epsc')
    saveas_eps('legend',[PARAMS.inter_dir '/AOC_fit/'])
else
    %     saveas(gcf, [PARAMS.inter_dir '\AOC_fit\legend'], 'epsc')
    saveas_eps('legend',[PARAMS.inter_dir '\AOC_fit\'])
end
end
%% this is all old stuff before
% %% get the stats
% for iType = 1:length(types)
%     mean_AOC_low.(types{iType}) = nanmedian(all_AOC_low.(types{iType}), 3)';
%     mean_AOC_high.(types{iType}) = nanmedian(all_AOC_high.(types{iType}),3)';
%     mean_AOC_con_low.(types{iType}) =  nanmedian(all_AOC_con_low.(types{iType}),3)';
%     mean_AOC_con_high.(types{iType}) =  nanmedian(all_AOC_con_high.(types{iType}),3)';
%     % try it relative ot the control condition
%     norm_mean_AOC_low.(types{iType}) = nanmedian(norm_all_AOC_low.(types{iType}), 3)';
%     norm_mean_AOC_high.(types{iType}) = nanmedian(norm_all_AOC_high.(types{iType}),3)';
%     norm_mean_AOC_con_low.(types{iType}) =  nanmedian(norm_all_AOC_con_low.(types{iType}),3)';
%     norm_mean_AOC_con_high.(types{iType}) =  nanmedian(norm_all_AOC_con_high.(types{iType}),3)';
%     %individual subjects
%     for iSub = 1:length(subjects)
%         Rats.(subjects{iSub}).all_AOC_low_mean.(types{iType}) = nanmedian( Rats.(subjects{iSub}).all_AOC_low.(types{iType}), 3)';
%         Rats.(subjects{iSub}).all_AOC_high_mean.(types{iType}) = nanmedian( Rats.(subjects{iSub}).all_AOC_high.(types{iType}), 3)';
%         Rats.(subjects{iSub}).all_AOC_low_con_mean.(types{iType}) = nanmedian( Rats.(subjects{iSub}).all_AOC_con_low.(types{iType}), 3)';
%         Rats.(subjects{iSub}).all_AOC_high_con_mean.(types{iType}) = nanmedian( Rats.(subjects{iSub}).all_AOC_con_high.(types{iType}), 3)';
%
%         Rats.(subjects{iSub}).norm_all_AOC_low_mean.(types{iType}) = nanmedian(Rats.(subjects{iSub}).norm_all_AOC_low.(types{iType}), 3)';
%         Rats.(subjects{iSub}).norm_all_AOC_high_mean.(types{iType}) = nanmedian(Rats.(subjects{iSub}).norm_all_AOC_low.(types{iType}), 3)';
%         Rats.(subjects{iSub}).norm_all_AOC_low_con_mean.(types{iType}) = nanmedian(Rats.(subjects{iSub}).norm_all_AOC_low.(types{iType}), 3)';
%         Rats.(subjects{iSub}).norm_all_AOC_high_con_mean.(types{iType}) = nanmedian(Rats.(subjects{iSub}).norm_all_AOC_low.(types{iType}), 3)';
%     end
%
%     %     rm_idx = strfind(sites, 'Piri'); Index = find(not(cellfun('isempty', rm_idx)));
%     %     AOC_low.(types{iType})([4, 6],:) = [];
%     %     AOC_high.(types{iType})([4, 6],:) = [];
%     %     AOC_con_low.(types{iType})([4, 6],:) = [];
%     %     AOC_con_high.(types{iType})([4, 6],:) = [];
%     %     sites{1,6} = [];sites{1,4} = [];
%     %     sites = sites(~cellfun('isempty',sites));
%     %     for iPhase = 1:size(AOC_low.(types{iType}),2)
%     %         norm_low.(types{iType})(:,iPhase) = AOC_low.(types{iType})(:,iPhase) ./ AOC_low.(types{iType})(:,5);
%     %         norm_high.(types{iType})(:,iPhase) = AOC_high.(types{iType})(:,iPhase) ./ AOC_high.(types{iType})(:,5);
%     %         norm_con_low.(types{iType})(:,iPhase) = AOC_con_low.(types{iType})(:,iPhase) ./ AOC_con_low.(types{iType})(:,5);
%     %         norm_con_high.(types{iType})(:,iPhase) = AOC_con_high.(types{iType})(:,iPhase) ./ AOC_con_high.(types{iType})(:,5);
%     %     end
%
%
%     %% get error bars
%     SEM_AOC_low.(types{iType}) = (nanstd(all_AOC_low.(types{iType}),[],3)./sqrt(size(all_AOC_low.(types{iType}),3)))';
%     SEM_AOC_high.(types{iType}) = (nanstd(all_AOC_high.(types{iType}),[],3)./sqrt(size(all_AOC_high.(types{iType}),3)))';
%     SEM_AOC_con_low.(types{iType}) =  (nanstd(all_AOC_con_low.(types{iType}),[],3)./sqrt(size(all_AOC_con_low.(types{iType}),3)))';
%     SEM_AOC_con_high.(types{iType}) =  (nanstd(all_AOC_con_high.(types{iType}),[],3)./sqrt(size(all_AOC_con_high.(types{iType}),3)))';
%
%     norm_SEM_AOC_low.(types{iType}) = (nanstd(norm_all_AOC_low.(types{iType}),[],3)./sqrt(size(norm_all_AOC_low.(types{iType}),3)))';
%     norm_SEM_AOC_high.(types{iType}) = (nanstd(norm_all_AOC_high.(types{iType}),[],3)./sqrt(size(norm_all_AOC_high.(types{iType}),3)))';
%     norm_SEM_AOC_con_low.(types{iType}) =  (nanstd(norm_all_AOC_con_low.(types{iType}),[],3)./sqrt(size(norm_all_AOC_con_low.(types{iType}),3)))';
%     norm_SEM_AOC_con_high.(types{iType}) =  (nanstd(norm_all_AOC_con_high.(types{iType}),[],3)./sqrt(size(norm_all_AOC_con_high.(types{iType}),3)))';
%
%     % individual subjects.
%     for iSub = 1:length(subjects)
%         Rats.(subjects{iSub}).all_AOC_low_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).all_AOC_low.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).all_AOC_low.(types{iType}),3)))';
%         Rats.(subjects{iSub}).all_AOC_high_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).all_AOC_high.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).all_AOC_high.(types{iType}),3)))';
%         Rats.(subjects{iSub}).all_AOC_con_low_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).all_AOC_con_low.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).all_AOC_con_low.(types{iType}),3)))';
%         Rats.(subjects{iSub}).all_AOC_con_high_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).all_AOC_con_high.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).all_AOC_con_high.(types{iType}),3)))';
%
%         Rats.(subjects{iSub}).norm_all_AOC_low_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).norm_all_AOC_low.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).norm_all_AOC_low.(types{iType}),3)))';
%         Rats.(subjects{iSub}).norm_all_AOC_high_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).norm_all_AOC_high.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).norm_all_AOC_high.(types{iType}),3)))';
%         Rats.(subjects{iSub}).norm_all_AOC_con_low_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).norm_all_AOC_con_low.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).norm_all_AOC_con_low.(types{iType}),3)))';
%         Rats.(subjects{iSub}).norm_all_AOC_con_high_SEM.(types{iType}) = (nanstd(Rats.(subjects{iSub}).norm_all_AOC_con_high.(types{iType}),[],3)./sqrt(size(Rats.(subjects{iSub}).norm_all_AOC_con_high.(types{iType}),3)))';
%     end
%     %% get actual stats
%     ks = []; ksh = [];
%     for iSite = 1:length(sites)
%         labels = {'ipsi', 'contra', 'control'};
%         h = kstest(norm_all_AOC_low.(types{iType})(:, iSite,:));
%         if h
%             disp('***************************************************************')
%             disp(['KS test FAIL for low ' sites{iSite}])
%             disp('***************************************************************')
%             ks = [ks ; 1];
%         end
%         ks = [ks; 0];
%
%
%         hh = kstest(norm_all_AOC_high.(types{iType})(:, iSite,:));
%         if hh
%             disp('***************************************************************')
%             disp(['KS test FAIL for high ' sites{iSite}])
%             disp('***************************************************************')
%             ksh = [ksh ; 1];
%         end
%         ksh = [ksh; 0];
%     end
%
%     %% tests for differnences
%     for iSite = 1:length(sites)
%         this_ipsi = squeeze(norm_all_AOC_low.(types{iType})(2, iSite,:));
%         this_ipsi(isnan(this_ipsi)) = [];
%         this_con= squeeze(norm_all_AOC_low.(types{iType})(3, iSite,:));
%         this_con(isnan(this_con)) = [];
%         this_ctr = squeeze(norm_all_AOC_low.(types{iType})(5, iSite,:));
%         this_ctr(isnan(this_ctr)) = [];
%
%         h_this_ipsi = squeeze(norm_all_AOC_high.(types{iType})(2, iSite,:));
%         h_this_ipsi(isnan(h_this_ipsi)) = [];
%         h_this_con= squeeze(norm_all_AOC_high.(types{iType})(3, iSite,:));
%         h_this_con(isnan(h_this_con)) = [];
%         h_this_ctr = squeeze(norm_all_AOC_high.(types{iType})(5, iSite,:));
%         h_this_ctr(isnan(h_this_ctr)) = [];
%
%
%         if sum(ks)>=1 || sum(ksh)>=1
%             [p_ip_con.(types{iType})(iSite), h_ip_con.(types{iType})(iSite)] = signrank(this_ipsi, this_con);
%             [p_ip_ctr.(types{iType})(iSite), h_ip_ctr.(types{iType})(iSite)] = signrank(this_ipsi, this_ctr);
%             [p_con_ctr.(types{iType})(iSite), h_con_ctr.(types{iType})(iSite)] = signrank(this_con, this_ctr);
%
%             [h_p_ip_con.(types{iType})(iSite), h_h_ip_con.(types{iType})(iSite)] = signrank(h_this_ipsi, h_this_con);
%             [h_p_ip_ctr.(types{iType})(iSite), h_h_ip_ctr.(types{iType})(iSite)] = signrank(h_this_ipsi, h_this_ctr);
%             [h_p_con_ctr.(types{iType})(iSite), h_h_con_ctr.(types{iType})(iSite)] = signrank(h_this_con, h_this_ctr);
%         else
%             disp('Using T-Test')
%             [h_ip_con(iSite).(types{iType}), p_ip_con.(types{iType})(iSite), ~, l_stats_ip_con.(types{iType})(iSite)] = ttest2(this_ipsi, this_con);
%             [h_ip_ctr(iSite).(types{iType}), p_ip_ctr.(types{iType})(iSite), ~,l_stats_ip_ctr.(types{iType})(iSite)] = ttest2(this_ipsi, this_ctr);
%             [h_con_ctr(iSite).(types{iType}), p_con_ctr.(types{iType})(iSite), ~,l_stats_con_ctr.(types{iType})(iSite)] = ttest2(this_con, this_ctr);
%
%             [h_h_ip_con(iSite).(types{iType}), h_p_ip_con(iSite).(types{iType}), ~, h_stats_ip_con.(types{iType})(iSite)] = ttest2(h_this_ipsi, h_this_con);
%             [h_h_ip_ctr(iSite).(types{iType}), h_p_ip_ctr(iSite).(types{iType}), ~,h_stats_ip_ctr.(types{iType})(iSite)] = ttest2(h_this_ipsi, h_this_ctr);
%             [h_h_con_ctr(iSite).(types{iType}), h_p_con_ctr(iSite).(types{iType}), ~,h_stats_con_ctr.(types{iType})(iSite)] = ttest2(h_this_con, h_this_ctr);
%         end
%     end
%     %%
%     % sites = sites'
%     if sum(ks) >=1
%         fprintf('\nWilcoxin Sign Rank test\n')
%         fprintf('low Gamma\n')
%         fprintf('                                      ')
%         for iSite = 1:length(sites)
%             fprintf([sites{iSite} '        '] )
%             fprintf(repmat('\b', 1, length(sites{iSite})-2))
%         end
%         fprintf(['\nIpsilateral   vs. Contralateral:    P:' num2str(p_ip_con.(types{iType}), '%10.4f') '\n' ])
%         fprintf(['Ipsilateral   vs. Control:          P:' num2str(p_ip_ctr.(types{iType}), '%10.4f') '\n' ])
%         fprintf(['Contralateral vs. Control:          P:' num2str(p_con_ctr.(types{iType}), '%10.4f') '\n' ])
%
%         fprintf('\nHigh Gamma\n')
%         fprintf('                                      ')
%         for iSite = 1:length(sites)
%             fprintf([sites{iSite} '        '] )
%             fprintf(repmat('\b', 1, length(sites{iSite})-2))
%         end
%         fprintf(['\nIpsilateral   vs. Contralateral:    P:' num2str(h_p_ip_con.(types{iType}), '%10.4f') '\n' ])
%         fprintf(['Ipsilateral   vs. Control:          P:' num2str(h_p_ip_ctr.(types{iType}), '%10.4f') '\n' ])
%         fprintf(['Contralateral vs. Control:          P:' num2str(h_p_con_ctr.(types{iType}), '%10.4f') '\n' ])
%     else
%         fprintf('\nPaired T-Test\n')
%         for iSite = 1:length(sites)
%             fprintf(['\nLow Gamma   ' sites{iSite} '\n'])
%             fprintf([sites{iSite} ' Ipsilateral   vs. Contralateral:   df(' num2str(l_stats_ip_con.(types{iType})(iSite).df) ')   t:' num2str(l_stats_ip_con.(types{iType})(iSite).tstat, '%4.4f') '  P:' num2str(p_ip_con.(types{iType})(iSite), '%4.4f') '\n' ])
%             fprintf([sites{iSite} ' Ipsilateral   vs. Control:         df(' num2str(l_stats_ip_ctr.(types{iType})(iSite).df) ')   t:' num2str(l_stats_ip_ctr.(types{iType})(iSite).tstat, '%4.4f') '  P:' num2str(p_ip_ctr.(types{iType})(iSite), '%4.4f') '\n' ])
%             fprintf([sites{iSite} ' Contralateral vs. Control:         df(' num2str(l_stats_con_ctr.(types{iType})(iSite).df) ')   t:' num2str(l_stats_con_ctr.(types{iType})(iSite).tstat, '%4.4f') '      P:' num2str(p_con_ctr.(types{iType})(iSite), '%4.4f') '\n' ])
%
%             fprintf('\nPaired T-Test\n')
%             fprintf(['High Gamma    ' sites{iSite} '\n'])
%             fprintf(['Ipsilateral   vs. Contralateral:   df(' num2str(h_stats_ip_con(iSite).df) ')   t:' num2str(h_stats_ip_con(iSite).tstat, '%4.4f') '  P:' num2str(h_p_ip_con(iSite), '%4.4f') '\n' ])
%             fprintf(['Ipsilateral   vs. Control:         df(' num2str(h_stats_ip_ctr(iSite).df) ')   t:' num2str(h_stats_ip_ctr(iSite).tstat, '%4.4f') '  P:' num2str(h_p_ip_ctr(iSite), '%4.4f') '\n' ])
%             fprintf(['Contralateral vs. Control:         df(' num2str(h_stats_con_ctr(iSite).df) ')   t:' num2str(h_stats_con_ctr(iSite).tstat, '%4.4f') '   P:' num2str(h_p_con_ctr(iSite), '%4.4f') '\n' ])
%         end
%     end
%
%
%     %% make a bar plot
%
%     for iFig = 1:2
%         if iFig ==1
%             F_id = 'Four';
%             s_idx = [1,2,3,5,7]; % corresponds to the PL, OFC, NAc, and CG
%             bar_names = {'PL', 'IL', 'OFC', 'NAc', 'CG'};
%
%         elseif iFig == 2
%             F_id = 'Piri';
%             s_idx = [3:6]; % corresponds to OFC, OFC_Piri, NAc, and NAc_Piri
%             bar_names = {'OFC', 'Piri OFC', 'NAc', 'Piri NAc'};
%         end
%         %             to_plot = [1,2,3,5,7];
%         bar_c_ord = linspecer(5);
%         %     error_bar_low = (nanstd(num_gamma.all.low,1,3)./sqrt(length(num_gamma.all.low)))';
%         %     error_bar_high = (nanstd(num_gamma.all.high,1,3)./sqrt(length(num_gamma.all.high)))';
%
%         if strcmp(cfg.plot_type, 'raw')
%             % shift control to the first column
%             bar_temp_Low = circshift(mean_AOC_low.(types{iType}),1,2);
%             bar_temp_High = circshift(mean_AOC_high.(types{iType}),1,2);
%             SEM_bar_temp_Low = circshift(SEM_AOC_low.(types{iType}),1,2);
%             SEM_bar_temp_High = circshift(SEM_AOC_high.(types{iType}),1,2);
%         elseif strcmp(cfg.plot_type, 'norm')
%             bar_temp_Low = circshift(norm_mean_AOC_low.(types{iType}),1,2);
%             bar_temp_High = circshift(norm_mean_AOC_high.(types{iType}),1,2);
%             SEM_bar_temp_Low = circshift(norm_SEM_AOC_low.(types{iType}),1,2);
%             SEM_bar_temp_High = circshift(norm_SEM_AOC_high.(types{iType}),1,2);
%         end
%
%         h_low = errorbar_groups(bar_temp_Low(s_idx,[1,3,4])',SEM_bar_temp_Low(s_idx,[1,3,4])', 'bar_colors', bar_c_ord, 'bar_names', bar_names, 'FigID', 100);
%         title(['All-low-' iComp{1} '-' types{iType} '-' cfg.plot_type])
%
%         %same for high gamma
%         h_high = errorbar_groups(bar_temp_High(s_idx,[1,3,4])',SEM_bar_temp_High(s_idx,[1,3,4])', 'bar_colors', bar_c_ord, 'bar_names', bar_names, 'FigID', 200);
%         title(['All-high-' iComp{1} '-' types{iType} '-' cfg.plot_type])
%
%         close all
%
%         %% save figures
%         mkdir(PARAMS.inter_dir, 'AOC_fit')
%         if isunix
%             saveas(h_low, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_low_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type])
%             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_low_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type], 'png')
%             saveas_eps(['AOC_Summary_low_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type],[PARAMS.inter_dir '/AOC_fit/'])
%         else
%             saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_low_' F_id '_' cfg.pot_trk '_' types{iType} '_' iComp{1} '_' cfg.plot_type])
%             saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_low_' F_id '_' cfg.pot_trk '_' types{iType} '_' iComp{1} '_' cfg.plot_type], 'png')
%             saveas_eps(['AOC_Summary_low_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type],[PARAMS.inter_dir '\AOC_fit\'])
%         end
%
%         % same for high gamma
%         if isunix
%             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_high_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type])
%             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_high_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type], 'png')
%             saveas_eps(['AOC_Summary_high_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type],[PARAMS.inter_dir '/AOC_fit/'])
%         else
%             saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_high_' F_id '_' cfg.pot_trk '_' types{iType} '_' iComp{1} '_' cfg.plot_type])
%             saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_high_' F_id '_' cfg.pot_trk '_' types{iType} '_' iComp{1} '_' cfg.plot_type], 'png')
%             saveas_eps(['AOC_Summary_high_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type],[PARAMS.inter_dir '\AOC_fit\'])
%         end
%
%         close all
%
%         %% to do.  make a plot for each subject.
%         %         for iSub = 1:length(subjects)
%         %             %            figure(iSub)
%         %             bar_temp_Low = circshift(Rats.(subjects{iSub}).all_AOC_low_mean.(types{iType}),1,2);
%         %             bar_temp_High = circshift(Rats.(subjects{iSub}).all_AOC_high_mean.(types{iType}),1,2);
%         %             SEM_bar_temp_Low = circshift(Rats.(subjects{iSub}).all_AOC_low_SEM.(types{iType}),1,2);
%         %             SEM_bar_temp_High = circshift(Rats.(subjects{iSub}).all_AOC_high_SEM.(types{iType}),1,2);
%         %
%         %             h_low = errorbar_groups(bar_temp_Low(to_plot,[1,3,4])',SEM_bar_temp_Low(to_plot,[1,3,4])', 'bar_colors', bar_c_ord, 'bar_names', bar_names, 'FigID', [iSub*100]);
%         %             title([subjects{iSub} '-low-' iComp{1} '-' types{iType}])
%         %             h_high = errorbar_groups(bar_temp_High(to_plot,[1,3,4])',SEM_bar_temp_High(to_plot,[1,3,4])', 'bar_colors', bar_c_ord, 'bar_names', bar_names, 'FigID',[iSub*1000]);
%         %             title([subjects{iSub} '-high-' iComp{1} '-' types{iType}])
%         %         end
%
%
%
%
%
%
%     end
%%
%         %%
%         ks = []; ksh = [];
%         for iPhase = 1:length(PARAMS.Phases);
%             for iSite = 1:length(sites)
%                 labels = { 'pre', 'ipsi', 'contra', 'post','control'};
%                 this_val = norm_all_AOC_low.(types{iType})(iPhase, iSite,:);
%                 this_val(isnan(this_val)) = [];
%                 [hl,p] = kstest(squeeze(this_val)');
%                 if hl
%                     disp('***************************************************************')
%                     disp(['KS test FAIL for low ' labels{iPhase}])
%                     disp('***************************************************************')
%                     ks = [ks ; 1];
%                 end
%                 ks = [ks; 0];
%
%                 this_val = norm_all_AOC_high.(types{iType})(iPhase, iSite,:);
%                 this_val(isnan(this_val)) = [];
%                 [hh,~] = kstest(squeeze(this_val));
%                 if hh
%                     disp('***************************************************************')
%                     disp(['KS test FAIL for high ' sites{iPhase}])
%                     disp('***************************************************************')
%                     ksh = [ksh ; 1];
%                 end
%                 ksh = [ksh; 0];
%             end
%         end

%         %% probably no t passing KS
%         for iSite  = 1:length(sites)
%             this_ipsi = squeeze(norm_all_AOC_low.(types{iType})(2, iSite,:));
%             this_ipsi(isnan(this_ipsi)) = [];
%             this_con= squeeze(norm_all_AOC_low.(types{iType})(3, iSite,:));
%             this_con(isnan(this_con)) = [];
%             this_ctr = squeeze(norm_all_AOC_low.(types{iType})(5, iSite,:));
%             this_ctr(isnan(this_ctr)) = [];
%
%             h_this_ipsi = squeeze(norm_all_AOC_high.(types{iType})(2, iSite,:));
%             h_this_ipsi(isnan(h_this_ipsi)) = [];
%             h_this_con= squeeze(norm_all_AOC_high.(types{iType})(3, iSite,:));
%             h_this_con(isnan(h_this_con)) = [];
%             h_this_ctr = squeeze(norm_all_AOC_high.(types{iType})(5, iSite,:));
%             h_this_ctr(isnan(h_this_ctr)) = [];
%
%
%             if sum(ks)>=1
%
%                 [ l_p_ip_con.(types{iType}).(sites{iSite}), l_h_ip_con.(types{iType}).(sites{iSite})] = signrank(this_ipsi, this_con);
%                 [ l_p_ip_ctr.(types{iType}).(sites{iSite}), l_h_ip_ctr.(types{iType}).(sites{iSite})] = signrank(this_ipsi, this_ctr);
%                 [ l_p_con_ctr.(types{iType}).(sites{iSite}), l_h_con_ctr.(types{iType}).(sites{iSite})] = signrank(this_con, this_ctr);
%
%                 [ h_p_ip_con.(types{iType}).(sites{iSite}), h_h_ip_con.(types{iType}).(sites{iSite})] = signrank(h_this_ipsi, h_this_con);
%                 [ h_p_ip_ctr.(types{iType}).(sites{iSite}), h_h_ip_ctr.(types{iType}).(sites{iSite})] = signrank(h_this_ipsi, h_this_ctr);
%                 [ h_p_con_ctr.(types{iType}).(sites{iSite}),h_h_con_ctr.(types{iType}).(sites{iSite})] = signrank(h_this_con, h_this_ctr);
%             else
%                 disp('Using T-Test')
%                 [h_ip_con, p_ip_con, ~ ,l_stats_ip_con] = ttest2(this_ipsi, this_con);
%                 [h_ip_ctr, p_ip_ctr, ~ ,l_stats_ip_ctr] = ttest2(this_ipsi, this_ctr);
%                 [h_con_ctr, p_con_ctr, ~ , l_stats_con_ctr] = ttest2(all_count_low(:, 4), all_count_low(:, 1));
%
%                 [h_h_ip_con, h_p_ip_con,~ , h_stats_ip_con] = ttest2(all_count_high(:, 3), all_count_high(:, 4));
%                 [h_h_ip_ctr, h_p_ip_ctr, ~, h_stats_ip_ctr] = ttest2(all_count_high(:, 3), all_count_high(:, 1));
%                 [h_h_con_ctr, h_p_con_ctr, ~, h_stats_con_ctr] = ttest2(all_count_high(:, 4), all_count_high(:, 1));
%             end
%         end

%% switch between normalized to control and raw
%     switch cfg.plot_type
%
%         case 'raw'
%             % plot
%             figure(iType)
%             subplot(4,1,1)
%             b= bar(AOC_low.(types{iType}));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             title(types{iType})
%             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
%             leg_val = PARAMS.Phases; leg_val{5} = 'Control';
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims])
%
%             subplot(4,1,2)
%             b = bar(AOC_high.(types{iType}));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims])
%
%             subplot(4,1,3)
%             b = bar(AOC_con_low.(types{iType}));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):25:cfg.ylims(2)])
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims])
%
%             subplot(4,1,4)
%             b = bar(AOC_con_high.(types{iType}));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             set(gca, 'xticklabel', sites, 'ytick', [cfg.ylims(1):25:cfg.ylims(2)])
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims])
%
%         case 'norm'
%             % plot
%             figure(iType)
%             subplot(4,1,1)
%             b= bar(norm_low.(types{iType})(:,1:4));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             title(types{iType})
%             set(gca, 'xticklabel', sites)
%             leg_val = PARAMS.Phases; leg_val{5} = 'Control';
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims_norm])
%
%             subplot(4,1,2)
%             b = bar(norm_high.(types{iType})(:,1:4));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             set(gca, 'xticklabel', sites)
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims_norm])
%
%             subplot(4,1,3)
%             b = bar(norm_con_low.(types{iType})(:,1:4));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             set(gca, 'xticklabel', sites)
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims_norm])
%
%             subplot(4,1,4)
%             b = bar(norm_con_high.(types{iType})(:,1:4));
%             for iPhase = 1:5
%                 set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%             end
%             set(gca, 'xticklabel', sites)
%             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%             ylim([cfg.ylims_norm])
%
%             %% save the figure
%             if isunix
%                 saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type])
%                 saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'png')
%             else
%                 saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type])
%                 saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type], 'png')
%             end
%         close all
%         %%%%%%%%%%%%%%%%%%% Multisite 4-site Vs cross Piriform figure %%%%%%%%%%%%%
%         sites{4} = 'PC_O_F_C'; sites{6} = 'PC_N_A_c';
%         cfg_plt1.pos = [600 50 560*1.4 560*1.8];
%         cfg_plt1.ft_size = 18;
%
%         for iFig = 1:2
%             if iFig ==1
%                 F_id = 'Four';
%                 s_idx = [1,2,3,5,7]; % corresponds to the PL, OFC, NAc, and CG
%             elseif iFig == 2
%                 F_id = 'Piri';
%                 s_idx = [3:6]; % corresponds to OFC, OFC_Piri, NAc, and NAc_Piri
%             end
%             %% Gerenate the summary figure for the PL, OFC, NAc, CG
%             switch cfg.plot_type
%
%                 case 'raw'
%                     % plot
%                     figure(iType)
%                     subtightplot(9,3,[2,3,5,6],0.1) % strange subplots make for nice figures.  Ignore white space...
%                     b= bar(AOC_low.(types{iType})(s_idx,:));
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['low gamma AOC (' num2str(cfg.power_ratio.gamma_freq(1,1)) '-' num2str(cfg.power_ratio.gamma_freq(1,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
%                     ylim([cfg.ylims])
%                     SetFigure(cfg_plt1, gcf)
%
%                     subtightplot(9,3,[8,9,11,12],0.1)
%                     b = bar(AOC_high.(types{iType})(s_idx,:));
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['high gamma AOC (' num2str(cfg.power_ratio.gamma_freq(2,1)) '-' num2str(cfg.power_ratio.gamma_freq(2,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
%                     ylim([cfg.ylims])
%                     SetFigure(cfg_plt1, gcf)
%
%                     subtightplot(9,3,[17,18,20,21],0.1)
%                     b = bar(AOC_con_low.(types{iType})(s_idx,:));
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['low control AOC (' num2str(cfg.power_ratio.contrast(1,1)) '-' num2str(cfg.power_ratio.contrast(1,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
%                     %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%                     ylim([cfg.ylims])
%                     SetFigure(cfg_plt1, gcf)
%
%                     subtightplot(9,3,[23,24,26,27],0.1)
%                     b = bar(AOC_con_high.(types{iType})(s_idx,:));
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['high control AOC (' num2str(cfg.power_ratio.contrast(2,1)) '-' num2str(cfg.power_ratio.contrast(2,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx), 'ytick', [cfg.ylims(1):50:cfg.ylims(2)])
%                     %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%                     ylim([cfg.ylims])
%
%                     SetFigure(cfg_plt1, gcf)
%                     %%
%                 case 'norm'
%                     % plot
%                     figure(iType)
%                     subtightplot(12,1,1:3,0.1)
%                     b= bar(norm_low.(types{iType})(s_idx,:), 'BaseValue', 1);
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['low gamma AOC (' num2str(cfg.power_ratio.gamma_freq(1,1)) '-' num2str(cfg.power_ratio.gamma_freq(1,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx))
%                     %                 ylim([cfg.ylims_norm])
%
%                     subtightplot(12,1,4:6,0.1)
%                     b = bar(norm_high.(types{iType})(s_idx,:), 'BaseValue', 1);
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['high gamma AOC (' num2str(cfg.power_ratio.gamma_freq(2,1)) '-' num2str(cfg.power_ratio.gamma_freq(2,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx))
%                     %                 ylim([cfg.ylims_norm])
%
%
%                     subtightplot(12,1,7:9,0.1)
%                     b = bar(norm_con_low.(types{iType})(s_idx,:), 'BaseValue', 1);
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['low control AOC (' num2str(cfg.power_ratio.contrast(1,1)) '-' num2str(cfg.power_ratio.contrast(1,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx))
%                     %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%                     %                 ylim([cfg.ylims_norm])
%
%                     subtightplot(12,1,10:12,0.1)
%                     b = bar(norm_con_high.(types{iType})(s_idx,:), 'BaseValue', 1);
%                     for iPhase = 1:5
%                         set(b(iPhase), 'FaceColor', c_ord(iPhase,:))
%                     end
%                     title(['high control AOC (' num2str(cfg.power_ratio.contrast(2,1)) '-' num2str(cfg.power_ratio.contrast(2,2)) 'Hz)'], 'fontweight', 'normal');
%                     set(gca, 'xticklabel', sites(s_idx))
%                     %             legend(leg_val, 'location', 'eastoutside', 'orientation', 'vertical');
%                     %                 ylim([cfg.ylims_norm])
%
%                     %                 cfg_plt1.pos = [600 50 560*1.4 560*1.8];
%                     cfg_plt1.ft_size = 18;
%                     SetFigure(cfg_plt1, gcf)
%             end
%             %%
%             mkdir(PARAMS.inter_dir, 'AOC_fit')
%             if isunix
%                 saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type])
%                 saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type], 'png')
%                 saveas_eps(['AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type],[PARAMS.inter_dir '/AOC_fit/'])
%                 %             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'epsc')
%
%             else
%                 saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType} '_' iComp{1} '_' cfg.plot_type])
%                 saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType} '_' iComp{1} '_' cfg.plot_type], 'png')
%                 %             saveas(gcf, [PARAMS.inter_dir '\AOC_fit\AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType} '_' cfg.plot_type], 'epsc')
%                 saveas_eps(['AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' iComp{1} '_' cfg.plot_type],[PARAMS.inter_dir '\AOC_fit\'])
%
%             end
%             close all
%
%         end
%     end
% end
% % create the legend values and add some space


