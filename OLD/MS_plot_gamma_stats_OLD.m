function stats = MS_gamma_stats(cfg_in, Events)
%%
%
%
%
%
%



%% set up defaults

cfg_def = [];
cfg_def.pot_trk = 'pot';

cfg = ProcessConfig2(cfg_def, cfg_in)

global PARAMS

%% collect all the gamma event counts
out.low = [];  out.high = [];
out_norm.low = []; out_norm.high= [];
rate.low =[]; rate.high = [];

sites = {'PL'    'IL'    'OFC'    'Piri_O'    'NAc'    'Piri_N'    'CG'};
bands = {'low', 'high'};
phases = [PARAMS.Phases, 'control',];
sub_list = fieldnames(Events);

for iBand = 1:2
    for iSite = 1:length(sites)
        for iPhase = 1:length(phases)
            all_events_len.(bands{iBand}).(sites{iSite}).(phases{iPhase}) =[];
        end
    end
end

for iBand = 1:length(bands)
    for iSub = 1:length(sub_list)
        sess_list = fieldnames(Events.(sub_list{iSub}));
        % pre allocate temp fields
        temp = NaN(length(PARAMS.Phases)+1, length(sites),4);
        temp_norm = NaN(length(PARAMS.Phases)+1, length(sites),4);
        single_subject.rate = NaN(length(PARAMS.Phases)+1, length(sites),4);

        for iSess = 1:length(sess_list);
            these_sites = fieldnames(Events.(sub_list{iSub}).(sess_list{iSess}));
            for iSite = 1:length(these_sites)
                       site_idx = strcmp(sites, these_sites{iSite}(1:end-4)); 
%                        if sum(site_idx) ==1
                           site_idx = find(site_idx ==1);
                           site_field = strcmp(these_sites, [sites{site_idx} '_' cfg.pot_trk]); % find the correct subfield;
                         if  sum(site_field) ==1
                           site_field = find(site_field==1);
                           
                           for iPhase = 1:4
                               temp(iPhase,site_idx,iSess) = length(Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).tstart);
                               single_subject.rate(iPhase,site_idx , iSess) = Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).rate;
                                
                               evt_lens = Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).tend - Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{iPhase}).(bands{iBand}).tstart;
                               all_events_len.(bands{iBand}).(sites{site_idx}).(phases{iPhase}) = cat(1,all_events_len.(bands{iBand}).(sites{site_idx}).(phases{iPhase}),evt_lens);
                           end
                           temp(5, site_idx, iSess) = nanmean([temp(1,site_idx,iSess),temp(4,site_idx,iSess)]); 
                           single_subject.rate(5,site_idx , iSess) = Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{5}).(bands{iBand}).rate;
                           evt_lens = Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{5}).(bands{iBand}).tend - Events.(sub_list{iSub}).(sess_list{iSess}).(these_sites{site_field}).(phases{5}).(bands{iBand}).tstart;
                           all_events_len.(bands{iBand}).(sites{site_idx}).(phases{5}) = cat(1,all_events_len.(bands{iBand}).(sites{site_idx}).(phases{5}),evt_lens);
                        
                           temp_norm(:,site_idx,iSess) = temp(:,site_idx,iSess)./temp(5,site_idx,iSess);
                          
                           single_subject.num = temp;
                           single_subject.num_norm = temp_norm;

                        
                         end
%                     if iSess == 1
%                         out.(bands{iBand})(iSite, :,iSess) = temp.(bands{iBand})(iSite, :,iSess);
%                         norm.(bands{iBand})(iSite, :,iSess) = temp_norm.(bands{iBand})(iSite, :,iSess);
%                     else
%                         out.(bands{iBand}) = [out.(bands{iBand}); temp.(bands{iBand}).(site_list{iSite})(iSess, :)];
%                         norm.(bands{iBand}) = [norm.(bands{iBand}); temp_norm.(bands{iBand}).(site_list{iSite})(iSess, :)];
%                     end
            end
%             temp = []; temp_norm;
        end
        all_subs.(sub_list{iSub}).(bands{iBand}) = single_subject;

        out.(bands{iBand}) = cat(3,out.(bands{iBand}), temp);
        out_norm.(bands{iBand}) = cat(3,out_norm.(bands{iBand}), temp_norm);
        rate.(bands{iBand}) = cat(3,rate.(bands{iBand}), single_subject.rate);
        temp = []; temp_norm =[];

    end
end

%% make sure there are no infs in the normalized version

out_norm.low(out_norm.low == inf) = NaN;
out_norm.high(out_norm.high == inf) = NaN;


%% collect the output 
% out_norm.low = circshift(

stats_file = fopen([PARAMS.stats_dir 'Count_stats.txt'], 'w');
types= {'Four', 'Piri'}; 
for iTypes = 1:length(types)
    if strcmp(types{iTypes}, 'Four')
        s_idx = [1 2 3 5 7];
    else strcmp(types{iTypes}, 'Piri')
        s_idx = [3 4 5 6];
    end
cfg_stats = [];
        cfg_stats.title = ['Count low gamma ' types{iTypes}] ;
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out.low)
close all

cfg_stats = [];
        cfg_stats.title = ['Count high gamma ' types{iTypes}];
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out.high)

%normalized
cfg_stats = [];
        cfg_stats.title = ['Count low gamma norm ' types{iTypes}];
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out_norm.low)
close all

cfg_stats = [];
        cfg_stats.title = ['Count high gamma norm ' types{iTypes}];
        cfg_stats.method= 'median';
        cfg_stats.row_names= {'PL'  'IL'  'OFC'  'Piri OFC'  'NAc'  'Piri NAc'  'CG'};
        cfg_stats.col_names= {'pre'  'ipsi'  'contra'  'post', 'control'};
        cfg_stats.s_idx= s_idx;
        cfg_stats.ft_size= 20;
        cfg_stats.save_dir= [PARAMS.inter_dir 'Count'];
        cfg_stats.stats_dir = stats_file;
MS_stats(cfg_stats,out_norm.high)
close all
end

%% Get the length and rate 
all_len.low = []; all_rate.low = [];
all_len.high = []; all_rate.high = [];
fprintf('\nNumer of events\n')
for iBand= 1:length(bands)
    fprintf(['\n-------- ' bands{iBand} '---------\n'])
    for iSite = 1:length(sites)
        fprintf([sites{iSite} ':       '])
        fprintf(repmat('\b', 1, length(sites{iSite})))
        for iPhase = 1:length(phases)
            all_len.(bands{iBand})(iSite, iPhase) = median(all_events_len.(bands{iBand}).(sites{iSite}).(phases{iPhase}));
            all_len_std.(bands{iBand})(iSite, iPhase) = nanstd(all_events_len.(bands{iBand}).(sites{iSite}).(phases{iPhase}))./sqrt(size(all_events_len.(bands{iBand}).(sites{iSite}).(phases{iPhase}),1));
            fprintf([phases{iPhase} ' mean= %.2f +/- %.2f   '], all_len.(bands{iBand})(iSite, iPhase)*1000, all_len_std.(bands{iBand})(iSite, iPhase)*1000)
        end
        fprintf('\n')
    end
end

% make a box plot
this_data = [all_events_len.low.PL.pre; all_events_len.low.PL.contra, ]
this_group = 
box([all_events_len.low.PL.pre; all_events_len.low.OFC.pre])

data = rand(20,24);
month = repmat({'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'},1,2);
simobs = [repmat({'sim'},1,12),repmat({'obs'},1,12)];
boxplot(data,{month,simobs},'colors',repmat('rb',1,12),'factorgap',[5 2],'labelverbosity','minor');

%% all sites
x = [all_events_len.(bands{iBand}).(sites{1}).(phases{1})', all_events_len.(bands{iBand}).(sites{1}).(phases{2})', all_events_len.(bands{iBand}).(sites{1}).(phases{3})', all_events_len.(bands{iBand}).(sites{1}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{2}).(phases{1})', all_events_len.(bands{iBand}).(sites{2}).(phases{2})', all_events_len.(bands{iBand}).(sites{2}).(phases{3})', all_events_len.(bands{iBand}).(sites{2}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{3}).(phases{1})', all_events_len.(bands{iBand}).(sites{3}).(phases{2})', all_events_len.(bands{iBand}).(sites{3}).(phases{3})', all_events_len.(bands{iBand}).(sites{3}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{4}).(phases{1})', all_events_len.(bands{iBand}).(sites{4}).(phases{2})', all_events_len.(bands{iBand}).(sites{4}).(phases{3})', all_events_len.(bands{iBand}).(sites{4}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{5}).(phases{1})', all_events_len.(bands{iBand}).(sites{5}).(phases{2})', all_events_len.(bands{iBand}).(sites{5}).(phases{3})', all_events_len.(bands{iBand}).(sites{5}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{6}).(phases{1})', all_events_len.(bands{iBand}).(sites{6}).(phases{2})', all_events_len.(bands{iBand}).(sites{6}).(phases{3})', all_events_len.(bands{iBand}).(sites{6}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{7}).(phases{1})', all_events_len.(bands{iBand}).(sites{7}).(phases{2})', all_events_len.(bands{iBand}).(sites{7}).(phases{3})', all_events_len.(bands{iBand}).(sites{7}).(phases{4})',...
    ];
x = x.*1000;
group = [repmat(1,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{1}))), repmat(2,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{2}))), repmat(3,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{3}))), repmat(4,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{4}))),...
    repmat(5,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{1}))), repmat(6,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{2}))), repmat(7,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{3}))), repmat(8,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{4}))),...
    repmat(9,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{1}))), repmat(10,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{2}))), repmat(11,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{3}))), repmat(12,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{4}))),...
    repmat(13,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{1}))), repmat(14,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{2}))), repmat(15,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{3}))), repmat(16,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{4}))),...
    repmat(17,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{1}))), repmat(18,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{2}))), repmat(19,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{3}))), repmat(20,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{4}))),...
    repmat(21,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{1}))), repmat(22,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{2}))), repmat(23,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{3}))), repmat(24,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{4}))),...
    repmat(25,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{1}))), repmat(26,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{2}))), repmat(27,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{3}))), repmat(28,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{4}))),... 
];
p_1 = 1:7;
p_2 = p_1+.20;
p_3 = p_2+.20;
p_4 = p_3+.20;
positions = zeros(1,length([p_1, p_2, p_3, p_4]));
positions(1:4:end) = p_1;
positions(2:4:end) = p_2;
positions(3:4:end) = p_3;
positions(4:4:end) = p_4;

boxplot(x,group, 'positions', positions, 'symbol', '');
ylim([0 250])
set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16)) mean(positions(17:20)) mean(positions(21:24)) mean(positions(25:28))])
set(gca,'xticklabel',cfg_stats.row_names)
c_ord = linspecer(4); 
color = repmat(flipud(c_ord), 7, 1);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5);
end

c = get(gca, 'Children');

legend(c(1:4), 'pre', 'ipsi', 'contra', 'post' );
ylabel('Event length (ms)')
SetFigure([], gcf)
cfg.save_dir = [PARAMS.inter_dir 'Count']; 
save_name = 'all_site_length'; 
if isunix
    %           fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '/' save_name '\n\n'])
    saveas(gcf, [ cfg.save_dir '/' save_name])
    saveas(gcf, [cfg.save_dir '/' save_name], 'png')
    saveas_eps(save_name,[cfg.save_dir '/'])
    %             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'epsc')
else
    %           fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '\' save_name '\n\n'])
    saveas(gcf, [cfg.save_dir '\' save_name])
    saveas(gcf, [cfg.save_dir '\' save_name], 'png')
    saveas_eps(save_name,[cfg.save_dir '\'])
end
close all
%% four sites only
x = [all_events_len.(bands{iBand}).(sites{1}).(phases{1})', all_events_len.(bands{iBand}).(sites{1}).(phases{2})', all_events_len.(bands{iBand}).(sites{1}).(phases{3})', all_events_len.(bands{iBand}).(sites{1}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{2}).(phases{1})', all_events_len.(bands{iBand}).(sites{2}).(phases{2})', all_events_len.(bands{iBand}).(sites{2}).(phases{3})', all_events_len.(bands{iBand}).(sites{2}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{3}).(phases{1})', all_events_len.(bands{iBand}).(sites{3}).(phases{2})', all_events_len.(bands{iBand}).(sites{3}).(phases{3})', all_events_len.(bands{iBand}).(sites{3}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{5}).(phases{1})', all_events_len.(bands{iBand}).(sites{5}).(phases{2})', all_events_len.(bands{iBand}).(sites{5}).(phases{3})', all_events_len.(bands{iBand}).(sites{5}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{7}).(phases{1})', all_events_len.(bands{iBand}).(sites{7}).(phases{2})', all_events_len.(bands{iBand}).(sites{7}).(phases{3})', all_events_len.(bands{iBand}).(sites{7}).(phases{4})',...
     ];
x = x.*1000;
group = [repmat(1,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{1}))), repmat(2,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{2}))), repmat(3,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{3}))), repmat(4,1,length(all_events_len.(bands{iBand}).(sites{1}).(phases{4}))),...
    repmat(5,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{1}))), repmat(6,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{2}))), repmat(7,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{3}))), repmat(8,1,length(all_events_len.(bands{iBand}).(sites{2}).(phases{4}))),...
    repmat(9,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{1}))), repmat(10,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{2}))), repmat(11,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{3}))), repmat(12,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{4}))),...
    repmat(13,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{1}))), repmat(14,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{2}))), repmat(15,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{3}))), repmat(16,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{4}))),...
    repmat(17,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{1}))), repmat(18,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{2}))), repmat(19,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{3}))), repmat(20,1,length(all_events_len.(bands{iBand}).(sites{7}).(phases{4}))),...
 ];
p_1 = 1:5;
p_2 = p_1+.20;
p_3 = p_2+.20;
p_4 = p_3+.20;
positions = zeros(1,length([p_1, p_2, p_3, p_4]));
positions(1:4:end) = p_1;
positions(2:4:end) = p_2;
positions(3:4:end) = p_3;
positions(4:4:end) = p_4;

boxplot(x,group, 'positions', positions, 'symbol', '');
ylim([0 250])
set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16)) mean(positions(17:20)) ])
set(gca,'xticklabel',{'PL'  'IL'  'OFC'   'NAc'   'CG'})
c_ord = linspecer(4); 
color = repmat(flipud(c_ord), 5, 1);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5);
end

c = get(gca, 'Children');

legend(c(1:4), 'pre', 'ipsi', 'contra', 'post' );
ylabel('Event length (ms)')
SetFigure([], gcf)

cfg.save_dir = [PARAMS.inter_dir 'Count']; 
save_name = 'four_site_length'; 
if isunix
    %           fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '/' save_name '\n\n'])
    saveas(gcf, [ cfg.save_dir '/' save_name])
    saveas(gcf, [cfg.save_dir '/' save_name], 'png')
    saveas_eps(save_name,[cfg.save_dir '/'])
    %             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'epsc')
else
    %           fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '\' save_name '\n\n'])
    saveas(gcf, [cfg.save_dir '\' save_name])
    saveas(gcf, [cfg.save_dir '\' save_name], 'png')
    saveas_eps(save_name,[cfg.save_dir '\'])
end
close all
%% Piri sites only
x = [all_events_len.(bands{iBand}).(sites{3}).(phases{1})', all_events_len.(bands{iBand}).(sites{3}).(phases{2})', all_events_len.(bands{iBand}).(sites{3}).(phases{3})', all_events_len.(bands{iBand}).(sites{3}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{4}).(phases{1})', all_events_len.(bands{iBand}).(sites{4}).(phases{2})', all_events_len.(bands{iBand}).(sites{4}).(phases{3})', all_events_len.(bands{iBand}).(sites{4}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{5}).(phases{1})', all_events_len.(bands{iBand}).(sites{5}).(phases{2})', all_events_len.(bands{iBand}).(sites{5}).(phases{3})', all_events_len.(bands{iBand}).(sites{5}).(phases{4})',...
    all_events_len.(bands{iBand}).(sites{6}).(phases{1})', all_events_len.(bands{iBand}).(sites{6}).(phases{2})', all_events_len.(bands{iBand}).(sites{6}).(phases{3})', all_events_len.(bands{iBand}).(sites{6}).(phases{4})',...
     ];
x = x.*1000;
group = [repmat(1,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{1}))), repmat(2,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{2}))), repmat(3,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{3}))), repmat(4,1,length(all_events_len.(bands{iBand}).(sites{3}).(phases{4}))),...
    repmat(5,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{1}))), repmat(6,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{2}))), repmat(7,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{3}))), repmat(8,1,length(all_events_len.(bands{iBand}).(sites{4}).(phases{4}))),...
    repmat(9,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{1}))), repmat(10,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{2}))), repmat(11,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{3}))), repmat(12,1,length(all_events_len.(bands{iBand}).(sites{5}).(phases{4}))),...
    repmat(13,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{1}))), repmat(14,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{2}))), repmat(15,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{3}))), repmat(16,1,length(all_events_len.(bands{iBand}).(sites{6}).(phases{4}))),...
 ];
p_1 = 1:4;
p_2 = p_1+.20;
p_3 = p_2+.20;
p_4 = p_3+.20;
positions = zeros(1,length([p_1, p_2, p_3, p_4]));
positions(1:4:end) = p_1;
positions(2:4:end) = p_2;
positions(3:4:end) = p_3;
positions(4:4:end) = p_4;

boxplot(x,group, 'positions', positions, 'symbol', '');
ylim([0 250])
set(gca,'xtick',[mean(positions(1:4)) mean(positions(5:8)) mean(positions(9:12)) mean(positions(13:16)) ])
set(gca,'xticklabel',{'OFC'  'Piri OFC'  'NAc'  'Piri NAc'})
c_ord = linspecer(4); 
color = repmat(flipud(c_ord), 4, 1);
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5);
end

c = get(gca, 'Children');

legend(c(1:4), 'pre', 'ipsi', 'contra', 'post' );
ylabel('Event length (ms)')
SetFigure([], gcf)

cfg.save_dir = [PARAMS.inter_dir 'Count']; 
save_name = 'piri_site_length'; 
if isunix
    %           fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '/' save_name '\n\n'])
    saveas(gcf, [ cfg.save_dir '/' save_name])
    saveas(gcf, [cfg.save_dir '/' save_name], 'png')
    saveas_eps(save_name,[cfg.save_dir '/'])
    %             saveas(gcf, [PARAMS.inter_dir '/AOC_fit/AOC_Summary_' F_id '_' cfg.pot_trk '_' types{iType}  '_' cfg.plot_type], 'epsc')
else
    %           fprintf(cfg.stats_dir,['\n\nSaving output in:      '  cfg.save_dir '\' save_name '\n\n'])
    saveas(gcf, [cfg.save_dir '\' save_name])
    saveas(gcf, [cfg.save_dir '\' save_name], 'png')
    saveas_eps(save_name,[cfg.save_dir '\'])
end
close all
%% print all the stats
fprintf('\n\n\nGamma Event Stats\n')
fileID = fopen(cat(2,PARAMS.stats_dir,'\Naris_stats_events.txt'),'w');
fprintf(fileID,'Gamma Event Stats\n')
fprintf(fileID, ['_________________________________________\n'])
fprintf(fileID, [date '\n'])

for iband = 1:length(bands)
%     for iSite =
        all_stats.(bands{iband}).avg_len = nanmean(all_stats.(bands{iband}).length);
        all_stats.(bands{iband}).std_len = nanstd(all_stats.(bands{iband}).length);
        all_stats.(bands{iband}).total = length(all_stats.(bands{iband}).length);
        all_stats.(bands{iband}).avg_rate = nanmean(all_stats.(bands{iband}).rate); % divide by the number of minutes in the session.
        all_stats.(bands{iband}).avg_nCycle = nanmean(all_stats.(bands{iband}).nCycles);
        fprintf(fileID, [bands{iband}  ':       mean length: ' num2str(all_stats.(bands{iband}).avg_len*1000)...
            ' std: ' num2str(all_stats.(bands{iband}).std_len*1000) ...
            ' nCycle: ' num2str(all_stats.(bands{iband}).avg_nCycle)...
            ' rate: ' num2str(all_stats.(bands{iband}).avg_rate)...
            ' total: ' num2str(all_stats.(bands{iband}).total) '\n']);
%     end
    end

    fclose(fileID);
end
