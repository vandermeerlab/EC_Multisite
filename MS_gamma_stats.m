function stats = MS_gamma_stats(cfg_in, Events)
%%
%
%
%
%
%



%% set up defaults

cfg_def = [];


cfg = ProcessConfig2(cfg_def, cfg_in)

global PARAMS

%% collect all the gamma event counts
out = [];
bands = {'low', 'high'};
phases = ['control', PARAMS.Phases];
sub_list = fieldnames(Events);
for iBand = 1:length(bands)
    for iSub = 1:length(sub_list)
        sess_list = fieldnames(Events.(sub_list{iSub}));
        for iSess = 1:length(sess_list);
            site_list = fieldnames(Events.(sub_list{iSub}).(sess_list{iSess}));
            for iSite = 1:length(site_list)
                for iPhase = 1:length(phases)
                    temp.(bands{iBand}).(site_list{iSite})(iSess, iPhase) = length(Events.(sub_list{iSub}).(sess_list{iSess}).(site_list{iSite}).(phases{iPhase}).(bands{iBand}).tstart);
                end
                temp.(bands{iBand}).(site_list{iSite})(iSess, :) = temp.(bands{iBand}).(site_list{iSite})(iSess, :)./ temp.(bands{iBand}).(site_list{iSite})(iSess, 1);
                if iSess == 1
                    out.(bands{iBand}).(site_list{iSite})(iSess, :) = temp.(bands{iBand}).(site_list{iSite})(iSess, :);
                else
                    out.(bands{iBand}).(site_list{iSite}) = [out.(bands{iBand}).(site_list{iSite}); temp.(bands{iBand}).(site_list{iSite})(iSess, :)];
                end
            end
            temp = [];
        end
    end
end

%% Normalize to the "control"

%% Stats and output
ks = [];
sites = fieldnames(out.low);
for iBand = 1:length(bands)
    for iSite = 1:length(sites)
        for ii = [3,4];
            h = kstest(out.(bands{iBand}).(sites{iSite})(:,ii));
            if h ~=1
                disp('***************************************************************')
                disp(['KS test FAIL for ' phases{ii}])
                disp('***************************************************************')
                ks.(bands{iBand}).(sites{iSite})(:,ii-2) = 1;
            else
                ks.(bands{iBand}).(sites{iSite})(:,ii-2) = 0;
            end
        end
    end
end
%%

% ipsi = reshape(comp_data.ipsi,1,numel(comp_data.ipsi));
% control = reshape(comp_data.control,1,numel(comp_data.control));
% contra = reshape(comp_data.contra,1,numel(comp_data.contra));
if sum(ks)>=1
    [p_ip_con, h_ip_con] = signrank(low_gamma.ipsi, low_gamma.contra);
    [p_ip_ctr, h_ip_ctr] = signrank(low_gamma.ipsi, low_gamma.control);
    [p_con_ctr, h_con_ctr] = signrank(low_gamma.contra, low_gamma.control);
    
    [h_p_ip_con, h_h_ip_con] = signrank(high_gamma.ipsi, high_gamma.contra);
    [h_p_ip_ctr, h_h_ip_ctr] = signrank(high_gamma.ipsi, high_gamma.control);
    [h_p_con_ctr, h_h_con_ctr] = signrank(high_gamma.contra, high_gamma.control);
else
    fprintf('\n\nUsing T-Test\n\n')
    [h_ip_con, p_ip_con, ~, l_stats_ip_con] = ttest(low_gamma.ipsi, low_gamma.contra);
    [h_ip_ctr, p_ip_ctr, ~,l_stats_ip_ctr] = ttest(low_gamma.ipsi, low_gamma.control);
    [h_con_ctr, p_con_ctr, ~,l_stats_con_ctr] = ttest(low_gamma.contra, low_gamma.control);
    
    [h_h_ip_con, h_p_ip_con, ~, h_stats_ip_con] = ttest(high_gamma.ipsi, high_gamma.contra);
    [h_h_ip_ctr, h_p_ip_ctr, ~,h_stats_ip_ctr] = ttest(high_gamma.ipsi, high_gamma.control);
    [h_h_con_ctr, h_p_con_ctr, ~,h_stats_con_ctr] = ttest(high_gamma.contra, high_gamma.control);
end
%low SEM
low_pre_SEM = std(low_gamma.pre)/(sqrt(length(low_gamma.pre)));
low_ipsi_SEM = std(low_gamma.ipsi)/(sqrt(length(low_gamma.ipsi)));
low_contra_SEM = std(low_gamma.contra)/(sqrt(length(low_gamma.contra)));
low_post_SEM = std(low_gamma.post)/(sqrt(length(low_gamma.post)));
low_control_SEM = std(low_gamma.control)/(sqrt(length(low_gamma.control)));
% high SEM
high_pre_SEM = std(high_gamma.pre)/(sqrt(length(high_gamma.pre)));
high_ipsi_SEM = std(high_gamma.ipsi)/(sqrt(length(high_gamma.ipsi)));
high_contra_SEM = std(high_gamma.contra)/(sqrt(length(high_gamma.contra)));
high_post_SEM = std(high_gamma.post)/(sqrt(length(high_gamma.post)));
high_control_SEM = std(high_gamma.control)/(sqrt(length(high_gamma.control)));



