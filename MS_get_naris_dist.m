function out = MS_get_naris_dist(cfg_in, all_Naris)
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

% pl distances
Pl_dist = repmat([NaN,    4.841  NaN   NaN  4.327  3.976  NaN],4,1)';
IL_dist = repmat([3.124   NaN    NaN   NaN  NaN    NaN    3.561],4,1)';
OFC_dist = repmat([NaN    1.887  1.265 1    0.825  1.166  0.894],4,1)';
NAc_dist = repmat([0.447  1.897  0.6   1    1.414  0.6    1.4],4,1)';
CG_dist = repmat([NaN     6.251  NaN   NaN  5.855  5.492  6.030],4,1)';

R102 = repmat([NaN 3.124 NaN 0.447 NaN],4,1)'; 
R104 = repmat([ 4.841 NaN 1.887 1.897 6.251],4,1)';
R122 = repmat([NaN NaN 1.265 0.6 NaN],4,1)';
R123 = repmat([NaN NaN 1 1  NaN],4,1)';
R107 = repmat([4.327 NaN 0.825 1.414 5.855],4,1)';
R108 = repmat([3.967 NaN 1.166 0.6 5.492],4,1)';
R112 = repmat([NaN 3.561 0.894 1.4 6.030],4,1)'; 


distance = cat(3,R102, R104, R107, R108, R112, R122, R123)
% distance = cat(3, Pl_dist, IL_dist, OFC_dist, NAc_dist, CG_dist);


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
    all_low.(rec_type{iRec}).Contrast_Pxx = [];
    all_high.(rec_type{iRec}).Contrast_Pxx = [];
    all_low.(rec_type{iRec}).Contrast_White_Pxx = [];
    all_high.(rec_type{iRec}).Contrast_White_Pxx = [];
end

% collect PSDs
subjects = fieldnames(all_Naris);
for iRec= 1:length(rec_type)
    for iSub = 1:length(subjects)
        sess_list = fieldnames(all_Naris.(subjects{iSub}));
        this_low.Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        this_high.Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        this_low.White_Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        this_high.White_Pxx = NaN(length(sites), length(PARAMS.Phases)+1, length(sess_list));
        this_low.Contrast_Pxx = NaN(length(sites), length(sess_list));
        this_high.Contrast_Pxx = NaN(length(sites), length(sess_list));
        this_low.Contrast_White_Pxx = NaN(length(sites), length(sess_list));
        this_high.Contrast_White_Pxx = NaN(length(sites), length(sess_list));
        for iSess = 1:length(sess_list);
            for iSite = 1:length(sites)
                if sum(ismember(fieldnames(all_Naris.(subjects{iSub}).(sess_list{iSess}).pre), [sites{iSite} '_' rec_type{iRec}])) >0
                    for iPhase = 1:length(PARAMS.Phases)
                        for iType = 1:length(types)
                            current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess) = all_Naris.(subjects{iSub}).(sess_list{iSess}).(PARAMS.Phases{iPhase}).([sites{iSite} '_' rec_type{iRec}]).psd.(types{iType});
                            %                             all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}).(sites{iSite}) = cat(3,all_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType}).(sites{iSite}),current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).(types{iType})(iSite,:,iSess));
                        end
                        this_F = current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).F(iSite,:,iSess);
                        this_low.Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                        this_high.Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).Pxx(iSite,nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                        this_low.White_Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                        this_high.White_Pxx(iSite, iPhase, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).(PARAMS.Phases{iPhase}).White_Pxx(iSite,nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                    end
                    % contrast 
                    temp_contrast = (10*log10(current_PSD.(rec_type{iRec}).ipsi.Pxx(iSite,:,iSess))) - (10*log10(current_PSD.(rec_type{iRec}).contra.Pxx(iSite,:,iSess)));
                    this_low.Contrast_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F)));
                    this_high.Contrast_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F)));

                    temp_contrast = (10*log10(current_PSD.(rec_type{iRec}).ipsi.White_Pxx(iSite,:,iSess))) - (10*log10(current_PSD.(rec_type{iRec}).contra.White_Pxx(iSite,:,iSess)));
                    this_low.Contrast_White_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F)));
                    this_high.Contrast_White_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F)));

                    % set up control with average of pre and post
                    for iType = 1:length(types)
                        current_PSD.(rec_type{iRec}).control.(types{iType})(iSite,:,iSess) = mean([all_Naris.(subjects{iSub}).(sess_list{iSess}).pre.([sites{iSite} '_' rec_type{iRec}]).psd.(types{iType}),all_Naris.(subjects{iSub}).(sess_list{iSess}).post.([sites{iSite} '_' rec_type{iRec}]).psd.(types{iType})],2);
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
        
        % collect the contra-ipsi constrast mean for each subject/session/phase
%             for iSite = 1:length(sites)
%                 if ~all(isnan(this_low.Pxx(iSite, 2,:)))
%                     all_comp.(rec_type{iRec}).Pxx.low(iSite, iSub, :) = this_low.Pxx(iSite, 2,:)./this_low.Pxx(iSite, 3,:);
%                     all_comp.(rec_type{iRec}).Pxx.high(iSite, iSub, :) = this_high.Pxx(iSite, 2,:)./this_high.Pxx(iSite, 3,:);
%                     all_comp.(rec_type{iRec}).White_Pxx.low(iSite, iSub, :) = this_low.Pxx(iSite, 2,:)./this_low.White_Pxx(iSite, 3,:);
%                     all_comp.(rec_type{iRec}).White_Pxx.high(iSite, iSub, :) = this_high.White_Pxx(iSite, 2,:)./this_high.White_Pxx(iSite, 3,:);
% 
%                 else
%                     all_comp.(rec_type{iRec}).Pxx.low(iSite, iSub, 1:4) =  NaN(1,4);
%                     all_comp.(rec_type{iRec}).Pxx.high(iSite, iSub, 1:4) = NaN(1,4);
%                     all_comp.(rec_type{iRec}).White_Pxx.low(iSite, iSub, 1:4) = NaN(1,4);
%                     all_comp.(rec_type{iRec}).White_Pxx.high(iSite, iSub, 1:4) = NaN(1,4);
%                 end
%             end
        
        all_low.(rec_type{iRec}).Pxx = cat(3,all_low.(rec_type{iRec}).Pxx, this_low.Pxx);
        all_high.(rec_type{iRec}).Pxx = cat(3,all_high.(rec_type{iRec}).Pxx, this_high.Pxx);
        all_low.(rec_type{iRec}).White_Pxx = cat(3,all_low.(rec_type{iRec}).White_Pxx, this_low.White_Pxx);
        all_high.(rec_type{iRec}).White_Pxx = cat(3,all_high.(rec_type{iRec}).White_Pxx, this_high.White_Pxx);
        
        % same for the contrast comparisons
        all_low.(rec_type{iRec}).Contrast_Pxx = cat(3,all_low.(rec_type{iRec}).Contrast_Pxx, this_low.Contrast_Pxx);
        all_high.(rec_type{iRec}).Contrast_Pxx = cat(3,all_high.(rec_type{iRec}).Contrast_Pxx, this_high.Contrast_Pxx);
        all_low.(rec_type{iRec}).Contrast_White_Pxx = cat(3,all_low.(rec_type{iRec}).Contrast_White_Pxx, this_low.Contrast_White_Pxx);
        all_high.(rec_type{iRec}).Contrast_White_Pxx = cat(3,all_high.(rec_type{iRec}).Contrast_White_Pxx, this_high.Contrast_White_Pxx);
        
        clear current_PSD;
        
    end
end
all_low.labels = sites;
all_high.labels = sites;
%% match distance to each subject/site/session for the difference between contra and ipsi power
Subjects = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7'};
sites = {'PL', 'IL', 'OFC', 'NAc', 'CG'};

% for now remove the two piri rows.  format should be site x sess x subject
this_power_mat = all_low.pot.Contrast_Pxx;

this_power_mat(6,:,:) = [];
this_power_mat(4,:,:) = [];

%% apply the distance to piriform in a corresponding array. 
c_ord = linspecer(length(sites));
figure
hold on
for iSub = 1:length(Subjects)
    for iSite = 1:length(sites)
        if ~isnan(distance(iSite,1,iSub))
        scatter(distance(iSite,:,iSub), this_power_mat(iSite,:,iSub),50,c_ord(iSite,:), 'filled' )
        end
    end
end
legend(sites, 'location', 'southeast')










