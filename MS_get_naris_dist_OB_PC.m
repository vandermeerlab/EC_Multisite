function out = MS_get_naris_dist_OB_PC(cfg_in, all_Naris)
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
    Contra_low.(rec_type{iRec}).Pxx = [];
    Contra_high.(rec_type{iRec}).Pxx = [];
    Contra_low.(rec_type{iRec}).White_Pxx =[];
    Contra_high.(rec_type{iRec}).White_Pxx =[];
    Ipsi_low.(rec_type{iRec}).Pxx = [];
    Ipsi_high.(rec_type{iRec}).Pxx = [];
    Ipsi_low.(rec_type{iRec}).White_Pxx =[];
    Ipsi_high.(rec_type{iRec}).White_Pxx =[];
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
        %% contrast only
        just_contra_low.Pxx = NaN(length(sites), length(sess_list));
        just_contra_low.White_Pxx = NaN(length(sites), length(sess_list));
        just_contra_high.Pxx = NaN(length(sites), length(sess_list));
        just_contra_high.White_Pxx = NaN(length(sites), length(sess_list));
        
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
                    temp_contrast = ((10*log10(current_PSD.(rec_type{iRec}).ipsi.Pxx(iSite,:,iSess))) - (10*log10(current_PSD.(rec_type{iRec}).contra.Pxx(iSite,:,iSess))))./((10*log10(current_PSD.(rec_type{iRec}).ipsi.Pxx(iSite,:,iSess))) + (10*log10(current_PSD.(rec_type{iRec}).contra.Pxx(iSite,:,iSess))));
                    this_low.Contrast_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F)));
                    this_high.Contrast_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F)));
                    
                    just_contra_low.Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).contra.Pxx(iSite, nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                    just_contra_high.Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).contra.Pxx(iSite, nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                    just_ipsi_low.Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).ipsi.Pxx(iSite, nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                    just_ipsi_high.Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).ipsi.Pxx(iSite, nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));

                    
                    temp_contrast = ((10*log10(current_PSD.(rec_type{iRec}).ipsi.White_Pxx(iSite,:,iSess))) - (10*log10(current_PSD.(rec_type{iRec}).contra.White_Pxx(iSite,:,iSess))))./((10*log10(current_PSD.(rec_type{iRec}).ipsi.Pxx(iSite,:,iSess))) + (10*log10(current_PSD.(rec_type{iRec}).contra.Pxx(iSite,:,iSess))));
                    this_low.Contrast_White_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F)));
                    this_high.Contrast_White_Pxx(iSite, iSess) = mean(temp_contrast(nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F)));

                    just_contra_low.White_Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).contra.White_Pxx(iSite, nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                    just_contra_high.White_Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).contra.White_Pxx(iSite, nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));
                    just_ipsi_low.White_Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).ipsi.White_Pxx(iSite, nearest_idx(cfg.filter(1,1), this_F):nearest_idx(cfg.filter(1,2), this_F),iSess)));
                    just_ipsi_high.White_Pxx(iSite, iSess) = mean(10*log10(current_PSD.(rec_type{iRec}).ipsi.White_Pxx(iSite, nearest_idx(cfg.filter(2,1), this_F):nearest_idx(cfg.filter(2,2), this_F),iSess)));

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
        %% pure contra values. 
        Contra_low.(rec_type{iRec}).Pxx = cat(3,Contra_low.(rec_type{iRec}).Pxx, just_contra_low.Pxx);
        Contra_high.(rec_type{iRec}).Pxx = cat(3,Contra_high.(rec_type{iRec}).Pxx, just_contra_high.Pxx);
        Contra_low.(rec_type{iRec}).White_Pxx = cat(3,Contra_low.(rec_type{iRec}).White_Pxx, just_contra_low.White_Pxx);
        Contra_high.(rec_type{iRec}).White_Pxx = cat(3,Contra_high.(rec_type{iRec}).White_Pxx, just_contra_high.White_Pxx);

        %% pure ipsi values.
        Ipsi_low.(rec_type{iRec}).Pxx = cat(3,Contra_low.(rec_type{iRec}).Pxx, just_contra_low.Pxx);
        Ipsi_high.(rec_type{iRec}).Pxx = cat(3,Contra_high.(rec_type{iRec}).Pxx, just_contra_high.Pxx);
        Ipsi_low.(rec_type{iRec}).White_Pxx = cat(3,Contra_low.(rec_type{iRec}).White_Pxx, just_contra_low.White_Pxx);
        Ipsi_high.(rec_type{iRec}).White_Pxx = cat(3,Contra_high.(rec_type{iRec}).White_Pxx, just_contra_high.White_Pxx);
        
        clear current_PSD;
        
    end
end
all_low.labels = sites;
all_high.labels = sites;

%% match distance to each subject/site/session for the difference between contra and ipsi power
Subjects = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7'};
sites = {'PL', 'IL', 'OFC', 'NAc', 'CG'};

% for now remove the two piri rows.  format should be site x sess x subject
% this_power_mat = all_low.pot.Contrast_Pxx;
for iBand = {'low', 'high'}
    if strcmp(iBand, 'low')
        this_power_mat = Contra_low.pot.White_Pxx;
    elseif strcmp(iBand, 'high')
        this_power_mat = Contra_high.pot.White_Pxx;
    else
        error('no band selected')
    end
this_power_mat(6,:,:) = [];
this_power_mat(4,:,:) = [];

%% fix an issue where R102 has the PL and IL mislabeled.  
% this_power_mat(2,:,1) = this_power_mat(1,:,1);
% this_power_mat(1,:,1) = NaN;
% % temporary R102 CG fix
% this_power_mat(5,:,1) = NaN;
%% apply the distance to OB in a corresponding array. 
c_ord = linspecer(length(sites));
m_ord = {'o', '+', '*', 'x', 's', 'd', 'p'};
figure(100)
hold on
for iSub = 1:length(Subjects)
    for iSite = 1:length(sites)
        if ~isnan(distance_ob(iSite,1,iSub))
        plot(distance_ob(iSite,:,iSub), this_power_mat(iSite,:,iSub),m_ord{iSub}, 'color', c_ord(iSite,:), 'markersize', 10)
        end
    end
end
xlabel('Distance from OB (mm)')
ylabel('Ipsi/Contra contrast index')
% annoying forced legend.  Avoids issue of markers and colors not working properly. 
hold on
h = zeros(length(sites), 1);
for iSite = 1:length(sites)
    h(iSite) = plot(NaN,NaN,'color', c_ord(iSite,:));
end
[~, hobj, ~, ~] = legend(h, sites, 'location', 'northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',3);
% legend(sites, 'location', 'southeast')
SetFigure([], gcf)
%% apply the distance to piriform in a corresponding array. 
c_ord = linspecer(length(sites));
m_ord = {'o', '+', '*', 'x', 's', 'd', 'p'};
figure(101)
hold on
for iSub = 1:length(Subjects)
    for iSite = 1:length(sites)
        if ~isnan(distance_pc(iSite,1,iSub))
        plot(distance_pc(iSite,:,iSub), this_power_mat(iSite,:,iSub),m_ord{iSub}, 'color', c_ord(iSite,:), 'markersize', 10)
        end
    end
end
xlabel('Distance from PC (mm)')
ylabel('Ipsi/Contra contrast index')
% annoying forced legend.  Avoids issue of markers and colors not working properly. 
hold on
h = zeros(length(sites), 1);
for iSite = 1:length(sites)
    h(iSite) = plot(NaN,NaN,'color', c_ord(iSite,:));
end
[~, hobj, ~, ~] = legend(h, sites, 'location', 'northeast');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',3);
% legend(sites, 'location', 'southeast')
SetFigure([], gcf)
%% make a regression plot

% convert the data into a distance and power 1d array
dist_1d_ob = reshape(distance_ob, 140,1);
dist_1d_pc = reshape(distance_pc, 140,1);

pow_1d = reshape(this_power_mat, 140,1);
sess_1d = reshape(repmat([1:4]',5,1,7),140,1);
for ii = 1:length(Subjects)
rat_ids(:,:,ii) = ones(5,4)*ii;
end
rat_1d = reshape(rat_ids, 140,1);

%% 
% remove NaN values
nan_idx = isnan(dist_1d_pc); %first for any NaNs in the distances which correspond to missed electrodes per subject
pow_1d(nan_idx) = [];
dist_1d_ob(nan_idx) =[]; 
dist_1d_pc(nan_idx) =[]; 
sess_1d(nan_idx) = [];
rat_1d(nan_idx) = [];

nan_idx = isnan(pow_1d); % second for anythin in power.  Can correspond to missing sites that don't match the ExpKeys.  will be fixed with updated ExpKeys for a few subjects.
pow_1d(nan_idx) = [];
dist_1d_ob(nan_idx) =[]; 
dist_1d_pc(nan_idx) =[]; 
sess_1d(nan_idx) = [];
rat_1d(nan_idx) = [];

% figure
% p = polyfit(dist_1d, pow_1d, 1);
% hold on
%  plot(p)
%%
clear D_power
D_power.tbl = table(rat_1d, sess_1d,dist_1d_pc, dist_1d_ob, pow_1d,'VariableNames',{'RatID','SessID', 'Distance_pc','Distance_ob', 'Power'});
D_power.tbl.RatID = nominal(D_power.tbl.RatID);
D_power.tbl.SessID = nominal(D_power.tbl.SessID);
D_power.lme = fitlme(D_power.tbl,'Power~1+Distance_pc+(1|RatID)+(1|SessID)');
D_power.lme_ob = fitlme(D_power.tbl,'Power~1+Distance_ob+(1|RatID)+(1|SessID)');
% not used
D_power.lme_2 = fitlme(D_power.tbl,'Power~1+Distance_pc+(1+Distance_pc|RatID)+(1|SessID)');
D_power.lme_2b = fitlme(D_power.tbl,'Power~1+Distance_pc+(1+Distance_pc|RatID)+(1+Distance_pc|SessID)');
D_power.lme_3 = fitlme(D_power.tbl,'Power~1+Distance_pc^2+(1|RatID)+(1|SessID)');
D_power.lme_4 = fitlme(D_power.tbl,'Power~1+Distance_pc^2+(1+Distance_pc|RatID)+(1|SessID)');
D_power.lme_red = fitlme(D_power.tbl,'Power~1+(1|RatID)+(1|SessID)');
D_power.lme_red_fixed = fitlme(D_power.tbl,'Power~1+SessID +(RatID)');
% D_power.lme_red_sess_by_sub = fitlme(D_power.tbl,'Power~Distance +(1|RatID)+(SessID-1|RatID)');

figure
plotResiduals(D_power.lme,'fitted')

% APP.lme_reduced = fitlme(APP.tbl,'Power~(1|RatID)+(1|SessID)');
D_power.comparisonpc_ob = compare(D_power.lme,D_power.lme_ob);

D_power.comparison1_v2 = compare(D_power.lme,D_power.lme_2);
D_power.comparison1_v3 = compare(D_power.lme,D_power.lme_3);
D_power.comparison1_v4 = compare(D_power.lme,D_power.lme_4);
D_power.comparison1_vR = compare(D_power.lme_red, D_power.lme);
D_power.comparison2_v3 = compare(D_power.lme_2,D_power.lme_3);
D_power.comparison2_V4 = compare(D_power.lme_2,D_power.lme_4);
D_power.comparison2_vR = compare(D_power.lme_2,D_power.lme_red);
D_power.comparison3_v4 = compare(D_power.lme_3,D_power.lme_4);
D_power.comparison3_vR = compare(D_power.lme_3,D_power.lme_red);
D_power.comparison4_vR = compare(D_power.lme_4,D_power.lme_red);

%% write the output
if exist(['LME_' iBand{1} '_' datestr(date, 'YY_mm_dd') '.txt'], 'file');
    delete(['LME_' iBand{1} '_' datestr(date, 'YY_mm_dd') '.txt'])
end
clc
diary('on')
diary(['LME_' iBand{1} '_' datestr(date, 'YY_mm_dd') '.txt'])
disp([' LME Out PC' iBand{1}])
disp(D_power.lme)
disp(' LME Anova Out')
anova(D_power.lme)
disp('Compare out')
compare(D_power.lme_red, D_power.lme)
disp('     ')
disp([' LME Out for OB' iBand{1}])
disp(D_power.lme_ob)
disp(' LME Anova Out')
anova(D_power.lme_ob)
disp('Compare out OB vs red')
compare(D_power.lme_red, D_power.lme_ob)
disp('Compare out OB to PC')
compare(D_power.lme_ob, D_power.lme)
diary('off')
movefile(['LME_' iBand{1} '_' datestr(date, 'YY_mm_dd') '.txt'], PARAMS.stats_dir);
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


