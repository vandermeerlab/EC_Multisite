function [pow, coh] = Quick_psd()

%% add the vdmlab codebase.  this function calls LoadExpKeys, Linspecer, 
clear all
% addpath(genpath('C:\Users\admin\Documents\GitHub\vandermeerlab\code-matlab\shared'))
%%
LoadExpKeys
line_width = 2;
%% Loop through channels specified in the ExpKeys and get the PSD.  
sites = ExpKeys.Chan_to_use;
cfg = [];
cfg.fc = sites;
csc = LoadCSC(cfg);
% get the psd
cfg_psd = [];
cfg_psd.hann_win = 2^16; % always make this in base 2 for speed
    cut = round(length(csc.data(1,:))/2); % this can be used for speed.
%     only takes in half the data.  
for iSite = 1:length(sites)
    [pow.(sites{iSite}(1:end-4)).pxx, pow.(sites{iSite}(1:end-4)).f] = pwelch(csc.data(iSite,:), hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*4 , csc.cfg.hdr{1}.SamplingFrequency);
    real_label{iSite} = ExpKeys.Chan_to_use_labels{iSite}; % get the matching channel label. 
end

%% Get the coherence between pairs of sites. 
if length(sites) >1
    site_comb = nchoosek(sites, 2); % get all combinations of sites
    labels = [];
    for iComb = 1:length(site_comb); % loop through combinations of sites.  
        label = [site_comb{iComb, 1}(1:end-4) '_' site_comb{iComb, 2}(1:end-4)];
        S1 = strfind(sites, site_comb{iComb,1}); % find the corresponding
        S2 = strfind(sites, site_comb{iComb,2}); % find the other corrsponding site idx
        S1 = find(not(cellfun('isempty', S1))); % gets the actual idx
        S2 = find(not(cellfun('isempty', S2)));
        % get teh coherence using the same paramters as in pwelch.  
        [coh.(label).p, coh.(label).f] =mscohere(csc.data(S1,:), csc.data(S2,:),hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*4 , csc.cfg.hdr{1}.SamplingFrequency);
        labels{iComb} = label;
        pair_labels{iComb} = [real_label{S1}(1:3) '-' real_label{S2}(1:3)] % get the name of the pair.  
    end
end

%     %%
%
%
% for iSite = 1:length(sites)-1;
%     label = [sites{iSite}(1:end-4) '_' sites{iSite+1}(1:end-4)];
%     [coh.(label).p, coh.(label).f] =mscohere(csc.data(iSite,1:cut), csc.data(iSite+1,1:cut),hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*4 , csc.cfg.hdr{1}.SamplingFrequency);
%     labels{iSite} = label;
%     pair_labels{iSite} = [real_label{iSite}(1:3) '-' real_label{iSite+1}(1:3)]
% end
%% make a plot
figure
c_ord = linspecer(length(sites)+length(site_comb)); % predefine some colours. 

% electrode plot (beta)
if isfield(ExpKeys, 'Target')
    subplot(2,3,1)
    targets = fieldnames(ExpKeys.Target);
    hold on
    for iT = 1:length(targets)
        flag = 0;
        for iLab = 1:length(ExpKeys.Chan_to_use_labels)
            if strcmp(ExpKeys.Chan_to_use_labels{iLab}(1:3), targets{iT})
                flag = 1;
            end
        end
        if flag ==1
            rectangle('Position',[ExpKeys.Target.(targets{iT})(2)-.1, ExpKeys.Target.(targets{iT})(1)-.25, .5, .5],'Curvature',[1 1], 'facecolor', [c_ord(iT,:), .3])
        else
            rectangle('Position',[ExpKeys.Target.(targets{iT})(2)-.1, ExpKeys.Target.(targets{iT})(1)-.25, .5, .5],'Curvature',[1 1])
        end
        text(ExpKeys.Target.(targets{iT})(2), ExpKeys.Target.(targets{iT})(1),-ExpKeys.Depths.(targets{iT})(1), targets{iT})
    end
    text(0, 0, 'Bregma')
    xlim([-2 4])
    ylim([-2 2])
    zlim([-10 0])
end

%% plot the position data.note fails if there is more than one ..nvt file (can happen in RR2)
subplot(2,3,4)
if exist('VT1.zip') % unzips the nvt if it is zipped.  
    unzip('VT1.zip')
end
cfg_pos = [];
pos = LoadPos(cfg_pos);
plot(pos.data(1,:), pos.data(2,:), '.')
axis off
%% plot the PSD for each channel on one plot.  
subplot(2,3,[1,2,3])
for iSite = 1:length(sites)
    hold on
    plot(pow.(sites{iSite}(1:end-4)).f, 10*log10(pow.(sites{iSite}(1:end-4)).pxx), 'color', c_ord(iSite,:),'linewidth', line_width);
end
xlim([0 120])
y_val = ylim;
rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [0 0.2 .8 0.2],'edgecolor', [0 0.2 .8 0.2]); %low gamma rectangle
rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [0 0.8 .2 0.2],'edgecolor', [0 0.8 .2 0.2]) %high gamma rectangle

legend(strrep(ExpKeys.Chan_to_use_labels, '_', ' '))%,  'location', 'EastOutside', 'orientation', 'vertical');




%% plot the coherence between pairs of sites. 
subplot(2,3,[4,5, 6])
for iPairs = 1:length(labels)
    hold on
    plot(coh.(labels{iPairs}).f, coh.(labels{iPairs}).p, 'color', c_ord(length(sites)+iPairs,:), 'linewidth', line_width)
end
xlim([0 120])
ylim([0 1])
legend(pair_labels)%, 'location', 'EastOutside', 'orientation', 'vertical')
maximize
d_name = strsplit(cd, '\'); % what to name the figure.  
title(strrep(d_name{end}, '_', ' ')) 
