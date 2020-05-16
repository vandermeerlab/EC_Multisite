% Multi-site spike sandbox
addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));

addpath(genpath('C:\Users\ecarm\Documents\GitHub\EC_Multisite'))
cd('M:\Multi_spike\Cat_Spikes')
% %%  try to cut them.
% restoredefaultpath
% addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\toolboxes\MClust-3.5'))
%
%
% %% try simple cut.
% restoredefaultpath
% addpath(genpath('C:\Users\ecarm\Documents\GitHub\simpleclust'))
% %% load the spikes
% cd('M:\Multi_spike\Cat_Spikes\R108-2017-08-01')
% addpath(genpath('/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'));
%%
cd('M:\Multi_spike\Cat_Spikes\R123_2017_02_26\')
cfg = [] ;
cfg.fc = {};
cfg.getTTnumbers = 0;
S = LoadSpikes(cfg);
firstSpike(S);
fprintf('Total recording time in S: %0.0fs\n', lastSpike(S) - firstSpike(S));
plot(S)

Evt = LoadEvents([]);

Starts = MS_extract_TS(Evt, 'Starting Recording');
Stops = MS_extract_TS(Evt, 'Stopping Recording');

for ii = 1:1:8
    disp((Stops.t{1}(ii)- Starts.t{1}(ii))/60)
end

%% [In reality this will be done to match the LFP data with a strict cut at 10mins] restrict to specific periods (need to grab cat events files for start/stops).
S_pre = restrict(S, Starts.t{1}(1), Stops.t{1}(1));

S_ipsi = restrict(S, Starts.t{1}(3), Stops.t{1}(3));

S_contra= restrict(S, Starts.t{1}(5), Stops.t{1}(5));

S_post = restrict(S, Starts.t{1}(7), Stops.t{1}(7));

%% get the mean firing rate and ISI

for iCell = length(S.t):-1:1
    ISI{iCell} = mean(diff(S.t{iCell})); 
    Rate{iCell} = 1/mean(diff(S.t{iCell}));
    % get mean ISI
%     ISI_pre{iCell} = median(diff(S_pre.t{iCell}));
%     ISI_ipsi{iCell} = median(diff(S_ipsi.t{iCell}));
%     ISI_contra{iCell} = median(diff(S_contra.t{iCell}));
%     ISI_post{iCell} = median(diff(S_post.t{iCell}));
%     
%     % get firing rate. 
%     Rate_pre{iCell} = 1./mean(diff(S_pre.t{iCell}));
%     Rate_ipsi{iCell} = 1./mean(diff(S_ipsi.t{iCell}));    %(length(S_ipsi.t{iCell})) ./((lastSpike(S_ipsi) - firstSpike(S_ipsi)))
%     Rate_contra{iCell} = 1./mean(diff(S_contra.t{iCell}));  %(length(S_contra.t{iCell})) ./((lastSpike(S_contra) - firstSpike(S_contra)))
%     Rate_post{iCell} = 1./mean(diff(S_post.t{iCell}));    %(length(S_post.t{iCell})) ./((lastSpike(S_post) - firstSpike(S_post)))
end



%% cross cor
cfg_ccf = [];
cfg_ccf.max_t = 0.2; 
for iCell = 1:length(S_pre.t)
    figure(iCell)
    sgtitle(['Cell: ' S_pre.label{iCell}])
    subplot(4,3,3)
    [ccf_val,tvec] = ccf(cfg_ccf,S_pre.t{iCell},S_pre.t{iCell});
    plot(tvec, ccf_val)
    text(0.05, 0.08, ['Rate: ' num2str(ISI_pre{iCell},'%2.2f') 'Hz']);
    text(0.05, 0.05, ['ISI: ' num2str(ISI_pre{iCell}*1000,'%2.0f') 'ms']);
    xlabel('pre')
    ylim([0 .1])
    
    subplot(4,3,6)
    [ccf_val,tvec] = ccf(cfg_ccf,S_ipsi.t{iCell},S_ipsi.t{iCell});
    plot(tvec, ccf_val)
    xlabel('ipsi')
    ylim([0 .1])
        text(0.05, 0.08, ['Rate: ' num2str(ISI_ipsi{iCell},'%2.2f') 'Hz']);
    text(0.05, 0.05, ['ISI: ' num2str(ISI_ipsi{iCell}*1000,'%2.0f') 'ms']);
    
    subplot(4,3,9)
    [ccf_val,tvec] = ccf(cfg_ccf,S_contra.t{iCell},S_contra.t{iCell});
    plot(tvec, ccf_val)
    xlabel('contra')
    ylim([0 .1])
        text(0.05, 0.08, ['Rate: ' num2str(ISI_contra{iCell},'%2.2f') 'Hz']);
    text(0.05, 0.05, ['ISI: ' num2str(ISI_contra{iCell}*1000,'%2.0f') 'ms']);
    
    subplot(4,3,12)
    [ccf_val,tvec] = ccf(cfg_ccf,S_post.t{iCell},S_post.t{iCell});
    plot(tvec, ccf_val)
    xlabel('post')
    ylim([0 .1])
        text(0.05, 0.08, ['Rate: ' num2str(ISI_post{iCell},'%2.2f') 'Hz']);
    text(0.05, 0.05, ['ISI: ' num2str(ISI_post{iCell}*1000,'%2.0f') 'ms']);
    
end
%%
spk_t = Pre_S.t{1};
max_t = 1;
binsize = 0.01;
xbin_centers = -max_t-binsize:binsize:max_t+binsize; % first and last bins are to be deleted later
ac = zeros(size(xbin_centers));

for iSpk = 1:length(spk_t)
    
    relative_spk_t = spk_t - spk_t(iSpk);
    
    ac = ac + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.
    
end

xbin = xbin_centers(2:end-1); % remove unwanted bins
ac = ac(2:end-1);

ac = ac./max(ac); % normalize

%%
spk_t = S.t{1}; % spike times
isi = diff(spk_t); % interspike intervals

dt = 0.001; % in s, because spike times are in s
isi_edges = 0:dt:0.25; % bin edges for ISI histogram
isi_centers = isi_edges(1:end-1)+dt/2; % for plotting

isih = histc(isi,isi_edges);

bar(isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
set(gca,'FontSize',20,'XLim',[0 0.25]); xlabel('ISI (s)'); ylabel('count'); grid on;