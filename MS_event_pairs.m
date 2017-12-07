function [evts_out, mat_out] = MS_event_pairs(cfg_in, evts_in, mat)
%% MS_event_pairs: uses the detected events for each channel and extracts the IVs for the corresponding pair sites
%
%
%
%
%
%
%% default parameters


cfg_def = [];
cfg_def.wsize = 2^12; % keep in base 2 for speed. used for the entire session.  events are done using the length of the data.
cfg_def.f = [45 65; 70 90];
cfg = ProcessConfig2(cfg_def, cfg_in);

global PARAMS

%% dumb step to get all the permutations of the pairs of usable channels.
sites = mat.pre.ExpKeys.Chan_to_use_labels;
numS = 1:4;
pairs = {};
loop_num = 0;
for iSite = 1:length(sites)
    %     all_sites{iSite} = sites{iSite}(1:end-1);
    others = numS(numS~=iSite);
    for iOther = others
        loop_num = loop_num+1;
        pairs{loop_num} = [sites{iSite}(1:end-1) '_' sites{iOther}(1:end-1)];
    end
end

%% cycle through phases and recording cite pairs (uses only the "good pairs" specified in the ExpKeys)

types = {'_pot', '_trk'};
evts_out = evts_in;
phases = PARAMS.Phases;
bands = {'low', 'high'};

for iType = 1:length(types)
    for iPair = 1:length(pairs)
        sites = strsplit(pairs{iPair}, '_');
        if length(sites) >2
            newsites{1} = [sites{1} sites{2}];
            newsites{2} = sites{3};
            sites = newsites;
        end
        for iPhase = 1:length(PARAMS.Phases)
            for iBand = 1:length(bands)
                %                 evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).pairs = [];
                % collect the events using the current site
                if isfield(evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}), 'evts')
                    evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.pairs.(pairs{iPair}).(sites{1}) = evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts;
                    for iEvt = length(evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).tstart):-1:1
                        evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.pairs.(pairs{iPair}).(sites{2}){iEvt} = restrict(mat.(phases{iPhase}).([sites{2} types{iType}]), evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).tstart(iEvt),evts_in.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).tend(iEvt));
                    end
                else
                    evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.pairs.(pairs{iPair}).(sites{1}) = {};
                    evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.pairs.(pairs{iPair}).(sites{2}) = {};
                end
            end
        end
    end
end


%% compute the coherence for each gamma event for each pair of electrodes as specified in the ExpKeys.GoodPairs
% make a matrix for all the pairs
labels = mat.pre.ExpKeys.Chan_to_use_labels;
for ii = 1:length(labels)
    for jj = 1:length(labels)
        mat_out.labels{ii,jj} = [labels{ii}(1:end-1) '_' labels{jj}(1:end-1)];
    end
end
for iPhase = 1:4
    mat_out.(phases{iPhase}).sess_coh.low = NaN*zeros(size(mat_out.labels));
    mat_out.(phases{iPhase}).sess_coh.high = NaN*zeros(size(mat_out.labels));
    mat_out.(phases{iPhase}).evt_coh.low = NaN*zeros(size(mat_out.labels));
    mat_out.(phases{iPhase}).evt_coh.high = NaN*zeros(size(mat_out.labels));
    mat_out.(phases{iPhase}).sess_amp.low = NaN*zeros(size(mat_out.labels));
    mat_out.(phases{iPhase}).sess_amp.high = NaN*zeros(size(mat_out.labels));
    mat_out.(phases{iPhase}).evt_amp.low = NaN*zeros(size(mat_out.labels));
    mat_out.(phases{iPhase}).evt_amp.high = NaN*zeros(size(mat_out.labels));
end

%%
for iType = 1:length(types)
    for iPair = 1:length(pairs)
        sites = strsplit(pairs{iPair}, '_');
        if length(sites) >2
            newsites{1} = [sites{1} sites{2}];
            newsites{2} = sites{3};
            sites = newsites;
        end
        for iPhase = 1:length(PARAMS.Phases)
            for iBand = 1:length(bands)
                main = strsplit(pairs{iPair}, '_'); % identifies the "main" site for for coherence comparisons.
                
                %find the corresponding index in the mat_out matrix
                idx = strfind(mat_out.labels, pairs{iPair});
                [x_idx,y_idx] = find(not(cellfun('isempty', idx)));
                
                % set up the filter for the band of interest 
                                % set up the filter
                cfg_filter = [];
                cfg_filter.f = [cfg.f(iBand,1) cfg.f(iBand,2)];
                
                % filter the two signals
                d_f1 = FilterLFP(cfg_filter, mat.(phases{iPhase}).([sites{1} types{iType}]));
                d_f2 = FilterLFP(cfg_filter, mat.(phases{iPhase}).([sites{2} types{iType}]));
                %% get the coherence across the entire session for each channel
                [Coh_sess_temp.c, Coh_sess_temp.f] = mscohere(d_f1 .data,d_f2.data,hanning(round(cfg.wsize/2)),cfg.wsize/4,cfg.wsize,d_f1.cfg.hdr{1}.SamplingFrequency);
                
                f_idx = find(Coh_sess_temp.f > cfg.f(iBand,1) & Coh_sess_temp.f <= cfg.f(iBand,2));
                mat_out.(phases{iPhase}).sess_coh.(bands{iBand})(x_idx,y_idx) = mean(Coh_sess_temp.c(f_idx));
                clear Coh_sess_temp
                
                %% get the amplitude corr across all frequencies
               
                % get the envelope
                d_f1.data_env =  abs(hilbert(d_f1.data));
                d_f2.data_env =  abs(hilbert(d_f2.data));
                
                % get the amplitude coherence for the entire session (same parameters as phase coherence)
                [AMP_sess_temp.c, AMP_sess_temp.f] = mscohere(d_f1.data_env, d_f2.data_env,hanning(round(cfg.wsize/2)),cfg.wsize/4,cfg.wsize,d_f1.cfg.hdr{1}.SamplingFrequency);
                
                % find the values with the frequency range of interestet (cfg.f)
                f_idx = find(AMP_sess_temp.f > cfg.f(iBand,1) & AMP_sess_temp.f <= cfg.f(iBand,2));
                mat_out.(phases{iPhase}).sess_amp.(bands{iBand})(x_idx,y_idx) = mean(AMP_sess_temp.c(f_idx));
                clear d_f1 d_f2 AMP_sess_temp
                disp('Amp Cor')
                disp(mat_out.(phases{iPhase}).sess_amp.(bands{iBand}))
                %%    Get the coherence and amp corr on an event by event basis
                pairs_in = evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.pairs;  % keeps things a bit more clean
                if isempty(pairs_in.(pairs{iPair}).(main{1})) ~=1
                    for iEvt = length(pairs_in.(pairs{iPair}).(main{1})):-1:1
                        % get the coherence across all frequencies
                        win_size = length(pairs_in.(pairs{iPair}).(main{1}){iEvt}.data);
                        [Coh_evt_temp.c{iEvt}, Coh_evt_temp.f{iEvt}] = mscohere(pairs_in.(pairs{iPair}).(main{1}){iEvt}.data,pairs_in.(pairs{iPair}).(main{2}){iEvt}.data,...
                            hanning(round(win_size/2)),round(win_size/4),win_size,...
                            pairs_in.(pairs{iPair}).(main{1}){iEvt}.cfg.hdr{1}.SamplingFrequency);
                        %find the frequencies of interest and then
                        f_idx = find(Coh_evt_temp.f{iEvt} > cfg.f(iBand,1) & Coh_evt_temp.f{iEvt} <= cfg.f(iBand,2));
                        mat_out.(phases{iPhase}).evt_coh.(bands{iBand})(x_idx,y_idx,iEvt) = mean(Coh_evt_temp.c{iEvt}(f_idx));
                    end
                else
                    mat_out.(phases{iPhase}).evt_coh.(bands{iBand})(x_idx,y_idx,1) = NaN;
                end
                
                %% get the amplitude corr on an event by event basis
                if isempty(pairs_in.(pairs{iPair}).(main{1})) ~=1
                    %                     env_IV = TSDtoIV(
                    for iEvt = length(pairs_in.(pairs{iPair}).(main{1})):-1:1
                        % get the coherence across all frequencies
                        win_size = length(pairs_in.(pairs{iPair}).(main{1}){iEvt}.data);
                        [Coh_evt_temp.c{iEvt}, Coh_evt_temp.f{iEvt}] = mscohere(pairs_in.(pairs{iPair}).(main{1}){iEvt}.data,pairs_in.(pairs{iPair}).(main{2}){iEvt}.data,...
                            hanning(round(win_size/2)),round(win_size/4),win_size,...
                            pairs_in.(pairs{iPair}).(main{1}){iEvt}.cfg.hdr{1}.SamplingFrequency);
                        %find the frequencies of interest and then
                        f_idx = find(Coh_evt_temp.f{iEvt} > cfg.f(iBand,1) & Coh_evt_temp.f{iEvt} <= cfg.f(iBand,2));
                        mat_out.(phases{iPhase}).evt_amp.(bands{iBand})(x_idx,y_idx,iEvt) = mean(Coh_evt_temp.c{iEvt}(f_idx));
                    end
                else
                    mat_out.(phases{iPhase}).evt_amp.(bands{iBand})(x_idx,y_idx,1) = NaN;
                end
                mat_out.(phases{iPhase}).evt_amp.(bands{iBand})
                
                
                %%
                % collect output in a matrix that contains all the poissible pairs.
                disp([sites{1} ' ' types{iType} ' ' phases{iPhase} ' ' bands{iBand} ' ' pairs{iPair}])
                evts_out.([sites{1} types{iType}]).(phases{iPhase}).(bands{iBand}).evts.(pairs{iPair}).coh = mat_out;
                
            end
        end
    end
end
end


%% temp test
% figure
% hold on
% samples = [];
% for iSamp = 4:0.5:7
%     cfg.wsize = round(2^iSamp);
%     [Coh_temp.c{iEvt}, Coh_temp.f{iEvt}] = mscohere(pairs_in.(pairs{iPair}).(main{1}){iEvt}.data,pairs_in.(pairs{iPair}).(main{2}){iEvt}.data,hanning(cfg.wsize),cfg.wsize/2,2*cfg.wsize,pairs_in.(pairs{iPair}).(main{1}){iEvt}.cfg.hdr{1}.SamplingFrequency);
%
%     plot(Coh_temp.f{iEvt}, Coh_temp.c{iEvt})
%     samples = [samples, cfg.wsize];
% end
% xlim([0 120])
% legend(num2str(samples'))

%% temp Mat plot

% figure
% % mat = mat_out;
% for iPhase = 1:4
%     plot_mat = tril(mat_out.(phases{iPhase}).sess_amp.low,-1) +triu(mat_out.(phases{iPhase}).sess_amp.high,1);
%     s=size(plot_mat,1);
%     plot_mat(1:s+1:s*s) = NaN;
%
%     subplot(1,4,iPhase)
%     h = nan_imagesc_ec(plot_mat);
%     add_num_imagesc(h, plot_mat)
%     caxis([0 0.5])
%     Square_subplots
%     set(gca, 'xticklabel', (labels),'ytick', 1:length(mat.pre.ExpKeys.Chan_to_use_labels), 'yticklabel',(labels), 'xaxisLocation','top');
%     title(phases{iPhase})
% end

