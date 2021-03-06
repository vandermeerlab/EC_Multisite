function MS_amp_xcorr_session%(cfg_in, data)

%% temp
cfg_in = [];
global PARAMS

load([PARAMS.inter_dir 'R107_Data.mat'])
% data = data.R102_2016_09_24;
data = data.R107_2017_08_04;

% load([PARAMS.inter_dir 'R102_Events.mat'])

%% MS_amp_xcorr: gets the amplitude cross correlation for eahc event across a range of frequencies specified in cfg.freq in steps of cfg.freq_step (default is 1Hz) for a given pair of channels.
%
%
%   - amplitude x_cor: measure of the timing and correlation between the
%   power envelope of two signals.  Can reveal limited lead/lag
%   relationships
%
%
%
% Inputs:
%   - cfg_in: [struct]   contains any configuration paramters that deviated from the defaults.
%   - data:   [struct]   this is for a single subject for a single session.  data structure containing the tvec, data, labels, cfg
%   (with hdr) which results from MS_load_data_fast.
%   - events: [struct]   events IV structure for a single subject for a single
%   session. this will also contain the time intervals for all of the
%   events of interest.  In this case for the "low" and "high" gamma
%   events.
%
% Outputs:
%   phase_vals: [struct]
%
%%  Set up defaults
global PARAMS
cfg_def = [];
cfg_def.type = '_pot'; % can also be '_trk'.  For now only works with pot or track and not both.  to get both just run it again with the other cfg.type
cfg_def.check = 1; % used to create a plot for verification.
cfg_def.check_dir = [PARAMS.inter_dir 'phase_check']; %where to save the check figure.

cfg_def.cfg_filter = [];
% cfg_def.cfg_filter.freq = [3 100];
% cfg_def.cfg_filter.freq_step = 3;
cfg_def.cfg_filter.freq = [3 100];
cfg_def.cfg_filter.freq_step = 1;


%cfgs for amp x-corr
cfg_def.cfg_amp = [];
cfg_def.cfg_amp.count = 100;
cfg_def.cfg_amp.nShuffle = 100;

% cfgs for amp-corr-ogram
cfg_def.cfg_amp_cor_ogram = [];
cfg_def.cfg_amp_cor_ogram.dT = .1;

cfg = ProcessConfig2(cfg_def, cfg_in);
%% rename any "piri_x_..." data labels.
for iPhase = 1:length(PARAMS.Phases)
    f_names = fieldnames(data.(PARAMS.Phases{iPhase}));
    for iF = 1:length(f_names)
        if strcmp(f_names{iF}, 'Piri_O_pot')
            data.(PARAMS.Phases{iPhase}).PiriO_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_O_pot');
        elseif strcmp(f_names{iF}, 'Piri_O_trk')
            data.(PARAMS.Phases{iPhase}).PiriO_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_O_trk');
        elseif strcmp(f_names{iF}, 'Piri_N_pot')
            data.(PARAMS.Phases{iPhase}).PiriN_pot = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_N_pot');
        elseif strcmp(f_names{iF}, 'Piri_N_trk')
            data.(PARAMS.Phases{iPhase}).PiriN_trk = data.(PARAMS.Phases{iPhase}).(f_names{iF});
            data.(PARAMS.Phases{iPhase}) = rmfield(data.(PARAMS.Phases{iPhase}), 'Piri_N_trk');
        end
    end
end

%% get only the specific type.
site = fieldnames(data.post);
expStr = ['*' cfg.type];
regStr = ['^',strrep(strrep(expStr,'?','.'),'*','.{0,}'),'$'];
starts = regexpi(site, regStr);
iMatch = ~cellfun(@isempty, starts);
idx = find(iMatch);
site_list = site(idx);

pairs = {};
loop_num = 0;
numS =1:length(site_list);
for iSite = 1:length(site_list)
    others = numS(numS~=iSite);
    for iOther = others
        loop_num = loop_num+1;
        pairs{loop_num} = [site_list{iSite}(1:end-4) '_' site_list{iOther}(1:end-4)];
    end
end
% temporary
% pairs = {'PL_NAc', 'PL_CG', 'PL_PiriO','NAc_CG', 'NAc_PiriO', 'CG_PiriO'};



%% loop over pairs for each frequency band
            Freq_list = cfg.cfg_filter.freq(1): cfg.cfg_filter.freq_step:cfg.cfg_filter.freq(end);

for iPhase = 1:length(PARAMS.Phases)
        iPair = 5;
%    S = strsplit(pairs{iPair}, '_');
    completed_pairs = {};
%     for  iPair = 1:length(pairs)
        Same_flag = 0;
        %     iPair = 5;
        S = strsplit(pairs{iPair}, '_');
        % make sure oyu only have one copy of each pair.
        for iC = 1:length(completed_pairs)
            if strcmp(completed_pairs{iC}, [S{2} '_' S{1}])
                Same_flag = 1;
            end
        end
        
        if Same_flag ~=1
            % set up an NaN array
            
            times = data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).tvec(1):cfg.cfg_amp_cor_ogram.dT:data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).tvec(end);
            
            amp_ac.(pairs{iPair}).(PARAMS.Phases{iPhase}) =NaN(length(Freq_list), length(times));
            
            % S = strsplit(pairs{iPair}, '_');
            fprintf(['Processing: Frequency (Hz)  '])
            Amp.f.(pairs{iPair}).(PARAMS.Phases{iPhase}) = NaN(1,length(Freq_list));
            Amp.ac.(pairs{iPair}).(PARAMS.Phases{iPhase}) = NaN(1,length(Freq_list));
            Amp.lag.(pairs{iPair}).(PARAMS.Phases{iPhase}) = NaN(1,length(Freq_list));
            
            for iF = 1:length(Freq_list)
                
                this_F = Freq_list(iF);
                
                
                if iF>1
                    for j=0:log10(this_F-1)
                        fprintf('\b'); % delete previous counter display
                    end
                end
                fprintf('%d', this_F);
                %%
                cfg_filter= [];
                cfg_filter.display_filter = 0;
                cfg_filter.f = [this_F-1 this_F+1];
%                 cfg_filter.wp = cfg_filter.f*(2/data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).cfg.hdr{1}.SamplingFrequency);
%                 cfg_filter.ws = [cfg_filter.f(1)-2.5 cfg_filter.f(2)+2.5]*(2/data.(PARAMS.Phases{iPhase}).([S{1} '_pot']).cfg.hdr{1}.SamplingFrequency);
%                 cfg_filter.R = 0.5;
                cfg_filter.order = 4;%cheb2ord(cfg_filter.wp, cfg_filter.ws, cfg_filter.R,20);
                cfg_filter.verbose = 0;
                cfg_filter.type = 'fdesign';%'cheby1';

                FILT_1 = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).([S{1} '_pot']));
                %%
                FILT_2 = FilterLFP(cfg_filter, data.(PARAMS.Phases{iPhase}).([S{2} '_pot']));
                
                AMP_1 = FILT_1;
                AMP_2 = FILT_2;
                
                AMP_1.data = abs(hilbert(FILT_1.data));
                AMP_2.data = abs(hilbert(FILT_2.data));
                
                [ac, lag] = xcov(AMP_1.data,AMP_2.data,200,  'coeff');
                lag = lag * 1/AMP_1.cfg.hdr{1}.SamplingFrequency;
                [ac_max,idx] = max(ac);
                lag_max = lag(idx);
                
                Amp.f.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF) = this_F;
                Amp.ac.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF) = ac_max;
                Amp.lag.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF) = lag_max;
                
                
                %% create a Amp-xcorr-ogram
                
                for iT = 1:length(times)-1
                    
                    D_1_temp = restrict(AMP_1, times(iT), times(iT+1));
                    D_2_temp = restrict(AMP_2, times(iT), times(iT+1));
                    %                 Dxx1 = restrict(FILT_1, times(iT), times(iT+1)); checks for debugging
                    %                 Dxx2 = restrict(FILT_2, times(iT), times(iT+1));
                    
                    
                    [ac, lag] = xcov(D_1_temp.data,D_2_temp.data,200,  'coeff');
                    max_ac = ac(lag==0);
                    lag = lag * 1/AMP_1.cfg.hdr{1}.SamplingFrequency;
                    [m_ac,idx] = max(ac);
                    lag_max = lag(idx);
                    amp_lag.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF, iT) = lag_max;
                    amp_ac.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF, iT) = max_ac;
					                    amp_m_ac.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF, iT) = m_ac;

                    %             amp_ac_t.(pairs{iPair}).(PARAMS.Phases{iPhase})(iF, iT) = times(iT);
                    
                end
                amp_ac_t.(pairs{iPair}).(PARAMS.Phases{iPhase}) = times;
                
            end;
            fprintf(['\n' S{1} '-' S{2} '...complete\n'])
            completed_pairs{end+1} = pairs{iPair};
        else
            fprintf(['\n' 'Skipping ' S{1} '-' S{2} '  because that pair has already been processed\n'])
        end
%     end
end

%% check plot
% pair_list = fieldnames(Amp.ac);
% c_ord = linspecer(length(pair_list));
% for iPair = 1:length(pair_list)
%     
%     hold on
%     plot(Amp.f.(pair_list{iPair}).contra,Amp.ac.(pair_list{iPair}).contra, 'color', c_ord(iPair,:), 'linewidth', 3)
%     
% end
% for iPair = 1:length(pair_list)
%     
%     hold on
%     plot(Amp.f.(pair_list{iPair}).ipsi,Amp.ac.(pair_list{iPair}).ipsi, '--', 'color', c_ord(iPair,:), 'linewidth', 3)
%     
% end
% legend(pair_list)
% ylim([0 1])
% title('- - ipsi     - contra ')
% ylabel('Amplitude xcorr')
% xlabel('Frequency')
% SetFigure([], gcf)
%%
% for iPair = 1:length(pair_list)
%     
%     hold on
%     plot(Amp.f.(pair_list{iPair}).post,Amp.ac.(pair_list{iPair}).post, '-.', 'color',c_ord(iPair,:), 'linewidth', 3)
%     
% end
% 
% for iPair = 1:length(pair_list)
%     
%     hold on
%     plot(Amp.f.(pair_list{iPair}).pre,Amp.ac.(pair_list{iPair}).pre, '.', 'color',c_ord(iPair,:), 'linewidth', 3)
%     
% end
% title('-. pre - - ipsi - contra .. post')
% ylabel('Amplitude xcorr')
% xlabel('Frequency')
% SetFigure([], gcf)
% 
% 
%% try the amp-cor-ogram
% 
% figure(1111)
% imagesc(amp_ac.NAc_PiriO.pre)
% set(gca, 'yticklabel', Freq_list(2:end), 'xticklabel', amp_ac_t.NAc_PiriO.pre)
% axis xy

% %% append data
% amp_ac = amp_out.amp_ac;
% Freq_list = amp_out.freq;
% cfg.c_ord = linspecer(4);
% this_pair = 'OFC_NAc';
% all_sess = [amp_ac.(this_pair).pre, amp_ac.(this_pair).ipsi, amp_ac.(this_pair).contra, amp_ac.(this_pair).post];
% 
% 
% if cfg.moving_avg
%     if strcmp(cfg.moving_avg, '1D')
%         
%         for ii = size(all_sess, 1):-1:1
%             all_sess(ii,:) = conv(all_sess(ii,:), ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
%         end
%     end
%     if strcmp(cfg.moving_avg, '2D')
%         
%         all_sess = conv2(all_sess, ones(cfg.m_avg_nSample,1)/cfg.m_avg_nSample, 'same');
%     end
% end
% 
% 
% 
% 
% ax1 = imagesc(all_sess);
% axis xy
% set(gca, 'ytick',[1 size(all_sess,1)], 'yticklabel',[Freq_list(1) Freq_list(end)])%, 'xtick', 1:length(times), 'xticklabel', times(1):times(end))
% 
% 
% 
% set(ax1, 'AlphaData', ~isnan(all_sess))
%     %     title(['Coherogram: ' S{1} ' ' S{2}] )
%     xlabel('time (s)')
%     ylabel('frequency')
%     %         caxis([0.2 1])
%     
%     set(gca,'FontSize',20);
%     xlabel('time (s)'); ylabel('Frequency (Hz)');
%     yL = get(gca, 'ylim');
%     % add lines to distinguish phases
%     ylim([-4.5 yL(2)]);
%     hold on
%     h(1) =rectangle('position', [0 -4.5 length(amp_ac.(this_pair).pre)  5.5], 'facecolor',cfg.c_ord(1,:), 'edgecolor', cfg.c_ord(1,:));
%     h(2) =rectangle('position', [length(amp_ac.(this_pair).pre) -4.5 length(amp_ac.(this_pair).ipsi) 5.5], 'facecolor',cfg.c_ord(2,:), 'edgecolor', cfg.c_ord(2,:));
%     h(3) =rectangle('position', [(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)) -4.5 length(amp_ac.(this_pair).contra) 5.5], 'facecolor',cfg.c_ord(3,:), 'edgecolor', cfg.c_ord(3,:));
%     h(4) = rectangle('position', [(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra)) -4.5 length(amp_ac.(this_pair).post) 5.5], 'facecolor',cfg.c_ord(4,:), 'edgecolor', cfg.c_ord(4,:));
% 
%     text(length(amp_ac.(this_pair).pre)/2, -1.5, 'pre', 'fontsize', 20, 'horizontalAlignment', 'center')
%     text((length(amp_ac.(this_pair).pre)+(length(amp_ac.(this_pair).ipsi)/2)), -1.5, 'ipsi', 'fontsize', 20, 'horizontalAlignment', 'center')
%     text((length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+(length(amp_ac.(this_pair).contra)/2)), -1.5, 'contra', 'fontsize', 20, 'horizontalAlignment', 'center')
%     text((length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).ipsi)+(length(amp_ac.(this_pair).post)/2)), -1.5, 'post', 'fontsize', 20, 'horizontalAlignment', 'center')
% 
%     
% %     plot(Tp(1):Tp(end), zeros(length(Tp(1):Tp(end)),1)-off_set,'color', cfg.c_ord(1,:), 'linewidth', 1); % pre
% % %     plot(Tp(end):(Tp(end)+Ti(end)), zeros(length(Tp(end):(Tp(end)+Ti(end))),1)-off_set,'color', cfg.c_ord(2,:), 'linewidth', 10); % ipsi
% %     plot((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end)), zeros(length((Tp(end)+Ti(end)):(Tp(end)+Ti(end)+Tc(end))),1)-off_set,'color', cfg.c_ord(3,:), 'linewidth', 10); % contra
% %     plot((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end)), zeros(length((Tp(end)+Ti(end)+Tc(end)):(Tp(end)+Ti(end)+Tc(end)+Tpo(end))),1)-off_set,'color', cfg.c_ord(4,:), 'linewidth', 10); % post
% %     h = vline([length(amp_ac.(this_pair).pre),(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)),(length(amp_ac.(this_pair).pre)+length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).post))], {'k', 'k', 'k'} );
%     %     if cfg.lfp_on
%     %        plot(D1.tvec_norm,(D1.data*10000)+20, 'color', [1 1 1])
%     %        plot(D2.tvec_norm, (D2.data*10000)+20, 'color', [0.5 1 1])
%     %     end
%     ax =gca;
%     ax.Clipping = 'off';
%     set(h(:), 'linewidth', cfg.linewidth);
%     colormap('PARULA');
% 
% 
% 
% 
% 
% h = vline([length(amp_ac.(this_pair).pre), (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)),...
%     (length(amp_ac.(this_pair).pre) +length(amp_ac.(this_pair).ipsi)+length(amp_ac.(this_pair).contra))], {'k','k', 'k', 'k'}) ;
%         set(h(:), 'linewidth', 3)
%         SetFigure([], gcf)
%         
        
        
        
%%
amp_out.amp_ac = amp_ac;
amp_out.amp_lag = amp_lag;
amp_out.amp_ac_t = amp_ac_t;
amp_out.amp_m_ac = amp_m_ac;
amp_out.AMP = Amp;
amp_out.cfg = cfg;
amp_out.freq = Freq_list;

save([PARAMS.inter_dir 'R107_amp_p1Hz.mat'], 'amp_out','-v7.3')
