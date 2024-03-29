% Multi_Master_preprocess:
%   Master control script for multi-site loading, analysis, and statistics.
%   Requirements:
%        - van der Meer lab codebase
%        - EC_Multisite functions (vandermeerlab/EC_Multisite on github)
%


% terminology:
%
%    - session: recording day which included four phases
%        - phase : recording phases of the naris protocol as "pre", "ipsi",
%                  "contra", "post"
%

%% make a log
fprintf(PARAMS.log, date);
% Extract the data from each recroding phase within each session and separate pot vs track sections
global PARAMS


for iSub = 4:length(PARAMS.Subjects)
   %% speed loading
%    load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'])
%    d_t = data;
%    clear data
%    data.(PARAMS.Subjects{iSub}) = d_t;

%% 
    if isunix
        cd([PARAMS.data_dir  PARAMS.Subjects{iSub}])
    else
        cd([PARAMS.data_dir  PARAMS.Subjects{iSub}])
    end

    dir_files = dir(); % get all the sessions for the current subject
    dir_files(1:2) = [];
    sess_list = [];
    for iDir = 1:length(dir_files)
        if dir_files(iDir).isdir == 1 && strcmp(dir_files(iDir).name(1:4), PARAMS.Subjects{iSub}) % ensure that the session folders have the correct names
            sess_list = [sess_list;dir_files(iDir).name];  % extract only the folders for the seesions
        end
    end
    sess_list = cellstr(sess_list);
    % load the data for each session within the current subject.
    for iSess = 1:length(sess_list)
        cfg_loading = [];
        cfg_loading.fname = sess_list{iSess};
        fprintf(['\n' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        fprintf(PARAMS.log,['\nLoading ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);

        [data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')), cfg_loading] = MS_load_data_fast(cfg_loading);
        fprintf(PARAMS.log, '...complete');
    end
    % ensure the correct number of sessions exist per rat
    if length(fieldnames(data.(PARAMS.Subjects{iSub}))) ~=4
        error('too many or too few sessions for multisite experiment.  Should only contain 4 per rat')
    end

% % 
%% get the gamma events per recording phase
fprintf(PARAMS.log,'\n\nCollecting Events');

    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess  = 1:length(sess_list)
        fprintf(PARAMS.log,['\nEvents ' PARAMS.Subjects{iSub} '  ' sess_list{iSess}]);
        cfg = []; cfg.f_bandpass = {[40 55],[70 100]}; 
        [Events.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_'))] = MS_extract_gamma([],data.(PARAMS.Subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')));
        fprintf(PARAMS.log, '...complete');
    end
    
    
% end % subjects


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% EVENT PROCESSING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Event Processing
% for iSub = 4%:length(PARAMS.Subjects)
%     load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Events.mat'])
%      load([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'])
% 	d_t = data;
    
    sess_list = fieldnames(data.(PARAMS.Subjects{iSub}));
    for iSess = 1:length(sess_list)
        fprintf(['\nRunning Phases Analyses on ' sess_list{iSess} '....\n']);
        
        cfg_in = [];
        cfg_in.Subject_id = sess_list{iSess};
        %	cfg.check = 1;
        
        Phase_mat.(sess_list{iSess}) =  MS_Phase_Analyses(cfg_in, data.(PARAMS.Subjects{iSub}).(sess_list{iSess}), Events.(PARAMS.Subjects{iSub}).(sess_list{iSess}));
        
        fprintf('...Complete\n');
    end
    
    fprintf('/nSaving Phase_mat output...');
    mkdir([PARAMS.inter_dir '/Phase_outputs/']); 
    save([PARAMS.inter_dir '/Phase_outputs/' (PARAMS.Subjects{iSub}) '_phase_out.mat'], 'Phase_mat', '-v7.3')
    fprintf('...Complete');


    
% end
%% save the intermediate files
fprintf(PARAMS.log,'\n\nSaving intermediates');
mkdir(PARAMS.data_dir, 'temp');
data = data.(PARAMS.Subjects{iSub}); 
save([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Data.mat'], 'data', '-v7.3')
save([PARAMS.inter_dir PARAMS.Subjects{iSub} '_Events.mat'], 'Events', '-v7.3')

clearvars -except iSub PARAMS
end

