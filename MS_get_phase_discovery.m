function MS_get_phase_discovery(Subject)
%% run phase analyses on discovery
% get the phase information from MS_Phase_Analyses for a single subject. 
% these analyses are rather heavy so breaking it into subjects allows for
% simple parallel processing.  
%
% Inputs
%    Subject: [string] the subject ID: 'R102' or 'R112', ....
%
% outputs:
%  none.  the outputs from MS_Phase_Analyses is saved in a new folder in
%  the PARAMS.inter_dir called "Phase_outputs".  

MS_initialize_discovery

load([PARAMS.inter_dir Subject '_Data.mat']);
load([PARAMS.inter_dir Subject '_Events.mat']);

sess_list = fieldnames(data); 
for iSess = 1:length(sess_list)
   fprintf(['\nRunning Phases Analyses on ' sess_list{iSess} '....\n']);

	cfg_in = []
    cfg_in.Subject = Subject; 
%	cfg.check = 1;   
   Phase_mat.(sess_list{iSess}) =  MS_Phase_Analyses(cfg_in, data.(sess_list{iSess}), Events.(sess_list{iSess}));
    
   fprintf('...Complete\n');
end

fprintf('/nSaving Phase_mat output...');
save([PARAMS.inter_dir '/Phase_outputs/' Subject '_phase_out.mat'], 'Phase_mat', '-v7.3')
fprintf('...Complete');   Phase_mat.(sess_list{iSess}) =  MS_Phase_Analyses([], data.(sess_list{iSess}), Events.(sess_list{iSess}));

