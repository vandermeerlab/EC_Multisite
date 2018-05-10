function MS_collect_phase()
%% MS_collect_phase: consolidates all of the phase measures for each
%      subject into the same data sctruature as the individual subjects.
%      The output can then be used in the 'MS_plot_phase" script to get all
%      the plots.
%
%
% INPUTS:
%
%
%
% OUTPUTS:
%
%
%
%

%% for debugging
clear all
run('/Users/jericcarmichael/Documents/GitHub/EC_Multisite/MS_initialize.m')
global PARAMS
mat_out =[];
Sub_count = 1;
for Subjects = {'R102','R104','R107', 'R112','R122','R123'}
    load([PARAMS.inter_dir '/Phase_outputs/' Subjects{1} '_phase_out.mat']);
    
    if Sub_count ==1
        mat_out = Phase_mat; 
        
    else
        f_names = fieldnames(Phase_mat);
        for Fn = 1:length(f_names)
         disp(f_names{Fn})
         mat_out.(f_names{Fn}) = Phase_mat.(f_names{Fn});
        end
    end
    Sub_count = Sub_count+1;
    
    clear Phase_mat
end
%% save the output
 save([PARAMS.inter_dir '/Phase_outputs/all_phase_out'], 'mat_out', '-v7.3')
 
 %%

