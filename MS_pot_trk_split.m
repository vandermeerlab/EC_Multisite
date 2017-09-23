function [Naris_pot, Naris_trk] = MS_pot_trk_split(Naris_in)
%% MS_pot_trk_split: separates 
global PARAMS

subjects = fieldnames(Naris_in);
for iSub = 1:length(subjects)
 sess_list = fieldnames(Naris_in.(subjects{iSub}));
 for iSess  = 1:length(sess_list)
       for iPhase = 1:length(PARAMS.Phases);
           data_list = fieldnames(Naris_in.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}));
           for iData = 1:length(data_list)
               if strfind(data_list{iData}, 'pot')
                   Naris_pot.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData}) = Naris_in.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData});
               elseif strfind(data_list{iData}, 'trk')
                   Naris_trk.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData}) = Naris_in.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData});
               elseif strfind(data_list{iData}, 'ExpKeys')
                   Naris_pot.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData}) = Naris_in.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData});
                   Naris_trk.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData}) = Naris_in.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData});
               elseif strfind(data_list{iData}, 'pos')
                   Naris_pot.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData}) = Naris_in.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData});
                   Naris_trk.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData}) = Naris_in.(subjects{iSub}).(strrep(sess_list{iSess}, '-', '_')).(PARAMS.Phases{iPhase}).(data_list{iData});
               end
           end
       end
       
   end
 end