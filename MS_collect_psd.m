function MS_data_out = MS_collect_psd(cfg_in, MS_data_in)
%%MS_collect_psd: Collects the psd values for each phase and recording
%    site.  
%
%Inputs:
%   cfg [struct]: can be blank if you want to use the defaults.
%   MS_data_in [struct]: should be in the form of a csc from LoadCSC

%Outputs:
%   out [struct]: collection of the PSD appended to the data section for
%   each phase and site split into track and pot. 

% EC 2017-05-24

%% cycle through each phase/site and collect the PSDs

phases = fieldnames(MS_data_in); 
for iPhase = 1:length(phases)
    temp = MS_data_in.(phases{iPhase}); 
    temp = rmfield(temp, 'pos'); temp = rmfield(temp, 'ExpKeys'); % just to remove position and Expkeys 
    sites = fieldnames(temp); clear temp; 
    
    % get the psd for each site. 
    for iSite = 1:length(sites)
       fprintf('\n<strong>%s</strong> Processing Phase: %s Site: %s',mfilename, phases{iPhase}, sites{iSite})
       MS_data_out.(phases{iPhase}).(sites{iSite}).psd = MS_get_psd([], MS_data_in.(phases{iPhase}).(sites{iSite})); 
    end
    
end
fprintf('\n')