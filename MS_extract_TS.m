function TS_out = MS_extract_TS(TS_in, label)
%% MS_get_NLX_events: pulls out a spefic label/cell from an IV.  Used for isolating a cell in a Spike IV or an event in an .nev IV.
%
%
%
%    Inputs:
%    - TS_in [struct]   TS 'timestamp' structure produced by function 'ts'
%
%    - label [n x str]  str or cell array with multiple labels.
%
%
%
%    Outputs:
%    - TS_out  [struct]  new TS containing only the labels found in TS_in
%
%
%
%
% EC 2020-05-13   initial version
%
%
%
%%  Find the desired labels

if ischar(label) % determine if the label input is a string or a cell array.
    %     fprintf('<strong> %s</strong>: label is a string\n', mfilename);
    
    l_idx = find(ismember(TS_in.label, label));
    
    TS_out = ts(TS_in.t(l_idx));
    TS_out.label = TS_in.label{l_idx};
    
elseif iscell(label)
    %         fprintf('<strong> %s</strong>: label is a cell\n', mfilename);
    
    %empty the TS_out to make it easy to fill in.
    TS_out = TS_in;
    TS_out.t =[];
    TS_out.label = {};
    
    l_idx = find(ismember(TS_in.label, label));
    
    for iL = length(l_idx):-1:1
        TS_out.t{iL} = TS_in.t{l_idx(iL)};
        TS_out.label{iL} = TS_in.label{l_idx(iL)};
    end
end

end % function.