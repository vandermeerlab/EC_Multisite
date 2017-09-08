function out = MS_get_subject_info(data);
%% MS_get_subject_info:  collects informationf rom the ExpKeys of each 
% subject used in the data structure.  Helpful for identifying any effect
% of electrode location.
%
%
% Inputs:
%    - data [struct] contains the CSC, position, hdr, and ExpKeys for each
%    subject used in the MS analyses
%
% Ouputs:
%    - out [struct] contains the targets as well as their depths as well as
%    weights for each subject along with averages across all subjects
%
%% collect the number of good electrodes and the depth of the electrode
% NOTE: for the MS experiment the electrodes were kept at a constant depth
% across all recording days. 