%%
% This is an example of how to create a polar histogram in MATLAB&#174;.
% 
% Read about the <http://www.mathworks.com/help/matlab/ref/polarhistogram.html |polarhistogram|> function in the MATLAB documentation. This function is available in R2016b or newer.
%
% For more examples, go to <http://www.mathworks.com/discovery/gallery.html MATLAB Plot Gallery>
%
% Copyright 2017 The MathWorks, Inc.

% Check version
if verLessThan('matlab','9.1')
    error(['polarhistogram is available in R2016b or newer. ', ...
        'For older releases, use rose instead.'])
end

% Load sunspot data
load sunspotData sunspot

% Create a polar histogram with 24 sectors
figure
polarhistogram(sunspot, 24)

% Add title
title('Sunspot Frequency')
