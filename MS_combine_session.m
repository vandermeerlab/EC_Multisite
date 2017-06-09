clear all; close all;
in_fd{1} = 'G:\JK_recordings\Naris\R123\R123-2017-03-02\R123-2017-03-02_pre_pot';
in_fd{2} = 'G:\JK_recordings\Naris\R123\R123-2017-03-02\R123-2017-03-02_pre_trk';

out_fd = 'G:\JK_recordings\Naris\R123\R123-2017-03-02\R123-2017-03-02_pre';
out_fd_postfix = '';

%% get filenames
for iFD = 1:length(in_fd)
    
   cd(in_fd{iFD});
   tt{iFD} = sort(FindFiles('*.ntt'));
   csc{iFD} = sort(FindFiles('*.ncs'));
   st{iFD} = sort(FindFiles('*.nst'));
   evt{iFD} = sort(FindFiles('*.nev'));
   pos{iFD} = sort(FindFiles('*.nvt'));
end

%% match filenames
for iTT = 1:length(tt{1,1})
    
    fn1 = tt{1,1}{iTT};
    [~,fn1,~] = fileparts(fn1);
    
    for ifn2 = 1:length(tt{1,2})
        [~,fn2{ifn2},~] = fileparts(tt{1,2}{ifn2});
    end
    
    match_id = strmatch(fn1,fn2);
    
    if isempty(match_id)
        fprintf('No match found for %s...\n',fn1);
        continue;
    end
    
    % match found, do loading
    cd(in_fd{1});
    [Timestamps1, ScNumbers1, CellNumbers1, Features1, Samples1, Header1] = Nlx2MatSpike(tt{1,1}{iTT}, [1 1 1 1 1], 1, 1, [] );
    
    cd(in_fd{2});
    [Timestamps2, ScNumbers2, CellNumbers2, Features2, Samples2, Header2] = Nlx2MatSpike(tt{1,2}{match_id}, [1 1 1 1 1], 1, 1, [] );
    
    
    
    % merge
    Timestamps3 = cat(2,Timestamps1,Timestamps2);
    ScNumbers3 = cat(2,ScNumbers1,ScNumbers2);
    CellNumbers3 = cat(2,CellNumbers1,CellNumbers2);
    Features3 = cat(2,Features1,Features2);
    Samples3 = cat(3,Samples1,Samples2);
    
    % write
    out_fn = cat(2,fn1,out_fd_postfix,'.ntt');
    
    mkdir(out_fd)
    cd(out_fd);
    Mat2NlxSpike(out_fn, 0, 1, [], [1 1 1 1 1], Timestamps3, ScNumbers3, CellNumbers3, Features3, Samples3, Header1);
    
end

%% do the same for any stereotrodes ".nst"
for iTT = 1:length(st{1,1})
    
    fn1 = st{1,1}{iTT};
    [~,fn1,~] = fileparts(fn1);
    
    for ifn2 = 1:length(st{1,2})
        [~,fn2{ifn2},~] = fileparts(st{1,2}{ifn2});
    end
    
    match_id = strmatch(fn1,fn2);
    
    if isempty(match_id)
        fprintf('No match found for %s...\n',fn1);
        continue;
    end
    
    % match found, do loading
    cd(in_fd{1});
    [Timestamps1, ScNumbers1, CellNumbers1, Features1, Samples1, Header1] = Nlx2MatSpike(st{1,1}{iTT}, [1 1 1 1 1], 1, 1, [] );
    
    cd(in_fd{2});
    [Timestamps2, ScNumbers2, CellNumbers2, Features2, Samples2, Header2] = Nlx2MatSpike(st{1,2}{match_id}, [1 1 1 1 1], 1, 1, [] );
    
    
    
    % merge
    Timestamps3 = cat(2,Timestamps1,Timestamps2);
    ScNumbers3 = cat(2,ScNumbers1,ScNumbers2);
    CellNumbers3 = cat(2,CellNumbers1,CellNumbers2);
    Features3 = cat(2,Features1,Features2);
    Samples3 = cat(3,Samples1,Samples2);
    
    % write
    out_fn = cat(2,fn1,out_fd_postfix,'.nst');
    
    mkdir(out_fd)
    cd(out_fd);
    Mat2NlxSpike(out_fn, 0, 1, [], [1 1 1 1 1], Timestamps3, ScNumbers3, CellNumbers3, Features3, Samples3, Header1);
    
end

%% do the same for the CSCs
for iTT = 1:length(csc{1,1})
    
    fn1 = csc{1,1}{iTT};
    [~,fn1,~] = fileparts(fn1);
    
    for ifn2 = 1:length(csc{1,2})
        [~,fn2{ifn2},~] = fileparts(csc{1,2}{ifn2});
    end
    
    match_id = strmatch(fn1,fn2, 'exact');
    
    if isempty(match_id)
        fprintf('No match found for %s...\n',fn1);
        continue;
    end
    
    % match found, do loading
    cd(in_fd{1});
    [Timestamps1, ScNumbers1, SampleFreq1, nValidSamples1, Samples1, Header1] = Nlx2MatCSC(csc{1,1}{iTT}, [1 1 1 1 1], 1, 1, [] );
    
    cd(in_fd{2});
    [Timestamps2, ScNumbers2, SampleFreq2, nValidSamples2, Samples2, Header2] = Nlx2MatCSC(csc{1,2}{match_id}, [1 1 1 1 1], 1, 1, [] );
    
    
    
    % merge
    Timestamps3 = cat(2,Timestamps1,Timestamps2);
    ScNumbers3 = cat(2,ScNumbers1,ScNumbers2);
    SampleFreq3 = cat(2,SampleFreq1,SampleFreq2);
    nValidSamples3 = cat(2,nValidSamples1,nValidSamples2);
    Samples3 = cat(2,Samples1,Samples2);
    
    % write
    out_fn = cat(2,fn1,out_fd_postfix,'.ncs');
    
    mkdir(out_fd)
    cd(out_fd);
    Mat2NlxCSC(out_fn, 0, 1, [], [1 1 1 1 1 1], Timestamps3, ScNumbers3, SampleFreq3,  nValidSamples3, Samples3, Header1);
    
end

%% for the tracking data
for iTT = 1:length(pos{1,1})
    
    fn1 = pos{1,1}{iTT};
    [~,fn1,~] = fileparts(fn1);
    
    for ifn2 = 1:length(pos{1,2})
        [~,fn2{ifn2},~] = fileparts(pos{1,2}{ifn2});
    end
    
    match_id = strmatch(fn1,fn2, 'exact');
    
    if isempty(match_id)
        fprintf('No match found for %s...\n',fn1);
        continue;
    end
    
    % match found, do loading
    cd(in_fd{1});
    [Timestamps1, eX1,eY1, Angles1, Targets1, Points1, Header1] = Nlx2MatVT(pos{1,1}{iTT}, [1 1 1 1 1 1], 1, 1, [] );
    
    cd(in_fd{2});
    [Timestamps2, eX2,eY2, Angles2, Targets2, Points2, Header2] = Nlx2MatVT(pos{1,2}{match_id}, [1 1 1 1 1 1], 1, 1, [] );
    
    
    
    % merge
    Timestamps3 = cat(2,Timestamps1,Timestamps2);
    eX3 = cat(2,eX1,eX2);
    eY3 = cat(2,eY1,eY2);
    Angles3 = cat(2,Angles1,Angles2);
    Targets3 = cat(2,Targets1,Targets2);
    Points3 = cat(2,Points1,Points2);

    % write
    out_fn = cat(2,fn1,out_fd_postfix,'.nvt');
    
    mkdir(out_fd)
    cd(out_fd);
    Mat2NlxVT(out_fn, 0, 1, [], [1 1 1 1 1], Timestamps3, eX3,eY3, Angles3, Targets3, Points3, Header1);
    
end

%% for the events
for iTT = 1:length(evt{1,1})
    
    fn1 = evt{1,1}{iTT};
    [~,fn1,~] = fileparts(fn1);
    
    for ifn2 = 1:length(evt{1,2})
        [~,fn2{ifn2},~] = fileparts(evt{1,2}{ifn2});
    end
    
    match_id = strmatch(fn1,fn2, 'exact');
    
    if isempty(match_id)
        fprintf('No match found for %s...\n',fn1);
        continue;
    end
    
    % match found, do loading
    cd(in_fd{1});
    [Timestamps1, EventIDs1, TTLs1, Extras1, EventStrings1, Header1] = Nlx2MatEV(evt{1,1}{iTT}, [1 1 1 1 1], 1, 1, [] );
    
    % Multisite specific removal of the failed recording phase
    indx = strfind(EventStrings1,'Starting Recording'); 
    full_idx = find(~cellfun(@isempty,indx)); 
    
    if length(full_idx) > 1
        EventStrings1{full_idx(2)} = 'Failed Start Recording';
    end
    
        indx = strfind(EventStrings1,'Stopping Recording'); 
    full_idx = find(~cellfun(@isempty,indx)); 
    
    if length(full_idx) > 1
        EventStrings1{full_idx(2)} = 'Failed Stop Recording';
    end
    
    cd(in_fd{2});
    [Timestamps2, EventIDs2, TTLs2, Extras2, EventStrings2, Header2] = Nlx2MatEV(evt{1,2}{match_id}, [1 1 1 1 1], 1, 1, [] );

    % merge
    Timestamps3 = cat(2,Timestamps1,Timestamps2);
    EventIDs3 = cat(2,EventIDs1,EventIDs2);
    TTLs3 = cat(2,TTLs1,TTLs2);
    Extras3 = cat(2,Extras1,Extras2);
    EventStrings3 = [EventStrings1;EventStrings2];

    % write
    out_fn = cat(2,fn1,out_fd_postfix,'.nev');
    
    mkdir(out_fd)
    cd(out_fd);
    Mat2NlxEV(out_fn, 0, 1, [], [1 1 1 1 1], Timestamps3, EventIDs3, TTLs3, Extras3, EventStrings3, Header1);
end