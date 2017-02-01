function [] = startup_pmtk3()

    %% change to the directory storing this function, (should be PMTK3 root).
    w = which(mfilename()); 
    thisDir = fileparts(w);
    %cd(thisDir);
    addpath(thisDir);

    %% include directories

    include  = @(d)addpath(genpathPMTK(d)); 
    include(fullfile(thisDir, 'projects')); 
    include(fullfile(thisDir, 'pmtkTools')); 
    include(fullfile(thisDir, 'matlabTools')); 
    include(fullfile(thisDir, 'toolbox')); 
    include(fullfile(thisDir, 'demos')); 
    include(fullfile(thisDir, 'pmtkdataSmall')); 
    include(fullfile(thisDir, 'pmtkdataCopy'));             % may be initially empty
    include(fullfile(thisDir, 'pmtksupportCopy')); % may be initially empty
    if exist(fullfile(thisDir, 'docs', 'tutorial'), 'dir')
        include(fullfile(thisDir, 'docs', 'tutorial')); 
    end