%% Load the data
M = load('/Users/kuitang/Datasets/RealityMining/realitymining.mat');

%% Initial bookkeeping
% First, gather the indices of the 94 who completed the survey. Gather
% their start and end dates too.

participants = zeros(94, 1);
startDates   = zeros(94, 1);
endDates     = zeros(94, 1);
activeTime   = zeros(94, 1);

ip = 1;
for i = 1:length(M.s)
    if ~isempty(M.s(i).surveydata)
        assert(ip <= 94, 'Only 94 subjects in the study...');
        
        participants(ip) = i;
        startDates(ip)   = datenum(M.s(i).my_startdate);
        
        if isempty(M.s(i).my_enddate)
            endDates(ip) = datenum('1/1/2007'); % sentinel
        else        
            endDates(ip) = M.s(i).my_enddate;
        end
        
        activeTime(ip) = endDates(ip) - startDates(ip);
        % Some may be negative... that seems like a data error.
%         assert(activeTime(ip) > 0);
        
        ip = ip + 1;      
    end
end

figure;
hist(activeTime, 50);
title('#Active days across users');

% We now refer to participants an index from 1:94, looking up s through the
% participants vector. For each of the 94 participants, map their mac
% address to their index.

macToIx = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:length(participants)
    macToIx(M.s(i).mac) = i;
end

% Check that the mac mapping is a bijection.
macToIxRange = cell2mat(macToIx.values());
assert(all(unique(macToIxRange) == 1:94));

%% Reconstruct the proximity log.
% local  = the user from whose phone I read the data.
% remote = the user to which the data entry refers.
dates  = [];
local  = [];
remote = [];

macsNotFound = [];
iMacsNotFound = 1;

iLog = 1;

date1Jan = datenum('1/1/2014');
nIgnoredDates = 0;

for iLocal = 1:length(participants)    
    ix = participants(iLocal);
    nEvents = length(M.s(ix).device_date);
    
    fprintf('Processing local user %d with %d events...\n', ...
        iLocal, nEvents);
    
    for n = 1:nEvents
        thisDate = M.s(ix).device_date(n);
        thisMacs = M.s(ix).device_macs{n};
        
        if thisDate == 0
            assert(isempty(thisMacs), 'Scan with zero date but nonempty mac.');
            continue;
        end
        
        if thisDate == date1Jan
            fprintf('Throwing away reset event...');
            nIgnoredDates = nIgnoredDates + 1;
            continue;
        end
        
        assert(~isempty(thisMacs), 'Scan with nonzero date but empty mac.');
        
        for m = 1:length(thisMacs)
            % Check that the mac exists. Sometimes proximity events yield
            % mac addresses that are not there...
            if macToIx.isKey(thisMacs(m))
                iRemote = macToIx(thisMacs(m));

                dates(iLog)  = thisDate;
                local(iLog)  = iLocal;
                remote(iLog) = iRemote;

                iLog = iLog + 1;
            else
                macsNotFound(iMacsNotFound) = thisMacs(m);
                iMacsNotFound = iMacsNotFound + 1;
            end            
        end
    end
end

fprintf('%d proximity events were logged, but %d macs were not found.\n', ...
    iLog - 1, iMacsNotFound - 1);

proximityLog = table(dates, local, remote);

%% checkpoint
save -v7.3 cp/realitymine_proximity_log.mat proximityLog participants startDates endDates activeTime

%% Construct daily networks out of the proximity log.

% Problem: Not every user present every day.
% Find earliest day everybody in the network is present.
% - Ignore 1 
% Find the earliest day
