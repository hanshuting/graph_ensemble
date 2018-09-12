function [M, Fk, population, matchedPeople] = parseRoommateMatch(year)
% parseRoommateMatch  Return perfect unipartite matching of matched people
%
%   [Msub, matchedPeople] = parseRoommateMatch(matchFile)
%   Msub          -- Perfect matching subgraph
%   matchedPeople -- ordered list of people in the matching
%
%   TODO: Data cleanup pages should output "2010Unmatched.txt" or
%   something.

    entryIds = importdata('data/flatEntryIds.txt');    
    nUsers   = length(entryIds);
    
    Mij = importdata(sprintf('data/%dMatch.txt', year));
    Mall = sparse(Mij(:,1), Mij(:,2), ones(size(Mij, 1), 1), nUsers, nUsers);        
    % symmetrize
    Mall = Mall + Mall';
    
    % Remember to convert 0 (NumPy) to 1 (MATLAB) indexing!
    population = importdata(sprintf('data/%dPopulation.txt', year)) + 1;
    
    matchedPeople = find(sum(Mall, 2));
    % Indices should be sorted
    assert(all(vec(matchedPeople == sort(matchedPeople))));
    
    M = Mall(matchedPeople, matchedPeople);
    assert(all(sum(M, 1) == 1) && all(sum(M, 2)) == 1, 'Msub is not a perfect matching.');
    
    Fk = roommateFeatures('SurveyData.txt', population);    
end
