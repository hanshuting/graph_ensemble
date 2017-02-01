%% New roommate features.

%% Setup
surveyData = importdata('data/SurveyData.txt');

years = 2010:2012;
M = length(years);

cellOfFeatures = cell(M, 1);
Ys = cell(M, 1);

%% (break to prevent for loop execution)
for m = 1:M 
    %% Just one year    
    matchFile      = sprintf('data/%dMatch.txt', years(m));
%     populationFile = sprintf('data/%dPopulation.txt', years(m));
        
    % IMPORTANT: matchSparse was stored on disk as 0-index. Convert to
    % 1-index.
    
    matchList = importdata(matchFile) + 1;
    maxPerson = max(matchList(:));
    % Sorted in ascending order
    population = sort(matchList(:));    
    nMatched = length(population);
    
    bigMatch    = sparse(matchList(:,1), matchList(:,2), ...
        ones(size(matchList, 1), 1), ...
        maxPerson, ...
        maxPerson);
    
    % Filter out only the matched individuals. First select the rows, and
    % then the columns corresponding to the matched population.
    Yrows  = bigMatch(population,:);
    Ys{m} = Yrows(:,population);    
    
    % Symmetrize
    Ys{m} = Ys{m} + Ys{m}';
    
    % Check that we have a perfect matching
    assert(all(sum(Ys{m}, 1) == ones(1, nMatched)) && ...
           all(sum(Ys{m}, 2) == ones(nMatched, 1)), ...
           'Filtered adjacency matrix was not a perfect matching!');
                        
    %% (break to prevent feature creation)
    cellOfFeatures{m} = [ rSquareDiffFeats(surveyData, population), ...
                          rBiasFeat([], population) ];
end

%% Make into MxK cell array
K = length(cellOfFeatures{1});

features = cell(M, K);
for m = 1:M
    features(m,:) = cellOfFeatures{m};
end

%% Plot adjacency matrices
% for m = 1:M
%     %% One iteration
%     figure;
%     spy(Ys{m});
%     title(sprintf('%d Roommate Matches', years(m)));
%     figFile = sprintf('fig/roommates_%d.eps', years(m));
%     print('-depsc', figFile);
%     close;
% end

%% Save
save('-v7.3', 'cp/roommatesdata/sqdiff_bias.mat', 'features', 'Ys');
