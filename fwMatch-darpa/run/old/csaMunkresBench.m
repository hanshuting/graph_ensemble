%% Blah

D = 2;
T = 1;

munkresTimes = zeros(T, 1);
munkresOptTimes = zeros(T, 1);
csaTimes = zeros(T, 1);

%% Break
for t = 1:T
    %% One try
    C = rand(D);
    
    tic;
    munkresPerm = munkres(C);    
    munkresTimes(t) = toc;
    
    tic;
    munkresOptPerm = munkresOpt(C);
    munkresOptTimes(t) = toc;    
    
%     assertElementsAlmostEqual(munkresPerm, munkresOptPerm);
    
    tic;
    [spC, scale] = sparsifyAndRound(C);    
    edges = csaAssign(2*D, spC);
    csaTimes(t) = toc;
    
    csaPerm = edges(2,:) - D;    
            
    munkresCost = matchingCost(C, munkresPerm);
    csaCost     = matchingCost(C, csaPerm);
       
    assertElementsAlmostEqual(munkresCost, csaCost);
    
%     assertElementsAlmostEqual(munkresPerm, edgesPerm(2,:));    
end

figure; hist(munkresTimes); title('Munkres times');
figure; hist(munkresOptTimes); title('Munkres Opt times');
figure; hist(csaTimes), title('CSA Times');

