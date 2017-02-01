T = 10;
D = 20;
scale = 1e3;

biTimes  = zeros(T, 1);
mexBiTimes = zeros(T, 1);
uniTimes = zeros(T, 1);

for t = 1:T
    %% Blah blah    
    B = scale * randn(D,D);
    biBeliefs = 1/D * ones(D,D);
    TolFun = 1e-4;
    
    tic;
    [beliefsBi, objectivesBi] = fw_permanent_simple(biBeliefs, B, TolFun);
    biTimes(t) = toc;

    tic;
    [beliefsMexBi, objectivesMexBi] = fwBipartite_mex(biBeliefs, B, TolFun);
    mexBiTimes(t) = toc;
    
    tic;
    A = biadjacencyToAdjacency(B);
    uniBeliefs = biadjacencyToAdjacency(biBeliefs);
    [beliefsUni, objectivesUni] = fwUnipartite(uniBeliefs, A, TolFun);    
    uniTimes(t) = toc;

    recoveredBeliefs = adjacencyToBiadjacency(beliefsUni);
    figure; imagesc(abs(beliefsBi - recoveredBeliefs)); colorbar; title('fw vs uni');
    
    figure; imagesc(abs(beliefsBi - beliefsMexBi)); colorbar; title('fw vs mexFw');

    fprintf('Bipartite obj=%g (%d iters, %g s), Mex Bipartite obj=%g (%d iters, %g s), 0.5*Unipartite obj=%g (%d iters, %g s)\n', ...
        objectivesBi(end), length(objectivesBi), biTimes(t), ...
        objectivesMexBi(end), length(objectivesMexBi), mexBiTimes(t), ...
        .5*objectivesUni(end), length(objectivesUni), uniTimes(t));
end
