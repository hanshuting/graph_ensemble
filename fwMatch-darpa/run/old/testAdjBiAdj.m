T = 100;
for t = 1:T
    M = randi([1, 20]);
    N = randi([1, 20]);
    B = rand(M, N);
    assertElementsAlmostEqual(adjacencyToBiadjacency(biadjacencyToAdjacency(B)), B);
end
