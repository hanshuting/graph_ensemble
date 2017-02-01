%% Sample something
T = 10;
D = 5;
scale = 10;


minBetheTimes = zeros(T, 1);
bpTimes       = zeros(T, 1);
normBpTimes   = zeros(T, 1);

% Biadjacency matrix
for t = 1:T
    %% One iter
    B = scale * rand(D, D);

    tic;
    [muMin, FBethe, FBetheHist] = minBethe([], B);
    minBetheTimes(t) = toc;
       
    A = biadjacencyToAdjacency(B);

    tic;
    [logZbp, bpMuMin, bpMsgs] = bpPerfectMatching(A);
    biBpMuMin = adjacencyToBiadjacency(bpMuMin);    
    bpTimes(t) = toc;
    
    tic;
    Amax  = max(abs(A(:)));
    Anorm = A ./ Amax;    
    [logZNormBp, normBpMuMin, normBpMsgs] = bpPerfectMatching(Anorm);
    biNormBpMuMin = adjacencyToBiadjacency(normBpMuMin)    

    normBpTimes(t) = toc;
    
    logZNormBpRecovered = logZNormBp + D*Amax*log(2);
    
    figure; imagesc(abs(muMin - biBpMuMin));   colorbar; title('|muMin - biBpMuMin|');
    figure; imagesc(abs(muMin - biNormBpMuMin)); colorbar; title('|muMin - biNormBpMuMin|');

    fprintf('FBethe = %g (%g s); -logZbp = %g (%g s); -logZNormBpRecovered = %g (%g s);', ...
        FBethe, minBetheTimes(t), ...
        -logZbp, bpTimes(t), ...
        -logZNormBpRecovered, normBpTimes(t));
    
    fprintf('bp was %g times faster\n', minBetheTimes(t) / bpTimes(t));
end

%% Plots
figure; hist(minBetheTimes); title('Frank Wolfe Runtimes');
figure; hist(bpTimes);       title('BP Runtimes');

