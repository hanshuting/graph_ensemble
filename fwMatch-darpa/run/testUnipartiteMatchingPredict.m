%% Load data
data = load('cp/roommatesdata/data.mat');


%% Load results; Construct score and label vectors
res  = load('expt/roommateExplore/result3.mat');

thetaEnd = res.paramsHist{end};

% The kth entry of this vector is the return value of
% thetaWavg(k).
K = size(data.features, 2);
T = length(res.paramsHist);
thetaWavgAtK = zeros(K, T);

for k = 1:T
    for t = 1:k
        thetaWavgAtK(:,t) = thetaWavgAtK(:,t) + t*res.paramsHist{t};
    end
    
    thetaWavgAtK(:,k) = 2/(k*(k+1)) * thetaWavgAtK(:,k);
end

thetaWavgEnd = thetaWavgAtK(:,end);

% third year...
M = 3;
Y = data.Ys{3};
yVec = logical(vecUT(Y, true));

maxIter = 20;

%% Run just the end
% Actually, we just needed final iterate and average.
objEnd  = UnipartiteMatchingPredict(thetaEnd, data.features(M,:), 'YTrue', Y);
probEnd = BatchFW(objEnd, 'MaxIter', maxIter, 'printInterval', 1);
probEnd.run();
[aucEnd,  fprEnd,  tprEnd]  = fastAUC_m(yVec, objEnd.x);

%% Run Wavg
objWavg  = UnipartiteMatchingPredict(thetaWavgEnd, data.features(M,:), 'YTrue', Y);
probWavg = BatchFW(objWavg,'MaxIter', maxIter, 'printInterval', 1);
probWavg.run();
[aucWavg, fprWavg, tprWavg] = fastAUC_m(yVec, objWavg.x);

aucEnd
aucWavg
