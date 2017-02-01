lambdas  = [0.0001 0.0003 0.001 0.003 .01 .03 .1 .3 1 3 10 30 100];
rhos     = 0.5:0.05:1.0;
algo     = 'fmincon';
nLambdas = length(lambdas);
nRhos    = length(rhos);

tuneParams = {'MaxIter', 1000, 'MaxFunEvals', 50000};

startup;
matlabpool open 8;

vll = zeros(nLambdas, nRhos);
thetas = cell(nLambdas, nRhos);

parfor n = 1:nLambdas
%for n = 1:nLambdas
    for m = 1:nRhos
        lambda = lambdas(n);
        rho    = rhos(m);        
        [vll(n,m), thetas{n,m}] = runGraphMatchMinReg(gap, lambda, rho, algo, tuneParams);
    end
end

[bestVll, iBestVll]  = min(vll);
[nBestVll, mBestVll] = ind2sub([nLambdas, nRhos], iBestVll);

bestTheta = thetas{iBestVll};
bestLam = lambdas(nBestVll);
bestRho = rhos(mBestVll);

test = load(sprintf('cp/imagedata/%d_test.mat', gap));
[testLearnLoss, testNoLearnLoss] = evalGraphMatch(bestTheta, test);
clear test;

saveFile = sprintf('cp/results/crossval_gap_%d.mat', gap);
save('-v7.3', saveFile);
