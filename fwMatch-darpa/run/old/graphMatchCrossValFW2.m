lambdas  = [0.0001 0.0003 0.001 0.003 .01 .03 .1 .3 1 3 10 30 100];
rhos     = 1.0;
algo     = 'fw2';
nLambdas = length(lambdas);
nRhos    = length(rhos);

tuneParams = {'MaxIter', 10000, 'lineSearch', true};

startup;
matlabpool open 8;

vll = zeros(nLambdas, nRhos);
thetas = cell(nLambdas, nRhos);

parfor n = 1:nLambdas
%for n = 1:nLambdas
    for m = 1:nRhos
        lambda = lambdas(n);
        rho    = rhos(m);        
        [vll(n,m), thetas{n,m}] = runGraphMatchMinReg(baseName, gap, lambda, rho, algo, tuneParams);
    end
end

