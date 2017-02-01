baseNames = {'hotel' 'house'}
nBaseNames = length(baseNames);

gaps = [0 10 20 30 40 50 60 70 80 90]';
nGaps = length(gaps);
lambdas = [0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100]';
nLambdas = length(lambdas);

CV = cell(size(nBaseNames));

%outFileTemplate = '${baseName}_gap${gap}_lambda${lambda}_linear.config'
outFileTemplate = '%s_gap%d_lambda%g_linear.lossLog'

for nb = 1:nBaseNames
    validLearnLoss = zeros(nGaps, 1);
    trainLearnLoss = zeros(nGaps, 1);
    learnTime      = zeros(nGaps, 1);
    nIters         = zeros(nGaps, 1);
    lambda         = zeros(nGaps, 1);

    for ng = 1:nGaps
        nts = zeros(nLambdas, 1);
        lts = zeros(nLambdas, 1);
        vll = zeros(nLambdas, 1);
        tll = zeros(nLambdas, 1);
        for nl = 1:nLambdas
            lossLog = sprintf(outFileTemplate, baseNames{nb}, gaps(ng), lambdas(nl));
            row = importdata(lossLog);
            nts(nl) = row(1);
            lts(nl) = row(2);
            tll(nl) = row(3);
            vll(nl) = row(4);
        end
        [validLearnLoss(ng), i] = min(vll);
        trainLearnLoss(ng) = tll(i);
        lambda(ng) = lambdas(i);
        nIters(ng) = nts(i);
        learnTime(ng) = lts(i);
    end

    CV{nb} = dataset(gaps, validLearnLoss, trainLearnLoss, learnTime, nIters, lambda);
end

save ../../../cp/results/bmrm_crossval -v7.3 

for nb = 1:nBaseNames
    saveFile = sprintf('%s_lambdas.txt', baseNames{nb});
    fid = fopen(saveFile, 'w');
    for ng = 1:nGaps
        fprintf(fid, '%g\n', CV{nb}.lambda(ng));
    end
    fclose(fid);
end

