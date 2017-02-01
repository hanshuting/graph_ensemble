load ../../../cp/results/bmrm_crossval

for nb = 1:nBaseNames
    testLossLog = sprintf('%s_test_loss.txt', baseNames{nb});
    testLearnLoss = importdata(testLossLog);
    CV{nb}.testLearnLoss = testLearnLoss;

    % Reorder
    CV{nb} = CV{nb}(:, {'gaps', 'testLearnLoss', 'validLearnLoss', 'trainLearnLoss', 'learnTime', 'nIters', 'lambda'})
end

save ../../../cp/results/bmrm_crossval_full -v7.3
