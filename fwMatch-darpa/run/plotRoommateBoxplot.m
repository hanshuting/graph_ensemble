%% Load data
load expt/roommateTest/config1_cp.mat;

%% Normalized err
normalizedErrs = errs ./ denoms;

% But we really want accuracy
normalizedAcc = 1 - normalizedErrs;
randomAccMean = 1 - randomErrMeans;

%% Gain over random
overRandomAccs = bsxfun(@times, normalizedAcc, 1./ randomAccMean);

figure; boxplot(overRandomAccs);

%% Not a boxplot
ratio = mean(normalizedAcc, 1) ./ randomAccMean;
figure; plot(mean(normalizedAcc, 1) ./ randomAccMean); title('Mean Acc / Mean Random Acc');

%% Not a ratio
figure; boxplot(normalizedAcc); title('Normalized accuracies');
figure; plot(randomAccMean);    title('Random accuracies');