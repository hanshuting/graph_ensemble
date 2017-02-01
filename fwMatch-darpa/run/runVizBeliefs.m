%% Load models
test_modelZ = load('cp/justinHorseCache/test_modelZ.mat');
test_modelZ = test_modelZ.modelZ;

%% Configuration stuff (look up!)
image_id = 8; % read from configuration file.
iter_per_macro = 1000;
% In horseReproduceCorrect configurations, we have
% params.M           = 1000;


%% Ground truth image
fileSuffix = sprintf('%d.png', 200 + image_id - 1);
filePath   = ['data/JustinHorses/HorsesTest/image-' fileSuffix];
image = importdata(filePath);
figure;
imshow(image);
title('Raw Image');

copyfile(filePath, ['fig/horse_margs/raw-' fileSuffix]);

%% Ground truth segmentation
truth = importdata(['data/JustinHorses/HorsesTestLabels/mask-' fileSuffix]);
figure;
yTrue = truth == 2;
imshow(yTrue);
title(sprintf('Ground truth image %d', image_id));

imwrite(yTrue, ['fig/horse_margs/truth-' fileSuffix]);

%% Load my beliefs

% The iteration selects are approximately time-equal. I used the data
% from plotHorseCurves to check.

% 250 is the correct number

for macro_iter = [2 100 250]
    chunk_id   = floor((macro_iter-1) / 50) + 1;
    bMat = load(sprintf('expt/mWavgHorseViz/result%d.mat', chunk_id));

    bPos = bMat.b_is{macro_iter}(2,:);
    bMPM = bPos > 0.5;
    mHammingErr = mean(bMPM.' ~= vec(yTrue));

    ly = test_modelZ{image_id}.ly;
    lx = test_modelZ{image_id}.lx;

    iter = iter_per_macro * macro_iter;
    figure;
    bImg = reshape(bPos, ly, lx);
    imshow(bImg);
    title(sprintf('My Algorithm iter %d, image %d, hamming %g', iter, image_id, mHammingErr));

    imwrite(bImg, sprintf('fig/horse_margs/my-belief-iter%d-%s', iter, fileSuffix));

    fprintf('MY ALGORITHM ITER %d, HAMMING ERR %g\n', iter, mHammingErr);
end

%% Load Justin beliefs
for j_iter = [1 17 42 96]
    jMat = load('expt/justinHorseCorrectViz/result1.mat');
    
    jPos = jMat.b_is{j_iter}(2,:);
    jMPM = jPos > 0.5;
    jHammingErr = mean(jMPM.' ~= vec(yTrue));
    
    ly = test_modelZ{image_id}.ly;
    lx = test_modelZ{image_id}.lx;

    figure;
    jImg = reshape(jPos, ly, lx);
    imshow(jImg);
    title(sprintf('Their Algorithm iter %d, image %d, hamming %g', j_iter, image_id, jHammingErr));
    
    imwrite(jImg, sprintf('fig/horse_margs/justin-belief-iter%d-%s', j_iter, fileSuffix));
    
    fprintf('THEIR ALGORITHM ITER %d, HAMMING ERR %g\n', j_iter, jHammingErr);
end