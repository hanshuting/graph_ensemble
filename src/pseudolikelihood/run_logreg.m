addpath(genpath('~/bin/liblinear-1.91/matlab/'))
addpath('../');
addpath('../../thirdparty/blossom/');
setupData = false;
if(setupData)
    %load('../../cp/roommatesdata/data.mat');
    importanceWeight = 5;
    negSampleFraction = 0.1;
    
    %make the train set
    
    XtrainFull = [ConvertData(features(1,:));ConvertData(features(2,:))];
    Xtest = ConvertData(features(3,:));
    
    labels = [Ys{1}(:)' Ys{2}(:)'];
    observedInds = find(labels == 1);
    unobservedInds = find(labels == 0);
    
    negativeTrainIndicesToUseFull = unobservedInds(randperm(length(unobservedInds)));
    Ytest = Ys{3}(:);
    %now, subsample
    for negSampleFraction = [.001 .01 .1]
        numTake = floor(negSampleFraction*length(negativeTrainIndicesToUseFull));
        negativeTrainIndicesToUse = negativeTrainIndicesToUseFull(1:numTake);
        
        trainIndicesToUse = [observedInds negativeTrainIndicesToUse];
        Xtrain = XtrainFull(trainIndicesToUse,:);
        Ytrain = labels(trainIndicesToUse);
        trainWeights = weightsFull(trainIndicesToUse);
        save(['../../cp/roommatesdata/data-' num2str(10000*negSampleFraction)],'Xtrain','Ytrain','Xtest','Ytest');
    end
else
    negSampleFraction = 0.1;
    file = ['../../cp/roommatesdata/data-' num2str(10000*negSampleFraction) '.mat'];
    fprintf('loading from %s\n',file);
    load(file);
end
%todo: find indices that have observed arcs
%make the test set
disp('training');
%liblinear expects -1/+1 labels
Ytrain = 2*Ytrain - 1;
Ytest = 2*Ytest-1;

% Ytest = Ytrain';
% Xtest = Xtrain;

allPrecs = [];
allRecs = [];
useLocalClassification = true;
accuracies = [];
for weight = [500]
    w = train(Ytrain',Xtrain,['-s 6 -B 1 -w1 ' num2str(weight)]);
    
    %now, predict. can't do it all in one go because memory requirements are too big.
    chunkSize = 1000000;
    len = length(Ytest);
    start = 1;
    stop = chunkSize;
    
    predEnergy = [];
    while start < len
        i = start:stop;
        [pred, acc, probs] = predict(Ytest(i),Xtest(i,:),w);
        predEnergy = [predEnergy; probs];
        start = min(len,start + chunkSize);
        stop = min(len,stop + chunkSize);
    end
    
    if(useLocalClassification)
        
        %%%%Get the predictions using local classification
        predProbs = 1./(1 + exp(-predEnergy));
        threshRange = linspace(0,1,100);
        threshRange = threshRange(2:(end-1));
        testLabels = Ytest == 1;
        numTest = sum(testLabels);
        prec = [];
        rec = [];
        for thresh = threshRange
            pred = (predProbs > thresh);
            numPred = sum(pred);
            recoveredLinks = pred & testLabels;
            numRecovered = sum(recoveredLinks);
            prec = [prec numRecovered/numPred];
            rec = [rec numRecovered/numTest];
            
        end
        plot(prec,rec);
        
        allPrecs = [allPrecs prec'];
        allRecs = [allRecs rec'];
        n = sqrt(length(predEnergy));
        probMatrix = reshape(predProbs,n,n);
        labelMatrix = reshape(Ytest,n,n);
        idxs = [];
        for i = 1:n
            row = probMatrix(i,:);
            trueRow = labelMatrix(i,:);
            truth = find(trueRow == 1);
            [s ii] = sort(row,'descend');
            idx = ii(truth);
            idxs = [idxs idx];
        end
        
        %%%%%%%%%%%%%%%%%%%
    else
        %%%%%%%%%%%%%%%%%%%
        %Alternatively, get predictions using perfect matching
        n = sqrt(length(predEnergy));
        predProbs = 1./(1 + exp(-predEnergy));
        indsToIgnore = find(predProbs < 0.0001);
%         shift = 2*abs(min(predEnergy));
%         predEnergyShift = predEnergy + shift;
%         predEnergyShift(indsToIgnore) = 0;
%        scoreMatrix = reshape(predEnergyShift,n,n);
        predProbsSparse = predProbs;
        predProbsSparse(predProbs < 0.25) = 0;
        scoreMatrix = reshape(predProbsSparse,n,n);
        tic;
        predMatrix = perfectMatching(-scoreMatrix);
        el = toc;
        fprintf('%f to compute max matching\n',el);
        preds = predMatrix(:);
        testLabels = Ytest == 1;
        numRecoveredLinks = sum(preds & testLabels);
        accuracy = numRecoveredLinks/n;
        accuracies = [accuracies accuracy];
        %%%%%%%%%%%%%%%%%%%
    end
    
    
end
