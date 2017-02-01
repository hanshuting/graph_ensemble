%% Load data
data = load('cp/roommatesdata/data.mat');


%% Load results; Construct score and label vectors
res  = load('expt/roommateExplore/result3.mat');
theta = res.paramsHist{end};

mainFeatureIxs = 1:14;
thetaBaseline = zeros(size(theta));
% This was the FAIREST benchmark.
thetaBaseline(mainFeatureIxs) = -1;
% thetaBaseline(1:end) = -1;

K = size(data.features, 2);
m = 3;
A = zeros(size(data.features{m,1}));
Abaseline = zeros(size(A));
for k = 1:K
    A = A + theta(k) .* data.features{m,k};
    Abaseline = Abaseline + thetaBaseline(k) .* data.features{m,k};
end

yVec = vecUT(data.Ys{m}, true);
AVec = vecUT(A, true);
AbaselineVec = vecUT(Abaseline, true);

%% And plot

% TODO: Run this several times for different predictors, then overlay onto
% one plot.
%
% (Well, they take different datasets in columns, so we could just squeeze
% them together and use its function to plot.)


yMat = repmat(logical(yVec), 1, 2);
AMat = [AVec AbaselineVec];
[auc,fpr,tpr] = fastAUC_m(yMat, AMat, false);

% Add the pure random case -- a straight line
fpr = [fpr fpr(:,1)];
tpr = [tpr fpr(:,1)];

% Plot it ourseves
f = paperFig(0.5, 0.25);

plot(fpr(:,1), tpr(:,1), '--', ...
     fpr(:,2), tpr(:,2), '-', ...
     fpr(:,3), tpr(:,3), '-.');

xlabel('False Positive');
ylabel('True Positive');
legend(sprintf('MLE-Struct AUC=%.3f', auc(1)), ...
       sprintf('Baseline AUC=%.3f', auc(2)), ...
       'Random AUC=0.500', ...
       'Location', 'best');
auc

print -dpdf fig/roommate_roc.pdf

%% other one
