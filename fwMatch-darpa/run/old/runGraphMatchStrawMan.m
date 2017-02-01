%% blah blah

%pairFile = 'houses0_train.txt.simple';
%Fkm = parsePair(pairFile);
%[D, ~, K, M] = size(Fkm);
%% In the houses data, the ground truth is the identity permutation:
%% "Each frame in this sequence has been hand-labeled, with the same 30
%%  landmarks identified in each frame"
%Ys = repmat(1:D, M, 1);
%
%[thetaApprox, histApprox] = learnFeats(Ys, Fkm, @betheLikeFeats, 'TolFun', 1e-2, 'initStep', 5, 'lambda', 0.1);
%

theta = zeros(60, 1);

[~, noLearnTrainLoss] = evalGraphMatch(thetaApprox, 'houses0_train.txt.simple')
[~, noLearnTestLoss]  = evalGraphMatch(thetaApprox, 'houses0_test.txt.simple')

