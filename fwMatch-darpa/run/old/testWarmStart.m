%% blah blah

pairFile = 'houses10_train.txt.simple';
Fkm = parsePair(pairFile);
[D, ~, K, M] = size(Fkm);

Msub = 5;
Fkm = Fkm(:,:,:,1:Msub);

% In the houses data, the ground truth is the identity permutation:
% "Each frame in this sequence has been hand-labeled, with the same 30
%  landmarks identified in each frame"
Ys = repmat(1:D, Msub, 1);

nTrials = 5;


warmTimes(nTrials) = 0;
coldTimes(nTrials) = 0;
for t = 1:nTrials
    %[thetaExact, histExact]   = learnFeats(Ys, Fkm, @exactLikeFeats, 'TolFun', 1e-2, 'initStep', 5);
    tic;
    [thetaApprox, histApprox] = learnFeats(Ys, Fkm, @betheLikeFeats, 'TolFun', 1e-2, 'initStep', 1, 'warmStart', true);
    warmTimes(t) = toc;
    fprintf('Warm start took %g seconds\n', warmTimes(t));

    tic;
    [thetaApprox, histApprox] = learnFeats(Ys, Fkm, @betheLikeFeats, 'TolFun', 1e-2, 'initStep', 1, 'warmStart', false);
    coldTimes(t) = toc;
    fprintf('Cold start took %g seconds\n', coldTimes(t));
end

