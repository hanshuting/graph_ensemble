function [statuses theta ]= GDBPLearning
addpath('src/BPGD');
addpath('thirdparty/minFunc_2012/');
clear;
baseName = 'house';
lambda = 100;
gap = 50;

rho = 1;
train    = load(sprintf('cp/imagedata/%s/%d_train.mat', baseName, gap));
validate = load(sprintf('cp/imagedata/%s/%d_validate.mat', baseName, gap));
test     = load(sprintf('cp/imagedata/%s/%d_test.mat', baseName, gap));
[M, K] = size(train.Fsquare);
D = size(train.Fsquare{1}, 1);
rhoVec = rho * ones(1, D);

[Mtest, Ktest] = size(test.Fsquare);
assert(Ktest == K);
Ys = repmat(1:D, Mtest, 1);
testYs = cell(Mtest, 1);
for m = 1:Mtest
    testYs{m} = Ys(m,:);
end
testSet  = BipartiteMatchingTestSet2(testYs, test.Fsquare);

% In the houses data, the ground truth is the identity permutation:
% "Each frame in this sequence has been hand-labeled, with the same 30
%  landmarks identified in each frame"
Ys = repmat(1:D, M, 1);
YsCell = cell(M, 1);
for m = 1:M
    YsCell{m} = Ys(m,:);
end
obj = BipartiteMatching(YsCell, train.Fsquare, lambda);
YsCell2 = cell(M, 1);

for m = 1:M
    YsCell2{m} = expandPerm(Ys(m,:), @zeros);
end

initTheta = obj.computeParams;
features = train.Fsquare;
FXBar = computeFeatInnerProduct(features, YsCell2);

FUN = @(theta) bpLikeFeats(theta, YsCell, features, FXBar);
REG_FUN =  @(theta) add_regularizer(theta,FUN,lambda);
iter = 1;
thetas = zeros(size(initTheta,1),1000);
objectives = zeros(1,1000);
statuses = [];

    function [o g] = testFun(fun,theta)
        thetas(:,iter) = theta;
        [o g] = fun(theta);
        objectives(iter) = o;
        testErr = testSet.MAPErr(theta);
        fprintf('iter %d: %f %f grad: %f\n',iter,o,testErr,norm(g));
        status = struct('testErr',testErr,'fval',o);
        statuses = [statuses status];
        iter = iter +1;        
    end
my_fun = @(theta) testFun(REG_FUN,theta);
 myOpts.maxIters = 150;
 myOpts.stepSize = @(t) .0005;%1/(t+100);%@(t) .1/sqrt(t);
 [x,FVAL] = gradDescent(my_fun,initTheta,myOpts);
% 
% options = optimoptions(@fminunc,'GradObj','on','display','iter-detailed','algorithm','quasi-newton','hessupdate','bfgs');
% [x,FVAL] = fminunc(my_fun,initTheta,options);
%   options2.Method = 'lbfgs';
%   options2.Display = 'excessive';
%   [x,FVAL] = minFunc(my_fun,initTheta,options2);
%   theta = x;
%   plot([statuses.fval])
end
