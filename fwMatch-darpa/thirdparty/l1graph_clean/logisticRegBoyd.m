%% Comments by Henrique Gubert
% Trains logistic regressor with L1 regularizer. Identifies which features
% are relevant to the classifier according to some tolerance threshold.
%
% Parameters
% y,X: labels and features, where each row is one observation 
% lambda: regularization coefficient
% lcnst: UNUSED
% tolrcnst: tolerance=tolrcnst*lambda. Tolerance is how high a coefficient
%   must be for a feature to be selected as 'relevant'
%
% Return Values
% nbrlasso: indicator row vector of the relevant features
% beta: row vector with coefficients of the regressor [b0 b1 b2 ...]
%
%%
%INPUT: RESPONSE Y, PREDICTOR X, PENALTY LAMBDA
%OUTPUT: NEIGHBORHOOD NBRLASSO, COEFFICIENTS BETA 
function [nbrlasso,beta] = logisticRegBoyd(y,X,lambda,lcnst,tolrcnst)

n = size(X,1);
p = size(X,2);
%TOLR = 0.02;
TOLR = tolrcnst * lambda;

Xfile = tempname;
yfile = tempname;
modelfile = tempname;

mmwrite(Xfile,X);
mmwrite(yfile,y);

system('thirdparty/l1_logreg-0.8.2/c_src/l1_logreg_train -q -s %s %s %f %s',Xfile,yfile,lambda,modelfile);
model_lg = full(mmread(modelfile));
beta = model_lg(2:length(model_lg));

nbrlasso = abs(beta) >= TOLR;
%numnbr = sum(nbrlasso)
%minnbr = min(beta(abs(beta) >= TOLR))

