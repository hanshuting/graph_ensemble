%% Comments by Henrique Gubert
%
% Generate 100 graphs and try to recover graph structure using lasso
% regression. After lasso regression, tries 5 different methods for
% combining the relevant features of each regression.
%
% Parameters
%   numNodes: number of nodes of the graphs generated 
%	numSamples: number of samples in each experiment
%	lcnst: constant that defines lambda (regularizing coefficient), 
%       lambda = lcnst * sqrt(log(numNodes)/numSamples)
%	tolrcnst: constant that defines the tolerance threshold. This threshold
%       is how high a coefficient of a node must be to be selected as 
%       'relevant' to the node being predicted in the lasso. 
%       TOLR = tolrcnst x lambda
%
% Return Values
%   probsucc: fraction of the 100 trials and 5 methods (500 settings)
%       that managed to recover the graph structure perfectly
%   fracdisagree: average edge disagreement per trial

%%
%INPUT: 
%	numNodes: NUMBER OF VARIABLES 
%	numSamples: NUMBER OF SAMPLES 
%	lcnst: penalty lambda = lcnst x asymptotic-lambda 
%	tolrcnst: threshold TOLR = tolrcnst x lambda
%OUTPUT: 
%	probsucc: PROB. OF SUCCESS
%	fracdisagree:  FRACTION OF DISAGREEING EDGES

function [probsucc,fracdisagree] = structSearchInstance(numNodes,numSamples,lcnst,tolrcnst)

topology = 'star'
%degmax = ceil(log(numNodes));
%degmax = numNodes - 1;
%degmax = 5;
%dcoup = 5;

degmax = ceil(0.1 * numNodes)
dcoup = min(0.1,2.5/degmax)
edgeWeightSign = 'mixed'

%TONY
topology = 'complete';
degmax=ceil(0.3*numNodes);
dcoup = min(0.1,2.5/degmax);
edgeWeightSign = 'mixed'


mtype = java.util.Hashtable;
mtype.put('numNodes',numNodes);
mtype.put('degmax',degmax);
mtype.put('topology',topology);
mtype.put('edgeWeightSign',edgeWeightSign);
mtype.put('coupStrength',dcoup);

nexp = 1
lambdaAsymp = lcnst * sqrt(log(numNodes)/numSamples)

lambda = lambdaAsymp

%[graph,X] = getSearchInstance(mtype,numSamples,lambda);

graph = getGraphModel(mtype);
truegraph = (graph ~= 0);	

numTrials = 100;
numMethods = 5;

probsucc = zeros(numMethods,1);
fracdisagree = zeros(numMethods,1);
graphret = cell(numMethods,1);
for(m = 1:numMethods)
	graphret{m} = zeros(numNodes);
end
tmp = zeros(numMethods,1);
tmpos = zeros(numMethods,1);
tmpneg = zeros(numMethods,1);
betamat = cell(numTrials,1);


figure(4); imagesc(graph)
graph
sum(truegraph)


for(numT = 1:numTrials)
	numT
	
	X = getSamples(graph,numSamples,topology);
	
	disp('L1')
	graphL1 = zeros(numNodes);
	for s = 1:numNodes
		s
		Vbs = setdiff(1:numNodes,[s]);
		[thetanbr1,betamat{numT}] = logisticRegBoyd(X(:,s),X(:,Vbs),lambda,lcnst,tolrcnst);
		graphL1(Vbs,s) = thetanbr1;
	end	
	%Combining neighborhood estimates via AND, OR
	graphret{1} = ((graphL1 + graphL1') == 2);
	graphret{2} = ((graphL1 + graphL1') > 0);
	
    % X is a matrix of samples from the graph of size
    % numSamples times numNodes
    figure(1); imagesc(X)
    figure(2); imagesc(   ((graphL1+graphL1')==2) );
    figure(3); imagesc(   ((graphL1+graphL1')>0) );
    pause
        
	disp('PC')
	%graphret{3} = pcMatWrapper(X); % TONY, add your code here
	disp('CL')
	[graphret{4},graphret{5}] = chowliu(X,degmax);
		
	
	for(m = 1:numMethods)
		tmpos(m) = sum(sum(graphret{m} > truegraph));
		tmpneg(m) = sum(sum(graphret{m} < truegraph));
		tmp(m) = sum(sum(graphret{m} ~= truegraph));
		probsucc(m) = probsucc(m) + (tmp(m) == 0);
		fracdisagree(m) = fracdisagree(m) + tmp(m);
	end
end

tmpos
tmpneg
probsucc = probsucc/numTrials
fracdisagree = fracdisagree/numTrials
filestr = sprintf('outputx/ss-%d-%d-%.2f.mat',numNodes,numSamples,lcnst);
save(filestr,'probsucc','fracdisagree','betamat');

