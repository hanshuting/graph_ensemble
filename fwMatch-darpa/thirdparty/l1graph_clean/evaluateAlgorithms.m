
%INPUT: 
%	lcnst: penalty lambda = lcnst x asymptotic-lambda 
%	tolrcnst: threshold TOLR = tolrcnst x lambda
%OUTPUT: 
%	probsucc: PROB. OF SUCCESS
%	fracdisagree:  FRACTION OF DISAGREEING EDGES

function [probsucc,fracdisagree] = evaluateAlgorithms(graphs,Xs,trialSize,lcnst,tolrcnst,C,alpha)

lambdaAsymp = lcnst * sqrt(log(numNodes)/numSamples);
lambda = lambdaAsymp;
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

numExperiments = max(size(graphs));

for(numE = 1:numExperiments)
  graph = graphs{numE};
  XX = Xs{numE};
  numNodes=size(XX,2);
  numTrials = size(XX,1)/trialSize;
  
  for(numT = 1:numTrials)
    numT
    X = XX((numT-1)*trialSize+1:numT*trialSize,:);
	
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
    % TONY, add your code here
    [S,lambdas,iters]=multigraphical(X,C,alpha);
    graphret{3} = S>thresh;
    
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
end

tmpos
tmpneg
probsucc = probsucc/numTrials
fracdisagree = fracdisagree/numTrials
filestr = sprintf('outputx/ss-%d-%d-%.2f.mat',numNodes,numSamples,lcnst);
save(filestr,'probsucc','fracdisagree','betamat');

