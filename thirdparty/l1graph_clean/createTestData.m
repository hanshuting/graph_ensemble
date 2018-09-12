function [graphs,Xs] = createTestData(numNodes,numSamples,numTrials)
%
% function [graphs,Xs] = createTestData(numNodes,numSamples,numTrials)
%INPUT: 
%	numNodes: NUMBER OF VARIABLES 
%	numSamples: NUMBER OF SAMPLES 
%	numTrials: NUMBER OF TRIALS
%OUTPUT: 
%	graphs: the graph models
%	Xs: cell array of samples

graphs=cell(numTrials);
Xs=cell(numTrials);
for j=1:numTrials
  j
  topology = 'complete';
  degmax=ceil(0.3*numNodes);
  dcoup = min(0.1,2.5/degmax);
  edgeWeightSign = 'mixed';
  mtype = java.util.Hashtable;
  mtype.put('numNodes',numNodes);
  mtype.put('degmax',degmax);
  mtype.put('topology',topology);
  mtype.put('edgeWeightSign',edgeWeightSign);
  mtype.put('coupStrength',dcoup);
  graph = getGraphModel(mtype);
  imagesc(graph)
  X = getSamples(graph,numSamples,topology);
  graphs{j}=graph;
  Xs{j}=X;
end
save data.mat graphs Xs