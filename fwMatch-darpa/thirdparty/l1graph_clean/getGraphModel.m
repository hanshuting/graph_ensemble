%% Comments by Henrique Gubert
% This function generates a graphical model according to the parameters
% passed through a java.util.Hashtable.
%  
% Parameters
% numNodes: # of nodes in the model
% edgeWeightSign: ['positive','negative','mixed'], 'mixed' means randomly mixed.
% coupStrength: absolute value of the edge weights
% topology: ['complete','star','grid',...], 'complete' means all nodes 
%   connected to all nodes, 'grid' means each node connected to 4 neighbors.
% degmax: maximum number of edges in a node. Edges are set according to
%   topology and then are randomly removed to not be over degmax. Some
%   nodes may end up with less connection or no connection at all, due to
%   the simplistic algorithm.
% 
% Output
% graph: a dense N-by-N matrix with the edge weights between nodes, or zero
% if no edge. The matrix is symmetric and its diagonal is all zeros.
%
%%
function graph = getGraphModel(mtype)

numNodes = mtype.get('numNodes');
edgeWeightSign = mtype.get('edgeWeightSign');
dcoup = mtype.get('coupStrength');
topology = mtype.get('topology');
degmax = mtype.get('degmax');

if(sum(size(degmax)) == 0) degmax = (numNodes - 1); end
p = numNodes;
d = degmax;

switch edgeWeightSign

case 'positive'
	W = dcoup * ones(p);
case 'negative'
	W = -1 * dcoup * ones(p);
case 'mixed'
	W = dcoup * (2 * (rand(p) <= 0.5) - 1);
	for s = 1:p
		for t = (s+1):p
			W(t,s) = W(s,t);
		end
	end
end

%REMOVE DIAG TERMS
W = W - diag(diag(W));

%TOPOLOGY
switch topology

case 'complete'
	%Complete graph K_{n}
	%All Edges are present
	%whittle down according to sparsity -- construct an s-regular graph
	for s = 1:numNodes
		Vbs = setdiff(1:numNodes,[s]);
		nbrs = find(W(s,:));
		if(length(nbrs) > degmax)
			rednbrs = nbrs(randperm(length(nbrs)));
			rednbrs = rednbrs(1:degmax);
			nonnbrs = setdiff(nbrs,rednbrs);
			W(s,nonnbrs) = zeros(1,length(nonnbrs));
			W(nonnbrs,s) = zeros(length(nonnbrs),1);
		end
	end

	
case 'grid'
	%Grid Graph, 2D mesh

	gridSize = sqrt(numNodes);
	for i = 1:numNodes
		for j = (i+1):numNodes
			if(noGridEdge(i,j,gridSize,gridSize))
				W(i,j) = 0;
				W(j,i) = 0;
			end
		end
	end

case '8ngrid'
		%Grid Graph, 2D mesh, 8 nn

		gridSize = sqrt(numNodes);
		for i = 1:numNodes
			for j = (i+1):numNodes
				if(no8nGridEdge(i,j,gridSize,gridSize))
					W(i,j) = 0;
					W(j,i) = 0;
				end
			end
		end

case 'star'
 	graph = zeros(p);
	sig = eye(p);
	
	nbrs = 1 + randperm(p-1);
	nbrs = nbrs(1:d);
	graph(1,nbrs) = 1;
	graph(nbrs,1) = 1;
	rho = min(0.1,2.5/d);
	
	for(s = 1:p)
		if(graph(1,s) ~= 0)
			sig(1,s) = rho;
			sig(s,1) = rho;
			
			for(t = 1:p)
				if((t ~= s) && (graph(1,t) ~= 0))
					sig(s,t) = rho^2;
					sig(t,s) = rho^2;
				end
			end
		end
	end
	W = inv(sig) .* graph;


	%STAR GRAPH
	%nbrs = 1 + randperm(p-1);
	%nbrs = nbrs(1:d);
	%Wb = zeros(p);
	%Wb(1,nbrs) = 1;
	%Wb(nbrs,1) = 1;
	%W = W .* Wb;
	
case 'tree1'
	%a specific tree or a random tree
	gridSize = sqrt(numNodes);
	graph = zeros(numNodes,numNodes);
	for i = 1:gridSize
		for j = 2:gridSize
			cnode = gridSize * (i-1) + j;
			pnode = cnode - 1;		
			graph(cnode,pnode) = 1;
			graph(pnode,cnode) = 1;
		end
	end
	for i = 2:gridSize
		pnode = gridSize * (i-1);
		cnode = gridSize * i;
		graph(cnode,pnode) = 1;
		graph(pnode,cnode) = 1;
	end
	W = W .* graph;

case 'chain'
	Wb = zeros(numNodes,numNodes);
	for i = 1:(numNodes-1)
		Wb(i,i+1) = 1;
		Wb(i+1,i) = 1;
	end
	W = W .* Wb;

end

graph = W;

function nedge = no8nGridEdge(i,j,m,n)
	edge = false;
	if(mod(i,m) == 1)
		edge = (j == (i+1)) || (j == (i+m)) || (j == (i-m)) || (j == (i - m + 1)) || (j == (i +m + 1));
	else if(mod(i,m) == 0)
	 		edge = (j == (i-1)) || (j == (i+m)) || (j == (i-m)) || (j == (i - m - 1)) || (j == (i + m -1));
		else
			edge = (j == (i-1)) || (j == (i+1)) || (j == (i+m)) || (j == (i-m)) || (j == (i - m - 1)) || (j == (i + m -1)) || (j == (i - m + 1)) || (j == (i +m + 1));
		end
	end
	nedge = ~edge;

function nedge = noGridEdge(i,j,m,n)
	nedge = true;
	if(mod(i,m) == 1)
		nedge = ~((j == (i+1)) || (j == (i+m)) || (j == (i-m)));
	else if(mod(i,m) == 0)
	 		nedge = ~((j == (i-1)) || (j == (i+m)) || (j == (i-m)));
		else
			nedge = ~((j == (i+1)) || (j == (i-1)) || (j == (i+m)) || (j == (i-m)));
		end
	end
		

function [gph,rn] = intGraphTraversal(g,li,ui)
	
	gph = g;
	if(ui < li)
		rn = [];
		return;
	end

	if(ui == li)
		rn = li;
		return;
	end

	rn = li + ceil((ui - li)/2);
	[gph,rnl] = intGraphTraversal(gph,li,rn-1);
	[gph,rnr] = intGraphTraversal(gph,rn+1,ui);
	
	gph(rn,rnl) = 1;
	gph(rnl,rn) = 1;

	gph(rn,rnr) = 1;
	gph(rnr,rn) = 1;

