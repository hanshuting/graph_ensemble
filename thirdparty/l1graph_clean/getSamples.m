%% Comments by Henrique Gubert
% Get samples from a MRF. For graphs other than 'star' and 'chain' this
% function uses Gibbs sampling. For 'star' and 'chain' it uses belief
% propagation.
%
% Parameters
% graph: the edge weights matrix
% n:     number of samples
% gtype: graph topology ['complete', 'star', 'grid',...]
%
%%
function X = getSamples(graph,n,gtype)

p = size(graph,1);
X = zeros(n,p);
gridSize = sqrt(p);

%X = gibbs_mex(graph,p,n,600,400);

switch gtype

case {'star','chain'}

%BELIEF PROPAGATION
%[margn1,marge1] = belief_prop(graph);

[margn,marge] = starmarg(graph);

sorder = 1:p;
for l = 1:n
	for j = 1:p
		snode = sorder(j);
		nbrs = find(graph(snode,:));
		pnode = nbrs(find(sorder(nbrs) < sorder(snode)));
	
		if(length(pnode) > 1)
			disp('right elimination order already\n') 
		end
		if(length(pnode) == 0)
			
			X(l,snode) = (rand > margn(snode,1));
		else
			
			%sprob = marge(snode,pnode,1,X(l,pnode)+1)/margn(pnode,X(l,pnode)+1);
			sprob = marge(pnode,snode,X(l,pnode)+1,1)/margn(pnode,X(l,pnode)+1);
			
			X(l,snode) = (rand > sprob);
		end
	end
end
X = 2 * X - 1;
	
case {'grid','8ngrid','complete'}
	burnin = 200;
	sampint = 20;
	X = gibbs_mex(graph,p,n,burnin,sampint);

end %end of switch 

end
