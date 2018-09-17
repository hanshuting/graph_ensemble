function [M, mu] = clusterConsensusLouvainAdjmat(CIJ,gamma)
% Consensus clustering with Louvain method for community detection
% INPUT:
%     data: raster matrix of spikes (N_frames by N_neurons)
%     thr: threshold of connection weight (0-1, percentage)
%     gamma: resolution parameter
%     mu: mixing parameter
% OUTPUT:
%     M: vector of community identieis
%     CIJ: undirected weighted connectivity matrix
%     R: pairwise cell correlation

% perform Newman's spectral community detection multiple times
nrep = 100;
Q = zeros(nrep,1);
comm = zeros(length(CIJ),nrep);
for ii = 1:nrep
    [comm(:,ii),Q(ii)] = modularity_und(CIJ,gamma);
end

% Louvain method
Q = zeros(1000,1);
CL = zeros(length(CIJ),1000);
for i = 1:1000
    [CL(:,i),Q(i)] = community_louvain(CIJ,gamma);
end

% agreement matrix from clusters
D = agreement(CL,1000);  % agreement (consensus) matrix from clusters;

% Consensus clustering
CC = consensus_und(D,0.8,200);

% Modify the CC result by removing cells with small
% within-community degrees and groups with few cells;
Z = module_degree(CIJ,CC,0); % calculate the module degree of individual cells;
CC(Z<3) = nan; % each cell should connect with >= 3 cells;
mC = nanmax(CC);
a = histcounts(CC,mC);
b = find(a<4); % need to determine the threshold of group size;   
L = ismember(CC,b);
CC(L) = nan;    

% calculate the mixing parameter;
mu = mixingcoef(CIJ,CC);
CC(isnan(CC)) = 0;

% 
v = unique(nonzeros(CC));
[~,M] = ismember(CC,v);

end
     

    
    