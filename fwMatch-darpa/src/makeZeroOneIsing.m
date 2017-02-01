function [YN, YE, labels] = makeZeroOneIsing(Mones, Mzeros, Mrandoms, N, edges)
%makeZeroOneIsing Return (attractive) Ising model samples 

M = Mones + Mzeros + Mrandoms;

Nedges = size(edges, 1);

YNones    = zeros(N * Mones, 2);
YNones(:,2) = 1;

% Edge corresponding to ones (really label 2) is (2, 2) -> linear ind 4.
YEones      = zeros(Nedges * Mones, 4);
YEones(:,4) = 1;

YNzeros   = zeros(N * Mzeros, 2);
YNzeros(:,1) = 1;

% Edge corresponding to zeros (really label 1) is (1, 1) -> linear ind 1.
YEzeros      = zeros(Nedges * Mzeros, 4);
YEzeros(:,1) = 1;

randomLabels = randi(2, N * Mrandoms, 1);
YNrandom   = zeros(N * Mrandoms, 2);

for i = 1:(N * Mrandoms)
    YNrandom(i,randomLabels(i)) = 1;    
end

YErandom   = zeros(Nedges * Mrandoms, 4);
for m = 1:Mrandoms
    nodeBase = (m - 1)*N;
    edgeBase = (m - 1)*Nedges; 
    
    for e = 1:Nedges
        i   = edges(e,1);
        j   = edges(e,2);        
        
        vi  = randomLabels(nodeBase + i);
        vj  = randomLabels(nodeBase + j);
        
        col = sub2ind([2, 2], vi, vj);

        row = edgeBase + e;
        YErandom(row,col) = 1;
    end
end

YN = vertcat(YNones, YNzeros, YNrandom);
YE = vertcat(YEones, YEzeros, YErandom);

labels = cell(M, 1);

for m = 1:Mones
    labels{m} = 2*ones(N, 1);
end
for m = 1:Mzeros
    labels{Mones + m} = ones(N, 1);
end
for m = 1:Mrandoms
    nBegin = (m-1)*N + 1;
    nEnd   = m*N;    
    labels{Mones + Mzeros + m} = randomLabels(nBegin:nEnd);
end
