function G = makebmatchexpanded(W,b);

%makes a bipartite graphical model for B-matching

N = size(W,1);

vals = nchoosek(1:N,b);

A = ones(2*N);

A(1:N,1:N) = zeros(N);
A(N+1:2*N, N+1:2*N) = zeros(N);

G = makegraph(A, size(vals,1));

G.bcount = b;

G.W = cat(2, cat(1,zeros(size(W)),W'), cat(1,W, zeros(size(W))));

G.vals = vals;

%set up phi potentials

for i=1:N
    for j=1:G.k
        G.phi(i,j) = sqrt(W(i,G.vals(j)));
        G.phi(i+N,j) = sqrt(W(G.vals(j),i));
    end
end

%set up pairwise clique functions

for i=1:N
    for j=1:N
        table = ones(G.k);%logical(ones(G.k));
        for r=1:size(vals,1)
            for s=1:size(vals,1)
                if (xor(j==r, i==s))
                    table(r,s) = 0;
                end
            end
        end
                
        G.psi{i}{N+j}=table;
        G.psi{N+i}{j}=table;
        
        G.b2{i}{N+j}=ones(G.k);
        G.b2{N+i}{j}=ones(G.k);
        
        %initialize messages to random
        G.m{i}{N+j} = rand(G.k,1);
        G.m{N+j}{i} = rand(G.k,1);
    end
end

%initialize one belief to delta by setting all its corresponding messages
%to delta
%for i=1:N
%    G.m{i}{N+1} = zeros(G.k,1);
%    G.m{i}{N+1}(1) = 1;
%end

for i=1:length(G.R)
    source = G.R(i);
    dest = G.C(i);
    G.phi2{source}{dest} = G.psi{source}{dest}.*repmat(G.phi(source,:)',1,G.k).*repmat(G.phi(dest,:),G.k,1);
end

G.iter = 1;

G.counttree = zeros(2*N,2*N);