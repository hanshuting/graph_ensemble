function G = makegraph(A,k);

%G = makegraph(A,k);
%A is the NxN adjacency matrix
%k is the number of values each variable can take

G.A = A;
G.N = size(A,1);
G.k = k;

G.b = ones(G.N,k);

G.phi = ones(G.N,k);

[G.R,G.C] = find(A);

for i=1:length(G.R)
    %G.m{i}{j} is the message vector from x_i to x_j
    G.m{G.R(i)}{G.C(i)} = ones(k,1);
    G.psi{G.R(i)}{G.C(i)} = ones(k);
end