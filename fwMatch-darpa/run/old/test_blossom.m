dbmex;
startup
N = 4000;
A = rand(N,N);
mask = rand(N,N) > 0.5;
A(mask) = 0;
%[as, cs] = blossom_sparse_mex(sparse(A));
[as, cs] = blossom(A);
