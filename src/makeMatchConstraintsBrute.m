function st = makeMatchConstraintsBrute(N)
%makeReducedMatchingConstraintsBrute  Brute-force create lin.indep. contrs
%
%   N - number of elements matched
%   st - structure with Aeq, beq constraint; Aeq is the sparse matrix, beq is RHS.
%
%   NOTE: rref may not be necessary (and may be really slowing things down)

Nsq = N * N;
inds = reshape(1:Nsq, N, N);

% badAeq = zeros(2*N, Nsq);

% nnz...
badAeqIvec = zeros(2*Nsq, 1);
badAeqJvec = zeros(2*Nsq, 1);
badAeqWvec = ones(2*Nsq, 1);
badbeq = ones(2*N, 1);

% For each column j, \sum_{ij} B_{ij} = 1.
k = 1;
for j = 1:N
    nInds = length(inds(:,j));
    range = k:k+nInds-1;
    
    badAeqIvec(range) = j;
    badAeqJvec(range) = inds(:,j);
    
    k = range(end) + 1;
end

% For each row i, \sum_{ij} B_{ij} = 1
for i = 1:N
    nInds = length(inds(i,:));
    range = k:k+nInds-1;
    
    badAeqIvec(range) = N + i;
    badAeqJvec(range) = inds(i,:);
    
    k = range(end) + 1;
end

badAeq = sparse(badAeqIvec, badAeqJvec, badAeqWvec, 2*N, Nsq);
[Aeq, beq] = reduceEqualityConstraints(badAeq, badbeq);

st = var2struct(Aeq, beq);

end

