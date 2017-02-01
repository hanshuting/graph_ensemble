function H = HTRW(mu, rho, L, N)
% HTRW Compute the TRW entropy in the overcomplete parameterization.
%   H = HTRW(mu, rho, L, N, nEdges)
%
%   mu, rho : Overcomplete belief and reweighting vectors; same dim as mu.
%   L, N    : scalar; number of labels and nodes.

muLogMu = mu .* log(mu);
muLogMu(mu == 0) = 0;

negHnode = (1 - rho(1:L*N)).' * muLogMu(1:L*N);
negHedge = rho(L*N+1:end).'   * muLogMu(L*N+1:end);

H = -negHnode - negHedge;

end
