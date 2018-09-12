function g = gradHTRW(mu, rho, L, N)
% HTRW Compute gradient of TRW entropy in the overcomplete parameterization.
%   H = HTRW(mu, rho, L, N, nEdges)
%
%   mu, rho : Overcomplete belief and reweighting vectors; same dim as mu.
%   L, N    : scalar; number of labels and nodes.

g = zeros(size(mu));
onePlusLogMu = 1 + log(mu);
g(1:L*N)     = -(1 - rho(1:L*N)) .* onePlusLogMu(1:L*N);
g(L*N+1:end) = -rho(L*N+1:end)   .* onePlusLogMu(L*N+1:end);

end
