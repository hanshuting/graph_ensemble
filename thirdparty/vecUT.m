function v = vecUT(W, dense)
% Vectorize the upper triangular of W.
%
%   v = vecUT(W) returns a sparse vector v
%   v = vecUT(W, true) returns a dense vector v.

[N, N2] = size(W);
assert(N == N2, 'W must be square');

M = N * (N - 1) / 2;

% Important: must be column vector.
v = zeros(M, 1);
vidx = 1;
for col = 2:N
    lastRow = col - 1;
    vrange  = vidx:(vidx + lastRow - 1);
    v(vrange) = W(1:lastRow,col);
    vidx = vrange(end) + 1;
end

if nargin < 2 || ~dense
    v = sparse(v);
end

end