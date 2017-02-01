function P = expandPerm(perm, factory)
% expandPerm  Expand permutation list to permutation matrix (column).
%
%   P = expandPerm(perm, factory) expands the linear permutation perm to a
%   permutation matrix, whose zeros are created by factory.
%
%   expandPerm(p, @zeros) produces a double matrix (suitable for averaging)
%   while expandPerm(p, @false) produces a logical matrix (suitable for
%   indexing), and expandPerms(p, @sparse) produces a sparse matrix.
%
%   Default is @false (logical).

    if nargin == 1
        factory = @false;
    end

    D = length(perm);
    P = factory(D, D);
    inds = sub2ind([D D], perm, 1:D);

    P(inds) = 1;                

end

