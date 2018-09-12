function [g, scale] = sparsifyAndRound(C)
%sparsifyAndRound  Convert biadjacency cost matrix to csaAssign format

    [Cround, scale] = nonnegRoundInt(C, 'int32');
    g = sparsify(Cround);

end

