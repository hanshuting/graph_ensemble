function [ siVec, sjVec, swVec ] = findUT(W)
% findUT Find the upper triangular of a (sparse) matrix
%   [siVec, sjVec, swVec] = findUT(W) returns the same results as find
%   restricted to the upper triangular of W.

    [iVec, jVec, wVec] = find(W);
    % As always, take the upper triangular
    sel = iVec < jVec;
    siVec = iVec(sel);
    sjVec = jVec(sel);
    swVec = wVec(sel);

end

