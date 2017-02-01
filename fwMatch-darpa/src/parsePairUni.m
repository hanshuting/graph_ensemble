function [Fsquare, FnegExpNegSquare] = parsePairUni(simpleFile)
% parsePair  Construct unipartite feature tensor from simplified pairing file
%
%   [Fsquare, FexpSquare] = parsePair(simpleFile) returns a DxDxKxM feature tensor
%   (bipartite) for the simplified pairing file simpleFile. simpleFile has
%   one number per line, zero-prefixed padded, and every pair of lines
%   corresponds to a pair.
%
%   Fsquare is intended for use with learning; FexpSquare for use without learning
%   (set theta = 1 uniformly). Follows Caetano09.

    dataRoot = 'thirdparty/graphmatchBMRM/Data/house';
    pth = sprintf('%s/pairs/%s', dataRoot, simpleFile);

    fid = fopen(pth);
    C = textscan(fid, '%s');
    fclose(fid);    

    M = length(C{1}) / 2;
    pairs = reshape(C{1}, 2, M);
    
    K = size(importdata(sprintf('%s/house%s.scf', dataRoot, pairs{1, 1})), 2);
    
    Fsquare          = cell(M, K);
    FnegExpNegSquare = cell(M, K);
    for m = 1:M
        file1 = sprintf('%s/house%s.scf', dataRoot, pairs{1, m});
        file2 = sprintf('%s/house%s.scf', dataRoot, pairs{2, m});
        
        X1 = importdata(file1);
        X2 = importdata(file2);

        for k = 1:K
            Fsquare{m,k} = squareDiffs(X1, X2, k);
            FnegExpNegSquare{m,k} = -exp(-Fsquare{m,k});
        end
    end        
end

function Funi = squareDiffs(X1, X2, k)
    [D , K]  = size(X1);
    [D2, K2] = size(X2);
    assert(D == D2 && K == K2, 'X1 and X2 differ in size');
    
    Fbi = zeros(D, D);
    for i = 1:D
        for j = 1:D
            Fbi(i,j) = (X1(i,k) - X2(j,k)) .^ 2;            
        end
    end
    
    Funi = biadjacencyToAdjacency(Fbi);    
end    

