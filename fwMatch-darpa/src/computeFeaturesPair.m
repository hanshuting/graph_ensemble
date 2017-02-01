function [Fsquare, FnegExpNegSquare] = computeFeaturesPair(baseName, setName)
% parsePair  Construct (and save) bipartite features from simplified pairing file
%
%   [Fsquare, FexpSquare] = parsePair(setName) returns a MxK cell of DxD
%   (bipartite) feature matrices for the simplified pairing file simpleFile. simpleFile has
%   one number per line, zero-prefixed padded, and every pair of lines
%   corresponds to a pair.
%
%   setName is the <NUMBER>_<train|test|valid> part of the filename.
%
%   Fsquare is intended for use with learning; FexpSquare for use without learning
%   (set theta = 1 uniformly). Follows Caetano09.

    dataRoot = sprintf('thirdparty/graphmatchBMRM/Data/%s', baseName);
    pth = sprintf('%s/pairs/%ss%s.txt.simple', dataRoot, baseName, setName);

    fid = fopen(pth);
    C = textscan(fid, '%s');
    fclose(fid);    

    M = length(C{1}) / 2;
    pairs = reshape(C{1}, 2, M);
    
    K = size(importdata(sprintf('%s/%s%s.scf', dataRoot, baseName, pairs{1, 1})), 2);
    
    Fsquare          = cell(M, K);
    FnegExpNegSquare = cell(M, K);
    for m = 1:M
        file1 = sprintf('%s/%s%s.scf', dataRoot, baseName, pairs{1, m});
        file2 = sprintf('%s/%s%s.scf', dataRoot, baseName, pairs{2, m});
        
        X1 = importdata(file1);
        X2 = importdata(file2);

        for k = 1:K
            Fsquare{m,k} = bsxfun(@minus, X1(:,k), X2(:,k)') .^ 2;
            FnegExpNegSquare{m,k} = -exp(-Fsquare{m,k});
        end
    end        
    
    saveFile = sprintf('cp/imagedata/%s/%s.mat', baseName, setName);
    save(saveFile, 'Fsquare', 'FnegExpNegSquare');
    
end
