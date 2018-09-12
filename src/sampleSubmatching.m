function [Asub, YTrueSub] = sampleSubmatching(A, YTrue, n)
% sampleSubmatching  Return a submatrix A and YTrue containing the true
%                    matching.
%
%   [Asub, YTrueSub] = sampleSubmatching(A, YTrue, n)

    assert(mod(n, 2) == 0, 'n must be even.');

    [aiVec, ajVec, ~] = findUT(A);
    [yiVec, yjVec, ~] = findUT(YTrue);
    
    E = length(yiVec);
    
    eixs = randsample(E, n);    
    
    subEdges = [yiVec(eixs) yjVec(eixs)];    
    subNodes = sort(vec(subEdges));
%   ^ Mapping from i to iBig
    
    Asub     = zeros(2*n);
    YTrueSub = zeros(2*n);
    for j = 1:length(subNodes)
        jBig = subNodes(j);
        for i = 1:length(subNodes)
            iBig = subNodes(i);
            
            if iBig < jBig            
                Asub(i,j) = A(iBig, jBig);
                YTrueSub(i,j) = YTrue(iBig, jBig);
            end
        end
    end
    
    % Double check
    
    
%     assert(all(sum(YTrueSub, 1) == 1) && ...
%            all(sum(YTrueSub, 2) == 1));

end

