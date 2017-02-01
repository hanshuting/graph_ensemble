function [one, two] = unpackOvercomplete(x, L, N, nEdges)
% unpackOvercomplete  De-linearize the overcomplete representation
%
%   [one, two] = unpackOvercomplete(x, N, nEdges)
%
%   See sparseMargExtreme

    one = reshape(x(1:(L*N)), L, N);
    two = reshape(x((L*N+1):end), L, L, nEdges);    

end

