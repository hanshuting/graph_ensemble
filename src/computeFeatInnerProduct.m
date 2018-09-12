function FY = computeFeatInnerProduct(features, Y)

    [M, K] = size(features);
    
    % Precompute F^(m)\top Y^(m)    
    FY  = zeros(M, K);                                                
    
    if iscell(Y)
        % Y is a cell of adjacency matrices -- perfect matchings.
        assert(length(Y) == M, 'Length of Y and features mismatch');
        
        for m = 1:M
            % Verify (im)perfect matching
            assert(all(vec(Y{m} == Y{m}')), 'Y{m} is not symmetric');
            assert(all(sum(Y{m}, 1) <= 1),  'Y{m} row sums not <= 1');
            assert(all(sum(Y{m}, 2) <= 1),  'Y{m} col sums not <= 1');
            
            Ym = cast(Y{m}, 'logical');
            
            for k = 1:K                
                FY(m,k) = sum(vec(features{m,k}(Ym)));                
            end            
        end        
    else 
        if length(size(Y)) == 2           
            % Y is a matrix of permutations (each row is a permutation)        

            for m = 1:M
                P = expandPerm(Y(m,:));
                for k = 1:K                                
                    FY(m,k) = sum(vec(features{m,k}(P)));                
                end
            end    
        else
            assert(length(size(Y)) == 3);
            for m = 1:M
                for k = 1:K
                    FY(m,k) = sum(vec(features{m,k} .* Y(:,:,m)));
                end
            end
        end
    end 
    
end
