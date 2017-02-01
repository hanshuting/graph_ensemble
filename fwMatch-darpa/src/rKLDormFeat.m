function features = rKLDormFeat(dormMat, population)
% rKLDormFeat  Feature containing symmetrized KL divergence of dorm
% predictions
%
%   Assume 

N = length(population);

features{1} = zeros(N);

for j = 1:N
    jj = population(j);
    
    for i = 1:(j-1)        
        ii = population(i);        
        
        dijk = symmKL(dormMat(ii,:), dormMat(jj,:));
        features{1}(i,j) = dijk;
    end
end
