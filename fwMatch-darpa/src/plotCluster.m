function plotCluster(A, map)
[edges, withins] = summarizeCluster(A, map);

colors = {'r', 'g', 'b', 'y', 'k'};

K = length(edges);
for k = 1:K
    % color index starts with 2    
    c = k + 2;    
    E = size(edges{k}, 2);

    xs = [edges{k}(1,:) edges{k}(2,:)];
    ys = [edges{k}(2,:) edges{k}(1,:)];
    
    scatter(xs, ys, 
    
end

figure;
imagesc(A);
colorbar;
    
end