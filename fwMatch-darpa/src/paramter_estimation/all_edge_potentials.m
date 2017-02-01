function [thetaE] = all_edge_potentials(edge_potentials,numNodes)
    
    thetaE = zeros(1,numNodes*(numNodes-1)/2);
    [r,c] = find(edge_potentials ~= 0);
    validIndx = [r c];
    validIndx = validIndx(find(validIndx(:,1) < validIndx(:,2)),1:2);
    thetaE(getEdgeId(numNodes,validIndx)) = ...
        edge_potentials(sub2ind([numNodes numNodes],validIndx(:,1),validIndx(:,2)));

end

