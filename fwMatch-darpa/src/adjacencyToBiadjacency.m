function B = adjacencyToBiadjacency(A)

    
firstCol = find(A(1,:), 1);
    firstRow = find(A(:,1), 1);
    
    B = A(1:firstRow-1,firstCol:end);

end

