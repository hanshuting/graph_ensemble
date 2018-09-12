function e = energyOvercomplete(costNode, costEdge, xNode, xEdge)
    c = [costNode(:) ; costEdge(:)];
    x = [xNode(:) ; xEdge(:)];

    e = sum(c .* x);
end

