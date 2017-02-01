function [xNode, xEdge, eBelow] = solveQPBO(gNode, gEdge, edges)    
%    gEdge(2:3,:) = gEdge([3 2],:);

    % Coerce edges
    if ~isa(edges, 'int32')
        edges = int32(edges);
    end

    % Dispatch
    if isa(gNode, 'double') && isa(gEdge, 'double')
        [xNode, xEdge, eBelow] = QPBO_double_mex(gNode, gEdge, flipud(edges));
    elseif isa(gNode, 'single') && isa(gEdge, 'single')
        [xNode, xEdge, eBelow] = QPBO_single_mex(gNode, gEdge, flipud(edges));
    else
        error('gNode, gEdges must be both double or single (or int, could support)');
    end

%    keyboard;
    xEdge([2 3],:) = xEdge([3 2],:);

%    eCheck = frobProd(xNode, gNode) + frobProd(xEdge, gEdge);
%    assertElementsAlmostEqual(eBelow, eCheck);

end

