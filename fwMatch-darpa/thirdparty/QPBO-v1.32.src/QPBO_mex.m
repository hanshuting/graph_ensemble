% [nodeBelief, edgeBelief] = function QPBO_mex(nodeCost, edgeCost, edgeList)
%
%   Computes fractional node labeling of a pairwise binary MRF.
%   Parameters and labels are provided in the canonical overcomplete
%   parameterization, separated into node and edge matrices. In particular,
%   beliefs satisfy the marginalization and consistency constraints of the
%   local polytope.
%
%     nodeBelief : 2 x nNodes matrix of fractional labels; column is node id.
%     edgeBelief : 4 x nEdges matrix of fractional labels; column is edge id
%                  and rows are (x_i, x_j) = 00, 10, 01, 11, e.g. the
%                  column-major storage of a 2x2 matrix where COLUMNS index x_i
%                  and ROWS index x_j.
%     nodeCost, edgeCost : same format as beliefs.
%     edgeList   : 2 x nEdges matrix of edges. Each edge should only appear once.

