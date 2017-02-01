function [logZ] = solveConditionalDAI(F,G,edge_list,varargin)

% [mu, xi, twoMargs, energy, complexity, exact, logZ] = solveDAI(theta, W, varargin)
%
%   Same interface as solveBetheExact.

% solveDAI Find marginals of problem in Eq 1 by libDAI's methods
%   [ logZ, joint, oneMarginals, twoMarginals, daiTime ] = solveDAI(theta, W, method, daiOpts)
%
%   theta - unary potentials
%   W     - (Sparse) pairwise interactions; only upper triangular taken
%
%   logZ         - exact (JTREE only) or approximate log partition function
%   joint        - final joint distribution
%   oneMarginals - nNodes x 1 vector of P(x_n = 1)
%   twoMarginals - 2 x 2 x nEdges array of 2x2 matrices where
%                  M(qi,qj) = P(x_i = qi, x_j = qj).
%                  We only compute pairwise marginals for edges that appear
%                  in W.

    p = inputParser;
    p.addRequired('F', @isnumeric);
    p.addRequired('G', @isnumeric);
    p.addRequired('edge_list', @isnumeric);
    
    % TODO: Better dispatch?
    p.addParamValue('method',    'JTREE');
    p.addParamValue('JTREEOpts', '[updates=HUGIN,verbose=0]');
    p.addParamValue('BPOpts',    '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');    
    p.addParamValue('HAKOpts',   '[doubleloop=1,clusters=BETHE,init=UNIFORM,tol=1e-9,maxiter=10000]');
    
    %not working, maybe try to fix later, see https://github.com/probml/pmtk3/blob/4582f0eed9acd65691681256f9c895769b439f81/toolbox/GraphicalModels/inference/libdai/libdaiOptions.m
    %p.addParamValue('CBPOpts', '[max_levels=12,updates=SEQMAX,tol=1e-9,rec_tol=1e-9,maxiter=500,choose=CHOOSE_RANDOM,recursion=REC_FIXED,clamp=CLAMP_VAR,min_max_adj=1.0e-9,bbp_cfn=CFN_FACTOR_ENT,rand_seed=0,bbp_props=[tol=1.0e-9,maxiter=10000,damping=0,updates=SEQ_BP_REV],clamp_outfile=]');

    p.addParamValue('restarts', 1);
    p.addParamValue('parallelRestarts', false);
    
    p.parse(F,G,edge_list,varargin{:});    
    
    o = p.Results;

    vars = edge_list;
    nNodes = size(F,2);
    nEdges = size(G,2);
    
    % Convert the energy functional to a factor graph.
    psi     = cell(nNodes + nEdges, 1);
    varsets = cell(nNodes + nEdges, 1);
    for n = 1:nNodes        
        psi{n} = struct('Member', n, 'P', [exp(F(1,n)) exp(F(2,n))]');
        varsets{n} = n;
    end
    
    for ne = 1:nEdges               
                
        P = ones(2,2);
        P(1,1) = exp(G(1,ne));
        P(1,2) = exp(G(3,ne));
        P(2,1) = exp(G(2,ne));
        P(2,2) = exp(G(4,ne));
        
        psi{ne + nNodes} = struct('Member', vars(ne,:), 'P', P);
        varsets{ne + nNodes} = vars(ne,:);
    end
    
    % Call DAI    
    %daiMethod  = [ '[' o.method ']' ];
    daiOpts = o.([o.method 'Opts']);
    
    mds = cell(o.restarts, 1);
    oneMStructs = cell(o.restarts, 1);
    twoMStructs = cell(o.restarts, 1);

    if o.parallelRestarts
        parfor t = 1:o.restarts
            [logZs(t),qs{t},mds{t},oneMStructs{t},twoMStructs{t}] = dai(psi, o.method, daiOpts);
        end
    else
        for t = 1:o.restarts
            [logZs(t),qs{t},mds{t},oneMStructs{t},twoMStructs{t}] = dai(psi, o.method, daiOpts);
        end
    end

    [~, iBest] = max(logZs);

    logZ = logZs(iBest);
    
end

