function [theta bels objs iterHist exitflag] = fwNewtonMinRegFeats(Y, features, lambda, varargin)

% exitflag can be
%   - dualityGapMet
%   - iterCountMet
%   - stuck
%
% If termination is due to *both duality gap and iteration count
% constraints being met, then the flag will be set to dualityGapMet (the
% more informative condition).

p = inputParser;
p.KeepUnmatched = true; % Enable extraneous options.

p.addRequired('Y', @isnumeric);
p.addRequired('features', @iscell);
p.addRequired('lambda', @(x) isnumeric(x) && isscalar(x));

[M N] = size(Y);
Nsq = N^2;

p.addParamValue('beliefs', []);
p.addParamValue('reweight', ones(N, 1));
p.addParamValue('TolGap', 1e-10);
p.addParamValue('MaxIter', inf);
p.addParamValue('MaxFunEvals', 10000);
p.addParamValue('debug', false);
p.addParamValue('newtonStep', true);
p.addParamValue('lineSearch', false);
p.addParamValue('plotLine', false);
p.addParamValue('checkStuck', false);
p.addParamValue('errNegDualityGap', false);
p.addParamValue('awayStep', false);

% lineSearchOpts = optimset('TolX', 1e-12, 'Display', 'iter');
lineSearchOpts = optimset('TolX', 1e-9);

p.parse(Y, features, lambda, varargin{:});
o = p.Results;

if o.plotLine
    % TODO: Make it bigger
    pHt = figure;    
end

[M2 K] = size(features);
assert(M == M2, 'Y and features must have same number of examples');

% Construct the big G and Y matrices.

G = cell2mat(cellfun(@vec, features, 'UniformOutput', false));

Ymats = cell(M, 1);
for m = 1:M
    Ymats{m} = expandPerm(Y(m,:), @zeros);    
end

y = vertcat(cell2mat(cellfun(@vec, Ymats, 'UniformOutput', false)));

if isempty(o.beliefs)
    tau = 1/N * ones(M*N^2, 1);
else
    tau = o.beliefs(:);
end

Rsum = bsxfun(@plus, o.reweight, o.reweight');
rsum = repmat(vec(Rsum), M, 1);

function H = hessFcn(tt, ~)
    % lam is the vector of Lagrange multipliers. Not used since we have
    % only linear constraints.
        
    Dtt = diag((1 - rsum.*tt) ./ (tt .* (1 - tt)));
    H = 1/lambda * G*G' + Dtt;
    H2 = G * (1/lambda * eye(size(G, 2))) * G' + Dtt;
    
    assert(norm(H - H2, 1) == 0);
end

function H = negBetheEntropy(tt)
    H = sum(tt.*log(tt) - (rsum - 1).*(1 - tt).*log(1 - tt));
end

GtTimesY = G'*y;
function [f, g, GtTimesYMinusTau] = obj(tt)    
    GtTimesYMinusTau = GtTimesY - G'*tt;
    
    sq = 1/(2*lambda) * sum(GtTimesYMinusTau.^2);
    
    f  = sq + negBetheEntropy(tt);            
    g  = -(G * GtTimesYMinusTau)/lambda + rsum + log(tt) + (rsum - 1).*log(1 - tt);    
end

function [f, ng] = naiveNewtonStep(tt)
    [f, g] = obj(tt);
    ng = hessFcn(tt, []) \ g;
end

function [f, ng, GtTimesYMinusTau] = fastNewtonStep(tt)
    [f, g, GtTimesYMinusTau] = obj(tt);
    
    % Deal with multiplication by diagonal efficiently.
    % See writeup for explanation of what the symbols mean.
    invdtt = (tt .* (1 - tt)) ./ (1 - rsum.*tt);    
       
    % Dinvgrad is MN^2 x 1, DinvG is MN^2 x K.
    Dinvgrad = invdtt .* g;
    DinvG    = bsxfun(@times, invdtt, G);    
    C = lambda*eye(K) + G'*DinvG;
    
    % Parenthesize from right the left, which is optimal.
    bigTerm = DinvG * (C \ (G' * Dinvgrad));
    
    ng = Dinvgrad - bigTerm;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up away step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CAUTION: Think about the memory usage when running this for large
% problems.

% TODO: Smarter data structures (that nevertheless permit good vectorized
% Mv multiply -- think about CSC sparse matrix?)
%
% Or perhaps just dynamically resize as needed.
if o.awayStep    
    nFWSteps   = 0;
    nAwaySteps = 0;
    
    maxVerts = o.MaxIter + 1;
    
    verts  = zeros(M*N^2, maxVerts);
    active = false(1, maxVerts);
    % Used in the notation (vertInds(active))(i) if Matlab allowed
    % recursive indexing; it doesn't so we split into two lines of code.
    vertInds = 1:(M*N^2);        
        
    % \alpha in Lacoste-Julien
    weights = zeros(1, maxVerts);
    
    % The initial condition
    verts(:,1) = tau;
    active(1)  = true;
    weights(1) = 1;    
end

% This should be defined in the above block but Matlab is a sucky
% programming language and doesn't let you define lambdas inside control
% structures.

function [v, w, i] = solveAwayLP(gr)
    % [v, w, f] = solveAwayLP(gr) returns the vertex and its weight in an
    % away step.
    
    awayVals = verts(:,1:active)' * gr;
    [~, awayAVIdx] = min(awayVals);

    % Really want awayIdx = (vertInds(active))(awayAIdx) but we can't
    % use that notation.
    avIdxs = vertInds(active);
    i = avIdxs(awayAVIdx);

    v = verts(:,i);
    w = weights(:,i);    
end


t = 1;
dualityGap = inf;
while t <= o.MaxIter && dualityGap > o.TolGap    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute the objective and gradient.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if o.newtonStep
%         [fNaive, ngNaive] = naiveNewtonStep(tau);
        [f, g, GtTimesYMinusTau] = fastNewtonStep(tau);
        
        % assertElementsAlmostEqual(f, fNaive)                
        % assertElementsAlmostEqual(ng, ngNaive)
    else
        [f, g, GtTimesYMinusTau] = obj(tau);
    end
    tauMAP = zeros(M*N^2, 1);
    objs(t) = f;
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute MAP assignments.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Break up the gradient into per-sample gradients
    for m = 1:M
        mIdxBegin = Nsq * (m - 1) + 1;
        mIdxEnd   = Nsq * m;
        idxRange  = mIdxBegin:mIdxEnd;        
        
        grad      = reshape(g(idxRange), N, N);
        
%         mapAssign = csaAssignPerm(grad');
        permMat = csaAssignPermMat_mex(grad);
%         permMat = expandPerm(munkresOpt(grad'), @zeros);
        tauMAP(idxRange) = permMat(:);        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute duality gap (could be difficult for block-coordinate)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    dFW = tauMAP - tau;
    dualityGap = -sum(g .* dFW);
    
    if dualityGap < 0
        if o.errNegDualityGap
            error('fwNewtonMinRegFeats:negDualityGap', 'Duality gap is negative!');
        else
            warning('fwNewtonMinRegFeats:negDualityGap', 'Duality gap is negative!');
        end
    end

    iterHist.dualityGaps(t) = dualityGap;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute away step (Lacoste-Julien 2014; arXiv:1312.7864v2 Alg. 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    goAwayStep = false;
    if o.awayStep && nFWSteps > 0
        [awayVert, awayWeight, awayIdx] = solveAwayLP(g);
        
        % Compare projections of the FW and away directions        
        dAway = tau - awayVert;        
        
        if sum(g .* dFW) <= sum(g .* dAway)
            goAwayStep = true;
        end
    end
    
    if goAwayStep
        dir = dAway;
        maxStep = awayWeight;        
    else
        dir = dFW;
        maxStep = 1;        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Linesearch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    if o.lineSearch
%         % For linesearch. See writeup for what these mean.
        GtTimesYMinusTauSq    = sum(GtTimesYMinusTau .^ 2);
%         GtTimesTauMinusTauMAP = G' * (tau - tauMAP);
% 
%         linCoeff  = 2 * sum(GtTimesYMinusTau .* GtTimesTauMinusTauMAP);
%         quadCoeff = sum(GtTimesTauMinusTauMAP .^ 2);    
% 
%         % MATLAB forces me to write an anonymous one-line function.
%         ht = @(eta) 1/(2*lambda) * (GtTimesYMinusTauSq + linCoeff.*eta + quadCoeff.*eta.^2) + ...
%             negBetheEntropy(tau + eta*dir);
        
        % Linesearch in arbitrary direction dir.
        GtTimesDir   = G'*dir;
        dirLinCoeff  = -2 * sum(GtTimesYMinusTau .* GtTimesDir);
        dirQuadCoeff = sum(GtTimesDir .^ 2);
        
        dirHt = @(eta) 1/(2*lambda) * (GtTimesYMinusTauSq + dirLinCoeff.*eta + dirQuadCoeff.*eta.^2) + ...
            negBetheEntropy(tau + eta*dir);
        
%         htSlow = @(eta) obj((1 - eta)*tau + eta*tauMAP);
%         
%         xs = linspace(0, 0.99, 300);
%         yht = zeros(size(xs));
%         yhtSlow = zeros(size(xs));
%                 
%         for i = 1:length(xs)
%             yht(i) = ht(xs(i));
%             yhtSlow(i) = htSlow(xs(i));
%         end
%         
%         relErr = abs(yht - yhtSlow) ./ max(abs(yht), abs(yhtSlow));
%         figure; plot(xs, yht, xs, yhtSlow); title(sprintf('Iter %d Linesearch Values', t));
%         figure; plot(xs, relErr); title(sprintf('Iter %d Relative Errors', t));        
        % TODO: Track function calls (really just negBetheEntropy calls; we
        % can assume scalar parts are trivial)
%         [CHECKstepSz, CHECKhtMin, CHECKexitflag] = fminbnd(ht, 0, maxStep, lineSearchOpts);
        [stepSz, htMin, exitflag] = fminbnd(dirHt, 0, maxStep, lineSearchOpts);        
        assert(exitflag == 1, 'fminbnd did not converge.');
        
%         assertElementsAlmostEqual(stepSz, CHECKstepSz);
%         assertElementsAlmostEqual(htMin, CHECKhtMin);

        if o.checkStuck
            % Exception: Don't move.
            if ht(0) < htMin
                warning('ht(0) was smaller than htMin: %g < %g', ht(0), htMin);
                stepSz = 0;
                htMin = ht(0);
            end           
        end
                
        iterHist.stepSzs(t) = stepSz;                
        
        if o.plotLine        
            figure(pHt);
            subplot(1,2,1);
            ax = gca;
            hold on;
            ezplot(ax, ht, [0, 1]);            
            plot(ax, stepSz, htMin, '-.or');
%             plot(ax, stepSz, ht(stepSz), '-.or');
            title(sprintf('Iter = %d, min = (%g, %g)', t, stepSz, htMin));
%             title(sprintf('Iter = %d, min = (%g, %g)', t, stepSz, ht(stepSz)));           
            hold off;
            
            subplot(1,2,2);
            hist([iterHist.stepSzs], 100);
            title('Step size histogram')
        end                                
    else
        % Deterministic stepsize.
        stepSz = 2 / (2 + t);
        iterHist.stepSzs(t) = stepSz;        
    end
    
    % If we didn't actually move, quit.
    if o.checkStuck && stepSz == 0
        exitflag = 'stuck';
        tauFinal = tau;
        theta = 1 / lambda * G' * (y - tauFinal);
        bels = reshape(tauFinal, N, N, M);
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Actual update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % Update the dense iterate
    tau = tau + stepSz*dir;            
    
    % Update the away step bookkeepers
    if o.awayStep
        if goAwayStep
            if o.debug
                fprintf('Away Step: maxStep - stepSz = %g\n', maxStep - stepSz);
            end

            if stepSz == maxStep
                % Drop step: Remove the vertex we went in the away step
                % direction from the set of active vertices.
                active(awayIdx) = false;
            end

            % Update alphas
            weights(awayIdx) = (1 + stepSz)*weights(awayIdx) - stepSz;

            otherActives = active;
            otherActives(awayIdx) = false;
            weights(otherActives) = (1 + stepSz)*weights(otherActives);        
        else
            % Took the FW step. We will need to add a vertex to our active
            % set.
            %
            % By convention, the vertex added at iteration t is at column
            % t+1.
            
            if stepSz == 1
                warning('Went to a vertex; should not happen for entropy!');
                active(1:end) = false;
            end
            
            % Whether or not we went to all the way to the vertex tauMAP,
            % we moved in that direction somewhat, so add to our active
            % set.
%             active(t+1) = true;
%             weights(t+1) = (1 - stepSz)*weights(
        end
    end
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Linesearch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % We have computed tauMAP now, so 
    if o.awayStep
        % Record the current extreme point.
        verts(:,t) = tauMAP;
        active(t) = true;
        steps(t)  = stepSz;
        
%         % Compute the away step
%         awayObjs = verts(:,1:        
    end
    
        
    if o.debug
        fprintf('[ .M t=%d] objectives(end) = %g; dualityGap = %g, norm(ng, 1) = %g, stepSz = %g, TolGap = %g\n', t, objs(end), dualityGap, norm(g, 1), stepSz, o.TolGap);        
    end
    
    t = t + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Report stopping condition and final answers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

if dualityGap <= o.TolGap
    exitflag = 'dualityGapMet';
elseif t == o.MaxIter + 1
    exitflag = 'maxIterMet';
else
    error('fwNewtonMinRegFeat:terminate', 'Unreachable termination condition.');
end
    
tauFinal = tau;
theta = 1 / lambda * G' * (y - tauFinal);
bels = reshape(tauFinal, N, N, M);

end


