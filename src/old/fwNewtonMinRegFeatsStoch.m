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

t = 1;
dualityGap = inf;
while t <= o.MaxIter && dualityGap > o.TolGap
    
    % Compute the big Newton step
    
    % TODO: Possibly optimize to only compute gradient for one sample.
    [f, ng, ~] = obj(tau);    
    objs(t) = f;

    m = randi(M);    
    mIdxBegin = Nsq * (m - 1) + 1;
    mIdxEnd   = Nsq * m;
    idxRange  = mIdxBegin:mIdxEnd;        

    grad      = reshape(ng(idxRange), N, N);

%         mapAssign = csaAssignPerm(grad');
    permMat = csaAssignPermMat_mex(grad);
    tauMIdxsMAP = permMat(:);
    
    % Compute update
    if o.lineSearch
        % HACK: Quick tau solutions
        tauSub = zeros(size(tau));
        tauSub(idxRange) = tau(idxRange);
        tauMAP = zeros(size(tau));
        tauMAP(idxRange) = tauMIdxsMAP;                
        
        
        % For linesearch. See writeup for what these mean.
        [~, ~, GtTimesYMinusTau] = obj(tauSub);
        GtTimesYMinusTauSq    = sum(GtTimesYMinusTau .^ 2);
        GtTimesTauMinusTauMAP = G' * (tauSub - tauMAP);

        linCoeff  = 2 * sum(GtTimesYMinusTau .* GtTimesTauMinusTauMAP);
        quadCoeff = sum(GtTimesTauMinusTauMAP .^ 2);    

        % MATLAB forces me to write an anonymous one-line function.
        ht = @(eta) 1/(2*lambda) * (GtTimesYMinusTauSq + linCoeff.*eta + quadCoeff.*eta.^2) + ...
            negBetheEntropy((1 - eta).*tauSub + eta.*tauMAP);
        
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
        
        % TODO: Steal code from FW w/ Adrian to plot stuff here.
        % TODO: Track function calls (really just negBetheEntropy calls; we
        % can assume scalar parts are trivial)
        [stepSz, htMin, exitflag] = fminbnd(ht, 0, 1.0, lineSearchOpts);
        assert(exitflag == 1, 'fminbnd did not converge.');

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
%         stepSz = 2 / (2 + t);
        stepSz = 2*M / (2*M + t);
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
    
    % And do the updat.
    tau(idxRange) = (1 - stepSz)*tau(idxRange) + stepSz*tauMIdxsMAP;
        
    if o.debug && mod(t, 100) == 0
        fprintf('[ .M t=%d] objectives(end) = %g; stepSz = %g, TolGap = %g\n', t, objs(end), stepSz, o.TolGap);        
    end
    
    t = t + 1;
end

% Why did we terminate?

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


