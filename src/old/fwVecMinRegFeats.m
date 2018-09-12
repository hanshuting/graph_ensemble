function [theta beliefs objs ] = fwNewtonMinRegFeats(Y, features, lambda, varargin)

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
p.addParamValue('matchConstraints', []);
p.addParamValue('boundConstraints', true);
p.addParamValue('stopIfSingular', false)

p.parse(Y, features, lambda, varargin{:});
o = p.Results;

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

function H = hessFcn(tt, lam)
    % lam is the vector of Lagrange multipliers. Not used since we have
    % only linear constraints.
        
    Dtt = diag((1 - rsum.*tt) ./ (tt .* (1 - tt)));
    H = 1/lambda * G*G' + Dtt;
end

function [f, g] = obj(tt)
    Gytt = G'*(y - tt);
    
    sq = 1/(2*lambda) * sum(Gytt.^2);
    H  = sum(tt.*log(tt) - (rsum - 1).*(1 - tt).*log(1 - tt));
    f  = sq + H;
    
    g  = -1/lambda*G*Gytt + rsum + log(tt) + (rsum - 1).*log(1 - tt);    
end

function [f, ng] = naiveNewtonStep(tt)
    [f, g] = obj(tt);
    ng = hessFcn(tt, []) \ g;
end

function [f, ng] = fastNewtonStep(tt)
    [f, g] = obj(tt);
    
    % Deal with multiplication by diagonal efficiently.
    % See math for explanation of what the symbols mean.
    invdtt = (tt .* (1 - tt)) ./ (1 - rsum.*tt);
    invDtt = diag(invdtt);    
       
    % Dinvgrad is MN^2 x 1, DinvG is MN^2 x K.
    Dinvgrad = invdtt .* g;
    DinvG    = bsxfun(@times, invdtt, G);
    C = lambda + G'*DinvG;
    
    % Parenthesize from right the left, which is optimal.
    bigTerm = DinvG * (C \ (G' * Dinvgrad));
    
    ng = Dinvgrad - bigTerm;    
end

function updateTau(stepSz, idxRange, mapAssign)
    tau(idxRange) = (1 - stepSz) * tau(idxRange);

    for c = 1:N
        r = mapAssign(c);
        ind = sub2ind([N, N], r, c) + idxRange(1) - 1;
        
        tau(ind) = tau(ind) + stepSz;
    end    
end

t = 1;
while t < o.MaxIter
    stepSz = 2 / (2 + t);
    
    % Compute the big Newton step
%     [f, ng] = naiveNewtonStep(tau);
    [f, gg] = obj(tau);
    objs(t) = f;
    
%     [fFast, ngFast] = fastNewtonStep(tau);
        
%     assertElementsAlmostEqual(f, fFast);
%     assertElementsAlmostEqual(ng, ngFast);
        
    dualityGap = 0;
    % Break up the gradient into per-sample gradients
    for m = 1:M
        mIdxBegin = Nsq * (m - 1) + 1;
        mIdxEnd   = Nsq * m;
        idxRange  = mIdxBegin:mIdxEnd;
                
        grad      = reshape(gg(idxRange), N, N)';
        
        mapAssign = csaAssignPerm(grad');
        
        % Update the duality gap BEFORE THE UPDATE (since the objectives
        % were pre-update as well)
        % Notation in Jaggi (2)
        xGrad = sum(tau(idxRange) .* ng(idxRange));
        sGrad = sum(vec(grad(expandPerm(mapAssign))));
        
        dualityGap = dualityGap + xGrad - sGrad;        
        
        % Step.
        updateTau(stepSz, idxRange, mapAssign);                
    end    
    if o.debug
%         disp('grad');
%         gradient
%         disp('map');
%         map_assignment
        if t > 2
            fprintf('[ .M t=%d] objectives(end) = %g; dualityGap = %g, TolGap = %g\n', t, objs(end), dualityGap, o.TolGap);
        end        
    end
    
    t = t + 1;
end
    
tauFinal = tau;

theta = 1 / lambda * G' * (y - tauFinal);

beliefs = reshape(tauFinal, N, N, M);
end


