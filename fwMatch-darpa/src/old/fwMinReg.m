function [beliefs objectives] = fwMinReg(X, lambda, varargin)

[N M] = size(X);
assert(N == M, 'A must be square');

p = inputParser;
p.KeepUnmatched = true; % Enable extraneous options.

p.addRequired('X', @isnumeric);
p.addRequired('lambda', @(x) isnumeric(X) && isscalar(x));

p.addParamValue('beliefs', 1/N * ones(N));
p.addParamValue('reweight', ones(N, 1));
p.addParamValue('TolFun', 1e-6);
p.addParamValue('maxIter', inf);
p.addParamValue('debug', false);

p.parse(X, lambda, varargin{:});
o = p.Results;

beliefs = o.beliefs;

t = 1;

while t <= o.maxIter && (t <= 2 || abs(objectives(end-1) - objectives(end)) > o.TolFun)
    
    gradient = (beliefs - X) / lambda + modOneDEntropyDerivative(beliefs, o.reweight);
    map_assignment = munkres(gradient');

    % Test CSA in place    
    [spC, scale] = sparsifyAndRound(gradient');
    edges = csaAssign(2*N, spC);
    csa_map_assignment = edges(2,:) - D;
    
    assertElementsAlmostEqual(matchingCost(gradient', map_assignment), ...
                              matchingCost(gradient', csa_map_assignment));    
    % End CSA in-place test
    
    objectives(t) = evaluate_objective(beliefs, X, lambda, o.reweight);
    beliefs = update_beliefs_deterministic_step(beliefs, map_assignment, N, t);

    t = t + 1;

    if o.debug
        disp('grad');
        gradient
        disp('map');
        map_assignment
        if t > 2
            fprintf('[ .M t=%d] objectives(end-1) = %g; objectives(end) = %g; diff = %g, TolFun = %g\n', t, objectives(end-1), objectives(end), abs(objectives(end-1) - objectives(end)), o.TolFun);
        end        
    end
end

objectives(t) = evaluate_objective(beliefs, X, lambda, o.reweight);

end

function obj = evaluate_objective(beliefs, X, lam, reweight)
sq = sum(vec((X - beliefs).^2)) / (2 * lam);
H  = +sum(vec(modOneDEntropy(beliefs, reweight)));
obj = sq + H;
end

function ent = modOneDEntropy(beliefs, reweight)
reweightMat = bsxfun(@plus, reweight, reweight') - 1;
ent = (beliefs .* log(beliefs)) - (reweightMat .* (1 - beliefs) .* log(1 - beliefs));
ent(beliefs == 0 | beliefs == 1) = 0;
end

function dEnt = modOneDEntropyDerivative(beliefs, reweight)
riPlusRj = bsxfun(@plus, reweight, reweight');
dEnt = riPlusRj + log(beliefs) + (riPlusRj - 1) .* log(1 - beliefs);

% Technically this isn't right but whatever... maybe David Sontag's new
% paper has some way to deal with this?
if any(vec(beliefs == 0 | beliefs == 1))
    error('modOneDEntropyDerivative:boundary', 'Hit the boundary; unbounded gradient');
end

end

%%this just uses the 2/(t+2) rule
function beliefs = update_beliefs_deterministic_step(beliefs, map_assignment, n, t)
step = 2/(2 + t);
beliefs = (1 - step)*beliefs;

for j = 1:n
    i = map_assignment(j);
    beliefs(i,j) = beliefs(i,j) + step;
end

end

function beliefs = update_beliefs_linesearch(beliefs, map_assignment, X, lam, reweight)

expanded_map = expandPerm(map_assignment, @zeros);

f = @(st) evaluate_objective((1 - st)*beliefs + st*expanded_map, X, lam, reweight);

MAX_STEP = 1 - 1e-9;

step = fminbnd(f, 0, MAX_STEP);

beliefs = (1 - step) * beliefs + step * expanded_map;


end

