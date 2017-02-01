function [theta beliefs objectives] = fwMinRegFeats(Y, features, lambda, varargin)

p = inputParser;
p.KeepUnmatched = true; % Enable extraneous options.

p.addRequired('Y', @isnumeric);
p.addRequired('features', @iscell);
p.addRequired('lambda', @(x) isnumeric(Y) && isscalar(x));

[M N] = size(Y);

p.addParamValue('beliefs', []);
p.addParamValue('reweight', ones(N, 1));
p.addParamValue('TolGap', 1e-6);
p.addParamValue('maxIter', inf);
p.addParamValue('debug', true);

p.parse(Y, features, lambda, varargin{:});
o = p.Results;

[M2 K] = size(features);
assert(M == M2, 'Y and features must have same number of examples');

if isempty(o.beliefs)
    % TODO: Unify cell/matrix forms
    beliefs = 1 / N * ones(N, N, M);
else
    beliefs = o.beliefs;
end
    

sumGY = computeFeatInnerProduct(features, Y);
sumGT = computeFeatInnerProduct(features, beliefs);

% sumGY and sumGT are MxK arrays. Summing out the rows gives a
% K-vector which we use in computing the gradient.

momentDiffs = sum(sumGY - sumGT, 1);

t = 1;

beliefHist = {};

while t <= o.maxIter && (t <= 2 || dualityGap > o.TolGap)
    % Solve M *separate* matching problems    
    entropy = 0;
    dualityGap = 0;    
    for m = 1:M
        gradMomentDiffs = zeros(N);        
        
        for k = 1:K
            gradMomentDiffs = gradMomentDiffs - momentDiffs(k) * features{m,k};            
        end
        
        grad = gradMomentDiffs / lambda + modOneDEntropyDerivative(beliefs(:,:,m), o.reweight);
%         map_assignment = munkres(grad');
        % Test CSA in place    
        
        map_assignment = csaAssignPerm(grad');        
        
%         assertElementsAlmostEqual(matchingCost(grad', map_assignment), ...
%                                   matchingCost(grad', csa_map_assignment));    
        % End CSA in-place test                  
        
        beliefs(:,:,m) = update_beliefs_deterministic_step(beliefs(:,:,m), map_assignment, N, t);
        
        % Notation in Jaggi (2)
        xGrad = sum(vec(beliefs(:,:,m) .* grad));
        sGrad = sum(vec(grad(expandPerm(map_assignment))));
        
        dualityGap = dualityGap + xGrad - sGrad;
        
%        beliefs = update_beliefs_linesearch(beliefs, map_assignment, X, lambda, o.reweight);

        entropy = entropy + sum(vec(modOneDEntropy(beliefs(:,:,m), o.reweight)));
    end    
    
    sumGT = computeFeatInnerProduct(features, beliefs);
    momentDiffs = sum(sumGY - sumGT, 1);

    beliefHist{t} = beliefs;
    objectives(t) = 1 / (2*lambda) * sum(momentDiffs .^ 2) + entropy;

    t = t + 1;

    if o.debug
%         disp('grad');
%         gradient
%         disp('map');
%         map_assignment
        if t > 2
            fprintf('[ .M t=%d] objectives(end) = %g; dualityGap = %g, TolGap = %g\n', t, objectives(end), dualityGap, o.TolGap);
        end        
    end
end


[minObj, minIdx] = min(objectives);
minSumGT = computeFeatInnerProduct(features, beliefHist{minIdx});
theta    = sum(sumGY - minSumGT, 1)' / lambda;

fprintf('Minimum obj: %g\n', minObj);
fprintf('theta at the min obj:\n');
theta

% objectives(t) = evaluate_objective(beliefs, Y, lambda, o.reweight);

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
% step = 2/(2 + t^(1/2 * exp(-t/1000) + 1/2));
beliefs = (1 - step)*beliefs;

for j = 1:n
    i = map_assignment(j);
    beliefs(i,j) = beliefs(i,j) + step;
end

end

% function beliefs = update_beliefs_linesearch(beliefs, map_assignment, X, lam, reweight)
% 
% expanded_map = expandPerm(map_assignment, @zeros);
% 
% f = @(st) evaluate_objective((1 - st)*beliefs + st*expanded_map, X, lam, reweight);
% 
% MAX_STEP = 1 - 1e-9;
% 
% step = fminbnd(f, 0, MAX_STEP);
% 
% beliefs = (1 - step) * beliefs + step * expanded_map;
% 
% 
% end
% 
