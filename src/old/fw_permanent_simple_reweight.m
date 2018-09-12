function [beliefs objectives] = fw_permanent_simple(beliefs, A, reweight, TolFun, isDebug, maxIter)

[n m] = size(A);
assert(n == m,'input must be square');

%initialization: uniform
if isempty(beliefs)
    beliefs = (1/n)*ones(n,n);
end

if nargin < 3
    reweight = ones(n, 1);
else
    % Force to be column matrix for correct bsxfun
    reweight = reweight(:);
end

if nargin < 4
    TolFun = 1e-6;
end

if nargin < 5
    isDebug = false;
end

if nargin < 6
    maxIter = inf;
end

%%We seek min_{beliefs} -sum_ij beliefs_ij A_ij  + sum_ij beliefs_ij log(beliefs_ij) - sum_ij (1 - beliefs_ij) log(1 - beliefs_ij)

%objectives = zeros(1,numIters);
t=1;
%objectives(t) = evaluate_objective(beliefs,A,n);

while t <= maxIter && (t <= 2 || abs(objectives(end-1) - objectives(end)) > TolFun)
    
    gradient = -A + modOneDEntropyDerivative(beliefs, reweight);
    map_assignment = lap_mex(gradient');
    
    if isDebug        
        disp('grad');
        gradient
        disp('map');
        map_assignment
    end
    
%     assert(length(map_assignment) == n && all(map_assignment > 0));
    objectives(t) = evaluate_objective(beliefs, A, reweight);
    beliefs = update_beliefs_deterministic_step(beliefs, map_assignment, reweight, n, t);
    t = t + 1;

%     if t > 2
%         fprintf('[ .M t=%d] objectives(end-1) = %g; objectives(end) = %g; diff = %g, TolFun = %g\n', t, objectives(end-1), objectives(end), abs(objectives(end-1) - objectives(end)), TolFun);
%     end
end

objectives(t) = evaluate_objective(beliefs, A, reweight);

end

function obj = evaluate_objective(beliefs, A, reweight)
U = -sum(vec(beliefs.*A));
H = +sum(vec(modOneDEntropy(beliefs, reweight)));
obj = U + H;
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
function beliefs = update_beliefs_deterministic_step(beliefs, map_assignment, reweight, n, t)
step = 2/(2 + t);
beliefs = (1 - step)*beliefs;

% disp('post de-step, pre step');
% beliefs

for j = 1:n
    i = map_assignment(j);
    beliefs(i,j) = beliefs(i,j) + step;
end

% disp('post step');
% beliefs

end


