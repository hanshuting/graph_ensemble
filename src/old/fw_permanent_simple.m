function [beliefs, objectives, debugLog] = fw_permanent_simple(beliefs, A, TolFun, isDebug, maxIter)

[n m] = size(A);
assert(n == m,'input must be square');

%initialization: uniform
if isempty(beliefs)
    beliefs = (1/n)*ones(n,n);
end

if nargin < 3
    TolFun = 1e-6;
end

if nargin < 4
    isDebug = false;
end

if nargin < 5
    maxIter = inf;
end

%%We seek min_{beliefs} -sum_ij beliefs_ij A_ij  + sum_ij beliefs_ij log(beliefs_ij) - sum_ij (1 - beliefs_ij) log(1 - beliefs_ij)

%objectives = zeros(1,numIters);
t=1;
%objectives(t) = evaluate_objective(beliefs,A,n);

while t <= maxIter && (t <= 2 || abs(objectives(end-1) - objectives(end)) > TolFun)
    
    % for t = 2:numIters
    % gradient means -gradient
    % ???????
    grad = -A + modOneDEntropyDerivative(beliefs);
    map = munkres(grad');
    
    if isDebug
        debugLog(t) = var2struct(grad, map);
        disp('grad');
        grad
        disp('map');
        map
    end
    
%     assert(length(map) == n && all(map > 0));
    objectives(t) = evaluate_objective(beliefs,A,n);
    beliefs = update_beliefs_deterministic_step(beliefs,map,n,t,isDebug);
    t = t + 1;

    if t > 2 && isDebug
        fprintf('[ .M t=%d] objectives(end-1) = %g; objectives(end) = %g; diff = %g, TolFun = %g\n', t, objectives(end-1), objectives(end), abs(objectives(end-1) - objectives(end)), TolFun);
    end
end

objectives(t) = evaluate_objective(beliefs,A,n);

end



function obj = evaluate_objective(beliefs,A,n)
U = -sum(sum(beliefs.*A));
H = +sum(sum(modOneDEntropy(beliefs)));
obj = U + H;
end

function H = evaluate_entropy(beliefs,n)
H = +sum(sum(modOneDEntropy(beliefs)));
end

function ent = modOneDEntropy(beliefs)
ent = beliefs .* log(beliefs) - (1 - beliefs) .* log(1 - beliefs);
ent(beliefs == 0 | beliefs == 1) = 0;
end

function ent = modOneDEntropyDerivative(beliefs)
ent = 2 + log(beliefs) + log(1 - beliefs);
ent(beliefs == 0 | beliefs == 1) = 0;
end

%%this just uses the 2/(t+2) rule
function beliefs = update_beliefs_deterministic_step(beliefs,map,n,t,isDebug)
step = 2/(2 + t);
beliefs = (1 - step)*beliefs;

if isDebug
    disp('post de-step, pre step');
    beliefs
end

for j = 1:n
    i = map(j);
    beliefs(i,j) = beliefs(i,j) + step;
end

if isDebug
    disp('post step');
    beliefs
end

end


