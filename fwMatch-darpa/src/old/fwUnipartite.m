function [beliefs objectives] = fwUnipartite(beliefs,A,TolFun)
% fwUnipartite  Frank-Wolfe partition function for unipartite graphs
%
%   TODO: Support non-perfect matchings too (add extra entropy term)

[n m] = size(A);
assert(n == m,'input must be square');

%initialization: uniform
if isempty(beliefs)
    beliefs = (1/n)*ones(n,n);
end

if nargin == 2
    % a bit ridiculous
    TolFun = 1e-6;
end

%%We seek min_{beliefs} -sum_ij beliefs_ij A_ij  + sum_ij beliefs_ij log(beliefs_ij) - sum_ij (1 - beliefs_ij) log(1 - beliefs_ij)

%objectives = zeros(1,numIters);
t=1;
objectives(t) = evaluate_objective(beliefs,A,n); 

while t <= 2 || abs(objectives(end-1) - objectives(end)) > TolFun
    
   gradient = -A + modOneDEntropyDerivative(beliefs);
   % bmatch is a *maximization* routine.
   [M, cost] = perfectMatching(gradient);   

%    assert(~infeasible, 'bmatch returned an infeasible solution');
   
   objectives(t) = evaluate_objective(beliefs,A,n);
   beliefs = update_beliefs_deterministic_step(beliefs,M,n,t);
   t = t + 1;
end

final_objective = evaluate_objective(beliefs,A,n);
perm = exp(-final_objective);
final_entropy = evaluate_entropy(beliefs,n);

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
%     ent = arrayfun(@(x) x*log(x) - (1-x)*log(1-x),beliefs);
    ent = beliefs .* log(beliefs) - (1 - beliefs) .* log(1 - beliefs);
%     TODO:
    ent(beliefs == 0 | beliefs == 1) = 0;
end
function ent = modOneDEntropyDerivative(beliefs)
%     ent = arrayfun(@(x) 2 + log(x*(1-x)), beliefs);
    ent = 2 + log(beliefs) + log(1 - beliefs);
    ent(beliefs == 0 | beliefs == 1) = 0;
end

%%this just uses the 2/(t+2) rule
function beliefs = update_beliefs_deterministic_step(beliefs,M,n,t)
    step = 2/(2 + t);
    beliefs = update_beliefs_with_step(beliefs,M,step,n);
end

%% TODO: CHANGE THIS FOR UNIPARTITE
function beliefs = update_beliefs_with_step(beliefs,M,step,n)
    beliefs  = (1 - step)*beliefs + step*M;    
end
