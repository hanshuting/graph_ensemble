function [perm numRequiredIters objectives final_entropy] = fw_permanent_simple(A,numIters)

[n m] = size(A);
assert(n == m,'input must be square');


%initialization: uniform
beliefs = (1/n)*ones(n,n);


%%We seek min_{beliefs} -sum_ij beliefs_ij A_ij  + sum_ij beliefs_ij log(beliefs_ij) - sum_ij (1 - beliefs_ij) log(1 - beliefs_ij)

objectives = zeros(1,numIters);
t=1;
objectives(t) = evaluate_objective(beliefs,A,n); 
for t = 2:numIters
   gradient = -A + modOneDEntropyDerivative(beliefs);
   [map_assignment value] = hungarian(gradient);
   assert(length(map_assignment) == n && all(map_assignment > 0));
   
   beliefs = update_beliefs_deterministic_step(beliefs,map_assignment,n,t);

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
    ent = arrayfun(@(x) x*log(x) - (1-x)*log(1-x),beliefs);
end
function ent = modOneDEntropyDerivative(beliefs)
    ent = arrayfun(@(x) 2 + log(x*(1-x)), beliefs);
end


function assertCalibrated(beliefs,n)
    rowSums = sum(beliefs,2);
    colSums = sum(beliefs,1);
    for i = 1:n
        assertEqual(rowSums(i),1.0);
        assertEqual(colSums(i),1.0);
    end    
end



%%this just uses the 1/(t+2) rule
function beliefs = update_beliefs_deterministic_step(beliefs,map_assignment,n,t)
    step = 2/(2 + t);
    beliefs = update_beliefs_with_step(beliefs,map_assignment,step,n);
end

function beliefs = update_beliefs_with_step(beliefs,map_assignment,step,n)
    beliefs  = (1 - step)*beliefs;
    for j = 1:n
        i = map_assignment(j);
        beliefs(i,j) = beliefs(i,j) + step;%%todo: check this
    end
    
end


function assertEqual(a,b)
    msg = [num2str(a) ' and ' num2str(b)];
    assert(closeToZero(a - b),msg);
end

function b = closeToZero(a) 
    b = abs(a) < 0.0001;
end


