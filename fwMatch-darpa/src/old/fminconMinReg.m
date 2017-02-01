function [beliefs, obj] = fminconMinReg(beliefs, X, TolFun)

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

% extra debug
fprintf('fminconMinReg called with reweight = %s\n', num2str(o.reweight'));

beliefs = o.beliefs;

a = vec(A);

    function w = hessMult(x, lambda, v)
        % hessMult  Hessian at x times v with Lagrange multipliers v
        %
        % For the Bethe objective, the Hessian is diagonal, because the
        % objective is separable (sum of B_{ij}; no terms involve different
        % Bs) and the constraints are linear (so they contribute zero).
        
%         w = 1 ./ (1 - x);
        w = (1 - 2*x) ./ (1 - x);
    end

    function H = hessFcn(x, lambda)
        n = length(x);
        B = (1 - 2*x) ./ (x .* (1 - x));
        H = spdiags(B, 0, n, n);
    end

    function [f, g] = FBethe(x)
        f = sum(-x.*a + x.*log(x) - (1 - x).*log(1 - x));
        g = -a + 2 + log(x) + log(1 - x);
    end

prob.options = optimset('Algorithm', 'interior-point', ...    
    'DerivativeCheck', 'on', 'FinDiffType', 'central', ...    
    'GradObj', 'on', ...    
    'Hessian', 'user-supplied', ...
    'HessFcn', @hessFcn, ...
    'TolFun', TolFun);

%     'HessMult', @hessMult, 'SubproblemAlgorithm', 'cg', ...    
% prob.options = optimset('Algorithm', 'interior-point', ...
%     'GradObj', 'on', ...
%     'DerivativeCheck', 'on', 'FinDiffType', 'central', ...
%     'TolFun', TolFun);

% Let the computer compute vectorized indices
N = D * D;
inds = reshape(1:N, D, D);

% Equality constraints.
% For imperfect matchings, move into inequality constraints.
% TODO: Make sparse

prob.objective = @FBethe;

badAeq = zeros(2*D, N);
badbeq = ones(2*D, 1);

prob.Aeq = zeros(2*D, N);
prob.beq = ones(2*D, 1);

for j = 1:D
    % For each column j, \sum_{ij} B_{ij} = 1    
    badAeq(j, inds(:,j)) = 1;
end

for i = 1:D
    % For each row i, \sum_{ij} B_{ij} = 1
    badAeq(D + i, inds(i,:)) = 1;
end

[prob.Aeq, prob.beq] = reduceEqualityConstraints(badAeq, badbeq);

% Bound constraints. (0, 1)
prob.lb = zeros(N, 1);
prob.ub = ones(N, 1);

prob.x0 = vec(beliefs);
prob.solver = 'fmincon';

[x, obj, exitflag] = fmincon(prob);

% x = fmincon(@FBethe, vec(beliefs), [], [], prob.Aeq, prob.beq, prob.lb, prob.ub);

beliefs = reshape(x, D, D);

end

