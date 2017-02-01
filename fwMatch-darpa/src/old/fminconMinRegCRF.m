function [theta beliefs objs iterTaus iterOptimValues] = fminconMinRegFeats(Y, features, lambda, varargin)

p = inputParser;
p.KeepUnmatched = true; % Enable extraneous options.

p.addRequired('Y', @isnumeric);
p.addRequired('features', @iscell);
p.addRequired('lambda', @(x) isnumeric(x) && isscalar(x));

[M N] = size(Y);

p.addParamValue('beliefs', []);
p.addParamValue('reweight', ones(N, 1));
p.addParamValue('TolFun', 1e-10);
p.addParamValue('MaxIter', inf);
p.addParamValue('MaxFunEvals', 10000);
p.addParamValue('debug', false);
p.addParamValue('matchConstraints', []);
p.addParamValue('boundConstraints', true);
p.addParamValue('stopIfSingular', false)

p.parse(Y, features, lambda, varargin{:});
o = p.Results;

% extra debug
disp('fminconMinReg called with reweight:');
disp(num2str(o.reweight'));

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

function w = hessMult(tt, lam, v)
    % Diagonal of the diagonal part of Hessian, in vector form.

    dtt = (1 - rsum.*tt) ./ (tt .* (1 - tt));        
    w   = 1/lambda * G*(G'*v) + dtt.*v;
end

function [f, g] = obj(tt)
    Gytt = G'*(y - tt);
    
    sq = 1/(2*lambda) * sum(Gytt.^2);
    H  = sum(tt.*log(tt) - (rsum - 1).*(1 - tt).*log(1 - tt));
    f  = sq + H;
    
    g  = -1/lambda*G*Gytt + rsum + log(tt) + (rsum - 1).*log(1 - tt);    
end


objs = [];
iterOptimValues = {};
iterTaus = {};

function stop = outFun(x, optimValues, state)
    stop = false;
    
    objs(end+1) = optimValues.fval;    
    iterOptimValues{end+1} = optimValues;
    iterTaus{end+1} = x;    
end

prob.options = optimset('Algorithm', 'interior-point', ...
    'GradObj', 'on', ...
    'Hessian', 'user-supplied', ...
    'SubproblemAlgorithm', 'cg', ...
    'HessMult', @hessMult, ...
    'OutputFcn', @outFun, ...
    'TolFun', o.TolFun, ...
    'MaxIter', o.MaxIter, ...
    'MaxFunEvals', o.MaxFunEvals, ...
    'Display', 'iter', ...
    'TolX', 1e-14);

% prob.options = optimset('Algorithm', 'interior-point', ...
%     'GradObj', 'on', ...
%     'Hessian', 'user-supplied', ...
%     'SubproblemAlgorithm', 'cg', ...
%     'HessFcn', @hessFcn, ...    
%     'TolFun', o.TolFun, ...
%     'Display', 'iter', ...
%     'TolX', 1e-14);

% Let the computer compute vectorized indices.
% Same constraints for each sample, so really a block-diagonal.

if isempty(o.matchConstraints)
    o.matchConstraints = makeMatchConstraintsBrute(N);
end

allAeqs = repmat({o.matchConstraints.Aeq}, M, 1);

prob.Aeq = blkdiag(allAeqs{:});
prob.beq = repmat(o.matchConstraints.beq, M, 1);

prob.objective = @obj;

% Apparently they need not be necessary (given feasible start)
if o.boundConstraints
    % Bound constraints [0, 1]. Though 1 is redundant given equality
    % constraints.
    prob.lb = zeros(M * N^2, 1);
end

prob.x0 = tau;
prob.solver = 'fmincon';

if o.stopIfSingular
    singularS       = warning('error', 'MATLAB:singularMatrix');
    nearlySingularS = warning('error', 'MATLAB:nearlySingularMatrix');
end

try
    [tauFinal, objVal, exitflag, output] = fmincon(prob);
catch err
    if strcmp(err.identifier, 'MATLAB:singularMatrix') || ...
       strcmp(err.identifier, 'MATLAB:nearlySingularMatrix')
   
        % Restore warning states
        warning(singularS, 'MATLAB:singularMatrix');
        warning(nearlySingularS, 'MATLAB:nearlySingularMatrix');
   
        warning('fminconMinRegFeats', 'Singular matrix detected in fmincon; bailing out and returning partial solutions in objs, iterTaus, and iterOptimValues.');
        theta = [];
        beliefs = [];
        
        iterOptimValues = cell2mat(iterOptimValues);
        return;    
    else
        rethrow(err);
    end
end

iterOptimValues = cell2mat(iterOptimValues);

theta = 1 / lambda * G' * (y - tauFinal);

beliefs = reshape(tauFinal, N, N, M);

% For sanity checking
sumGY    = computeFeatInnerProduct(features, reshape(y, N, N, M));
minSumGT = computeFeatInnerProduct(features, beliefs);
thetaCheck    = sum(sumGY - minSumGT, 1)' / lambda;

assertElementsAlmostEqual(theta, thetaCheck);

fprintf('fmincon took %d iterations and %d function evaluations', ...
    output.iterations, output.funcCount);

output

end

