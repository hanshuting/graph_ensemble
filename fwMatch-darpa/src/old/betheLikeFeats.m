function [obj, dTheta, beliefs] = betheLikeFeats(theta, features, FXBar, varargin)
        
% betheLikeFeats  Bethe (pseudo)-likelihood with features
    [M, K] = size(features);
    D = size(features{1,1}, 1);

    p = inputParser;
    p.KeepUnmatched = true; % Enable extraneous options.
    p.addRequired('theta');
    p.addRequired('features', @iscell);    
    p.addRequired('FXBar', @isnumeric);
    
    % Initial values
    p.addParamValue('beliefs', []);    
    p.addParamValue('rho', ones(M, 1));
    p.addParamValue('X', []);    


    % TODO: FACTOR OUT
    fwTolFun = 1e-6;

    FBethe = zeros(M, 1);
    obj = 0;
    
    % No warm start
    if isempty(o.beliefs)
        o.beliefs = 1/D * ones(D, D, M);
    end

    beliefs = zeros(D, D, M);
    beliefEntropyRatios = zeros(M, 1);
    fwNIters = zeros(M, 1);
    lapRatios = zeros(M, 1);
    mexTimes = zeros(M, 1);
                
    unifEntropy = beliefEntropy(1/D * ones(D));
    
    for m = 1:M
        % Linear A
        Dm = size(features{m,1}, 1);        
        
        A = zeros(Dm, Dm);            

        for k = 1:K
            A = A - theta(k) .* features{m,k};
        end
        
        tic;
%        [beliefs(:,:,m), FBetheHist, lapTimes] = fwBipartite_mex_debug(o.beliefs(:,:,m), A, o.rho, fwTolFun);        
        [beliefs(:,:,m), FBetheHist, lapTimes] = fwBipartite_mex(o.beliefs(:,:,m), A, o.rho, fwTolFun);        
%         [beliefs(:,:,m), FBetheHist] = fwBipartite_mex_old(o.beliefs(:,:,m), A, o.rho, fwTolFun);        
        lapTimes = 0;
        mexTimes(m) = toc;

        FBethe(m)          = FBetheHist(end);
        beliefEntropyRatios(m) = beliefEntropy(beliefs(:,:,m)) / unifEntropy;
        fwNIters(m)        = length(FBetheHist);
        lapRatios(m)       = sum(lapTimes) / mexTimes(m);

        if iscell(o.X)
            Ym = o.X{m};
        else
            Ym = expandPerm(o.X(m,:));
        end

        for k = 1:K
            Fm = features{m,k};
            obj = obj - theta(k) * sum(Fm(Ym));
        end

        % negative log-like => + logZ => - FBethe
        obj = obj - FBethe(m);

        fprintf('+');
    end
    fprintf('\n');

    % Subgradient step. (Batch gradients)

    if nargin >= 2
        dTheta = zeros(K, 1);
        for k = 1:K            
            for m = 1:M
                dTheta(k) = dTheta(k) - FXBar(m,k) + sum(vec(features{m,k} .* beliefs(:,:,m)));
            end            
        end
    end
    
    % Debug diagnostics
    fwTimePerIter = mexTimes ./ fwNIters;
    diagnostics = dataset(fwNIters, mexTimes, fwTimePerIter, lapRatios, beliefEntropyRatios, FBethe);
    summary(diagnostics)

end

