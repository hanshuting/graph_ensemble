function [obj, dTheta] = bpLikeFeats(theta, X, features, FXBar)
% betheLikeFeats  Bethe (pseudo)-likelihood with features

    M = size(features,1);
    K = size(features,2);    
    
    obj = 0;
    beliefs = {};    
    %M = 2; %todo: remove
 
        tic;
        max(theta(:))
    for m = 1:M
        % Linear A
         K2 = K;
        [Dm, Dm2] = size(features{m});
        assert(Dm == Dm2 && K == K2, 'Dimensions of features{m} inconsistent');
        
        B = zeros(Dm, Dm);            

        for k = 1:K
            B = B + theta(k) .* features{m,k};
        end            
        
        % BP code cannot handle zero-weight edges. (Really? This doesn't
        % seem to make sense since we put zero-weight edges on the
        % off-diagonal anyways.)
%        B = B + eps;
        
        % TODO: Handle unipartite case...
%        B = B - min(B(:));
       % fprintf('Max B is %g\n', max(B(:)));
              A = B; 

       %A = biadjacencyToAdjacency(B);
       aa = exp(A);
       assert(~any(isinf(aa(:))));
 %      [logZbp2, bpMsgs2, bpMuMin2] = bpPerfectMatching(max(0.0000000001,exp(A)));
        [logZbp,iters,mytime,bpMuMin] = estper(max(0.0000000000001,exp(A)),'~/FrankenWolf/matlab/bp_perm/estper',10^(-10));
        
       assert(all(abs(sum(bpMuMin,1) -1) < 0.001));
       assert(all(abs(sum(bpMuMin,2) -1) < 0.001));

        %bpMuMin = adjacencyToBiadjacency(bpMuMin);
      
        beliefs{m} = bpMuMin;     
        assert(~any(isnan(bpMuMin(:))));
        if iscell(X)
            Ym = X{m};
        else
            Ym = expandPerm(X(m,:));
        end                        
        
        for k = 1:K
            Fm = features{m,k};
            obj = obj - theta(k) * sum(Fm(Ym));
        end

        % negative log-like => + logZ => - FBethe
        assert(~isnan(logZbp));
%        assert(~isinf(logZbp));
        obj = obj + logZbp;

        fprintf('+');
    end
    toc
    fprintf('\n');

    % Subgradient step. (Batch gradients)

    dTheta = zeros(K, 1);
    for k = 1:K            
        for m = 1:M
            assert(all(~isnan(beliefs{m}(:))));
            assert(all(~isnan(features{m,k}(:))));

            dTheta(k) = dTheta(k) - FXBar(m,k) + sum(vec(features{m,k}.*beliefs{m}));
        end            
    end
    fprintf('finished iter\n');

end


