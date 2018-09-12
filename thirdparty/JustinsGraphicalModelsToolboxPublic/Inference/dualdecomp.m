function [b_ijOUT b_iOUT A nreps] = dualdecomp(model,theta_ij,theta_i,rho)

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_i should have size [model.nvals^2 size(model.pairs,1)]');
end

% This function implements a RESTRICTED VERSION of dual decomposition
% We assume:
% A) Every node participates in every tree.
% B) Each pairtype has probability rho of appearing.  (so rho = 1/#trees)
% C) Every pair participates in EXACTLY ONE tree
%
% it is important that the model defines the array "treenum" which will
% be used to decompose the graph
%
% thus, in inference, the pair potentials CANNOT CHANGE, and we need to do
% optimize only over the univariate potentials  (So we don't even bother to 
% define gamma_ij)

nreps = 0;
%b_iSAVE = {};

    function [M dgamma_i] = eval(gamma_i)
        M = 0;
        dgamma_i = zeros(size(gamma_i));
        nreps = nreps + 1;
        if T==2
            b_ijOUT = zeros(size(theta_ij));
            b_iOUT  = zeros(size(theta_i));
            for t=1:T
                if t==1
                    sign = 1;
                else
                    sign = -1;
                end
                theta_i2  = theta_i/T+sign*gamma_i;
                theta_ij2 = theta_ijs{t};
                                
                maxiter = 2;
                [b_ij b_i A] = trw_fast(models{t},theta_ij2/rho,theta_i2/rho,1,maxiter,0,1e-10);
                
                % recompute entropy and such (notice we have no rho here)
                who_i  = b_i >0;
                who_ij = b_ij>0;
                % count number of times each variable occurs in a clique
                
                theta_dot_bt = sum(theta_ij2(who_ij).*b_ij(who_ij)) + ...
                               sum(theta_i2 (who_i ).*b_i (who_i ));
                
                H = -sum((1-noccur{t}(who_i)).*b_i(who_i).*log(b_i(who_i))) - ...
                     sum(b_ij(who_ij).*log(b_ij(who_ij)));
                
                M = M + theta_dot_bt + rho*H;
                
                dgamma_i = dgamma_i + sign*b_i;
                
                b_ijOUT(:,who{t}) = b_ij;
                b_iOUT  = b_iOUT + b_i/T;
                
                %b_iSAVE{nreps}  = b_iOUT;
            end
        else
            error('not implemented');
        end
    end



% make different trees
T = max(model.treenum);

if T==2
    for t=1:T
        models{t} = model;
        who{t}             = find(model.pairtype==t);
        models{t}.pairs    = model.pairs(who{t},:);
        models{t}.ncliques = size(models{t}.pairs,1);
        models{t}.treenum  = [];
        models{t}.pairtype = [];
        [models{t}.N1 models{t}.N2] = find_node2pair_arrays(models{t});
        theta_ijs{t} = theta_ij(:,who{t});
        
        noccur{t} = sum(models{t}.N1~=-1,2)+sum(models{t}.N2~=-1,2);
        noccur{t} = repmat(noccur{t}',models{t}.nvals,1);
    end
    
    gamma_i = zeros(size(theta_i));
    %options = optimset('Display','iter','GradObj','on','DerivativeCheck','off',...
    %    'LargeScale','off','TolFun',1e-10,'TolX',1e-4);
    %fminunc(@eval,gamma_i,options);
    options = [];
    options.Corr = 10;
    options.TolX = 1e-5;
    %options.TolX = 1e-4;
    %options.Method = 'newton0lbfgs';
    options.Display = 0;
        
%     options.Corr    = 100;
%     options.TolX    = 1e-50;
%     options.MaxIter = 1e5;
%     options.MaxFunEvals = 1e4;
%     options.Display = 1;
    minFunc(@eval,gamma_i,options);
else
    error('not implemented!');
    % need a more complex situation where the change is overparameterized
end

end
    

% set up parameters

% % evaluation function
% % for all trees, do inference, return M, dL/dgamma




% % calculate approx log-partition function
% who_i  = b_i >0;
% who_ij = b_ij>0;
% % count number of times each variable occurs in a clique
% noccur = sum(model.N1~=-1,2)+sum(model.N2~=-1,2);
% noccur = repmat(noccur',model.nvals,1);
% %
% A = sum(logpsi_ij(who_ij).*b_ij(who_ij)) + ...
%     sum(logpsi_i (who_i ).*b_i (who_i )) - ...
%     sum((1-rho*noccur(who_i)).*b_i(who_i).*log(b_i(who_i))) - ...
%     rho*sum(b_ij(who_ij).*log(b_ij(who_ij)));