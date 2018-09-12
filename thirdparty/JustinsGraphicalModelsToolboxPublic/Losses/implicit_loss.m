function [L b_ij b_i dpsi_ij dpsi_i] = implicit_loss(model, psi_ij, psi_i,...
    rho, maxiter, convthresh, x, lossname, dorec, mode, inference, extra_args)

% mode should either be 'bprop' or 'pert'

%if dorec==0 && (maxiter < 50 || convthresh > 1e-3)
%    warning('with dorec=0, at least maxiter>=50 and convthresh <= 1e-3 recommended')
%end

% x is a vector of observed values
% zeros correspond to hidden variables

% [L db_i db_ij] = loss(b_i, b_ij);
function [L db_i db_ij] = loss_ul(b_i, b_ij, x)
%     L = 0;
%     db_i  = 0*b_i;
%     db_ij = 0*b_ij;
%     for i=1:length(x)
%         if x(i)
%             L            = L - log(b_i(x(i),i));
%             db_i(x(i),i) = -1/b_i(x(i),i);
%         end
%     end    
    % same but no for loops
    db_i  = 0*b_i;
    db_ij = 0*b_ij;
    valid = find(x);
    who   = x(valid);
    where = sub2ind(size(b_i),who,valid);
    bs    = b_i(where);
    e     = 0*1e-5; % slop factor to deal with extreme beliefs
    L     = -sum(log(e+bs));
    db_i(where) = -1./(e+bs);
    
    %e     = 1e-5; % slop factor to deal with extreme beliefs
    %L     = L    + 1e-5*sum(log(e+b_i(:)));
    %db_i  = db_i + 1e-5*1./(e+b_i);
end

    function [L db_i db_ij] = loss_uquad(b_i, b_ij, x)
        db_i  = 0*b_i;
        db_ij = 0*b_ij;
        valid = find(x);
        who   = x(valid);
        where = sub2ind(size(b_i),who,valid);
        bs    = b_i(where);
        L     = sum(-2*(bs));
        db_i(where) = -2;
        
        bs    = b_i(:,valid);
        L     = L + sum(bs(:).^2);
        db_i(:,valid ) = db_i(:,valid) + 2*b_i(:,valid);
    end

    function [L db_i db_ij] = loss_cl(b_i, b_ij, x)
%         L = 0;
%         db_i  = 0*b_i;
%         db_ij = 0*b_ij;
%         for n=1:size(model.pairs,1)
%             i = model.pairs(n,1);
%             j = model.pairs(n,2);
%             if x(i) && x(j)
%                 who          = x(i) + (x(j)-1)*model.nvals;
%                 L            = L - log(b_ij(who,n));
%                 db_ij(who,n) = -1/b_ij(who,n);
%             end
%         end
        % same but no for loops
        validpairs = find(prod(x(model.pairs),2));
        is = model.pairs(validpairs,1);
        js = model.pairs(validpairs,2);
        who = x(is) + (x(js)-1)*model.nvals;
        wherepairs = sub2ind(size(b_ij),who,validpairs);
        bs = b_ij(wherepairs);
        L = -sum(log(bs));
        db_i  = 0*b_i;
        db_ij = 0*b_ij;
        db_ij(wherepairs) = -1./bs;
    end

    function [L db_i db_ij] = loss_cquad(b_i, b_ij, x)
        validpairs = find(prod(x(model.pairs),2));
        is = model.pairs(validpairs,1);
        js = model.pairs(validpairs,2);
        who = x(is) + (x(js)-1)*model.nvals;
        wherepairs = sub2ind(size(b_ij),who,validpairs);
        bs = b_ij(wherepairs);
        L = sum(-2*(bs));
        db_i  = 0*b_i;
        db_ij = 0*b_ij;
        db_ij(wherepairs) = -2;
        
        bs    = b_ij(:,validpairs);
        L     = L + sum(bs(:).^2);
        db_ij(:,validpairs) = db_ij(:,validpairs)  + 2*b_ij(:,validpairs);
    end

%     function [L db_i db_ij] = loss_uquad(b_i, b_ij, x)
%         db_i  = 0*b_i;
%         db_ij = 0*b_ij;
%         valid = find(x);
%         who   = x(valid);
%         where = sub2ind(size(b_i),who,valid);
%         bs    = b_i(where);
%         L     = sum(-2*(bs));
%         db_i(where) = -2;
%         
%         bs    = b_i(:,valid);
%         L     = L + sum(bs(:).^2);
%         db_i(:,valid ) = db_i(:,valid) + 2*b_i(:,valid);
%     end

    function [L db_i db_ij] = loss_acc(b_i, b_ij, x, lambd)
        L = 0;
        db_i  = 0*b_i;
        db_ij = 0*b_ij;
        sig  = @(x) 1./(1+exp(-lambd*x));
        sigp = @(x) lambd*exp(-lambd*x)/(exp(-lambd*x) + 1)^2;
        
        b_i2 = b_i;
        for i=1:length(x)
            if x(i)
                b_i2(x(i),i) = 0;
            end
        end
        [~,top] = max(b_i2,[],1);   

        for i=1:length(x)
            if x(i)
                L              = L + sig(b_i(top(i),i) - b_i(x(i),i));
                db_i(x(i)  ,i) = -sigp(b_i(top(i),i) - b_i(x(i),i));
                db_i(top(i),i) = db_i(top(i),i) + ...
                                 + sigp(b_i(top(i),i) - b_i(x(i),i));
            end
        end
    end

if strcmp(lossname,'ul')
    loss = @(b_i,b_ij) loss_ul(b_i, b_ij, x);
elseif strcmp(lossname,'cl')
    loss = @(b_i,b_ij) loss_cl(b_i, b_ij, x);
elseif strcmp(lossname,'uquad')
    loss = @(b_i,b_ij) loss_uquad(b_i, b_ij, x);
elseif strcmp(lossname,'cquad')
    loss = @(b_i,b_ij) loss_cquad(b_i, b_ij, x);
elseif length(lossname)>=3 && strcmp(lossname(1:3),'acc')
    lambd = str2double(lossname(4:end));
    loss = @(b_i,b_ij) loss_acc(b_i, b_ij, x, lambd);
else
    error('unsupported loss: %s',lossname)
end

if strcmp(mode,'bprop')
    if strcmp(inference,'trw')
        %[L b_ij b_i dpsi_ij dpsi_i] = trw_bprop(model,psi_ij,psi_i,rho,...
        %    maxiter,convthresh,loss,dorec);
        if nargout < 4
            dorec = 0; % no reason to record if aren't calculating derivatives
            [L b_ij b_i] = trw_bprop_fast(model,psi_ij,psi_i,rho,...
                maxiter,convthresh,loss,dorec);
        else
            [L b_ij b_i dpsi_ij dpsi_i] = trw_bprop_fast(model,psi_ij,psi_i,rho,...
                maxiter,convthresh,loss,dorec);
        end
    elseif strcmp(inference,'trwpll')
        [L b_ij b_i dpsi_ij dpsi_i] = trw_bprop_scheduled_fast(model,psi_ij,psi_i,rho,...
                maxiter,convthresh,loss,dorec);
    elseif strcmp(inference,'mnf')
        [L b_ij b_i dpsi_ij dpsi_i] = meanfield_bprop_fast(model,psi_ij,psi_i,...
            maxiter,convthresh,loss,dorec);
    else
        error('unsupported inference method: %s',inference);
    end
elseif strcmp(mode,'pert')
    pertmode = extra_args(1);
    emult    = extra_args(2);
    [L b_ij b_i dpsi_ij dpsi_i]   = implicit_perturbation(model,psi_ij,psi_i,rho,...
        maxiter,loss,convthresh,pertmode,inference,emult);
else
    error('mode should be bprop or pert')
end

% [L b_ij b_i dpsi_ij dpsi_i] = trw_bprop(model,psi_ij,psi_i,rho,...
%     maxiter,convthresh,loss,dorec);


end