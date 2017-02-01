function [L b_ij b_i dtheta_ij dtheta_i Z_ij] = meanfield_bprop(model,theta_ij,theta_i,...
    maxiter,convthresh,loss,dorec)

%[b_ij b_i A] = meanfield(model,psi_ij,psi_i,maxiter)

% you give a convergence threshold and a maximum number of iterations.
% forward propagation stops when EITHER of these is reached.
% backward propagation always uses the same number of iterations as
% forward.

% generally, one would use something like
% dorec=0, maxiter=inf, convthresh = 1e-5 -> Back belief propagation
% dorec=1, maxiter=5,   convthresh = 0    -> Truncated fitting
% but you can use both if you are a wierdo.  (In practice, one probably
% doesn't want to allow infinite iterations with BBP.)

% loss is a function object that will be called like so:
% [L db_i db_ij] = loss(b_i, b_ij);

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_ij should have size [model.nvals^2 size(model.pairs,1)]');
end

if nargin < 7
    dorec = 1; % means do record messages
end

if dorec==1 && maxiter > 50
    warning('dorec=%d with maxiter=%d',dorec,maxiter);
end

psi_i  = exp(theta_i);
psi_ij = exp(theta_ij);

% initial beliefs
b_i  = 1+0*psi_i;
b_ij = 1+0*psi_ij;

% all messages are stored here before being updated
if dorec
    wpred = maxiter*(2*model.nnodes-1)*model.nvals;
    mstor = zeros(wpred,1); w = 0;
else
    w     = [];
    mstor = [];
end

b_i0 = ones(model.nvals,1);

for reps=1:maxiter    
    b_i_old = b_i;
    
    for j=[1:model.nnodes model.nnodes-1:-1:1]
        if dorec, for yj=1:model.nvals, w=w+1; mstor(w) = b_i(yj,j); end; end
        
        for yj=1:model.nvals
            b_i0(yj) = psi_i(yj,j);
            for d=model.N1(j,:)
                if d==-1, continue; end
                i = model.pairs(d,2);
                for yi=1:model.nvals
                    index = yj + (yi-1)*model.nvals;
                    %b_i(yj,j) = b_i(yj,j)*psi_ij(index,d)^b_i(yi,i);
                    b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
                end
            end
            for d=model.N2(j,:)
                if d==-1, continue; end
                i = model.pairs(d,1);
                for yi=1:model.nvals
                    index = yi + (yj-1)*model.nvals;
                    %b_i(yj,j) = b_i(yj,j)*psi_ij(index,d)^b_i(yi,i);
                    b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
                end
            end
        end
        %b_i(:,j) = b_i0(:,j) / sum( b_i0(:,j) );
        b_i(:,j) = norm(b_i0);
    end
    
    conv = max(abs(b_i_old(:)-b_i(:)));
    if conv < convthresh
        % do backprop by same number of iterations
        maxiter = reps;
        break
    end
end

for c=1:model.ncliques
    for yi=1:model.nvals
        for yj=1:model.nvals
            i = model.pairs(c,1);
            j = model.pairs(c,2);
            index = yi + (yj-1)*model.nvals;
            b_ij(index,c) = b_i(yi,i)*b_i(yj,j);
        end
    end
end


% compute loss
[L db_i db_ij] = loss(b_i, b_ij);
% propagate back to messages

for c=1:model.ncliques
    for yi=1:model.nvals
        for yj=1:model.nvals
            i = model.pairs(c,1);
            j = model.pairs(c,2);
            index = yi + (yj-1)*model.nvals;
            db_i(yi,i) = db_i(yi,i) + db_ij(index,c)*b_ij(index,c)/b_i(yi,i);
            db_i(yj,j) = db_i(yj,j) + db_ij(index,c)*b_ij(index,c)/b_i(yj,j);
        end
    end
end

dpsi_ij = 0*psi_ij;
dpsi_i  = 0*psi_i;

for reps=1:maxiter    
    %for j=model.nnodes:-1:1
    for j=fliplr([1:model.nnodes model.nnodes-1:-1:1])
        % recover non-normalized marginals
        for yj=1:model.nvals
            b_i0(yj) = psi_i(yj,j);
            for d=model.N1(j,:)
                if d==-1, continue; end
                i = model.pairs(d,2);
                for yi=1:model.nvals
                    index = yj + (yi-1)*model.nvals;
                    b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
                end
            end
            for d=model.N2(j,:)
                if d==-1, continue; end
                i = model.pairs(d,1);
                for yi=1:model.nvals
                    index = yi + (yj-1)*model.nvals;
                    b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
                end
            end
        end
        % denormed derivatives
        db_i0 = denorm(db_i(:,j),b_i0);
        for yj=1:model.nvals
            for d=model.N1(j,:)
                if d==-1, continue; end
                i = model.pairs(d,2);
                for yi=1:model.nvals
                    index = yj + (yi-1)*model.nvals;
                    %b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
                    db_i(yi,i) = db_i(yi,i) + db_i0(yj)*b_i0(yj)*log(psi_ij(index,d));
                    dpsi_ij(index,d) = dpsi_ij(index,d) + db_i0(yj)*b_i0(yj)/psi_ij(index,d)*b_i(yi,i);
                end
            end
            for d=model.N2(j,:)
                if d==-1, continue; end
                i = model.pairs(d,1);
                for yi=1:model.nvals
                    index = yi + (yj-1)*model.nvals;
                    %b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
                    db_i(yi,i) = db_i(yi,i) + db_i0(yj)*b_i0(yj)*log(psi_ij(index,d));
                    dpsi_ij(index,d) = dpsi_ij(index,d) + db_i0(yj)*b_i0(yj)/psi_ij(index,d)*b_i(yi,i);
                end
            end
        end
        
        for yj=1:model.nvals
            dpsi_i(yj,j) = dpsi_i(yj,j) + db_i0(yj)*b_i0(yj)/psi_i(yj,j);
        end
        
        if dorec, for yj=model.nvals:-1:1, b_i(yj,j) = mstor(w); w=w-1; end; end
        db_i(:,j) = 0;
    end
end

dtheta_i  = dpsi_i .*psi_i;
dtheta_ij = dpsi_ij.*psi_ij;

end


    function b = norm(a)
        b = a/sum(a);
    end

    function d2 = denorm(d,a)
        s  = sum(a);
        d2 = d/s - d'*a/s^2;
    end