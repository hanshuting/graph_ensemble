function [L b_ij b_i dtheta_ij dtheta_i Z_ij] = trw_bprop_scheduled(model,theta_ij,theta_i,rho,...
    maxiter,convthresh,loss,dorec)

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
    error('theta_i should have size [model.nvals^2 size(model.pairs,1)]');
end

if nargin < 8
    dorec = 1; % means do record messages
end

if dorec==1 && maxiter > 50
    warning('dorec=%d with maxiter=%d',dorec,maxiter);
end

% initialize messages
n  = ones(model.nvals,model.nnodes  );
m1 = ones(model.nvals,model.ncliques); % messages to first variable in clique
m2 = ones(model.nvals,model.ncliques); % messages to second variable

tree_ncliques = sum(double(model.tree2clique>0),1);
ntree = length(tree_ncliques);

% all messages are stored here before being updated
if dorec
    %wpred = 2*maxiter*model.ncliques*model.nvals;
    %mstor = zeros(wpred,1);
    %w = 0;
    
    % same as before but create an array
    wpred = 2*maxiter*tree_ncliques*model.nvals;
    mstor = zeros(max(wpred),ntree);
    w     = zeros(ntree,1);
else
    w     = [];
    mstor = [];
end

b_i = n;

psi_i   = exp(theta_i);
psi_ij0 = exp(theta_ij);

psi_ij = psi_ij0.^(1/rho);


S  = zeros(model.nvals,model.nvals);
m0 = zeros(model.nvals,1);

for reps=1:maxiter
    for block=1:size(model.treeschedule,2)
        for treenum=1:size(model.treeschedule,1)
            tree = model.treeschedule(treenum,block);
            if tree==0, continue; end
            
            ncliques = tree_ncliques(tree);
            for c0=1:2*ncliques;
                if c0 <= ncliques
                    c = model.tree2clique(c0,tree);
                    mode = 1;
                else
                    c = model.tree2clique(2*ncliques-c0+1,tree);
                    mode = 2;
                end
                
                i = model.pairs(c,1);
                j = model.pairs(c,2);
                
                if mode==1
                    % compute n(y_i)
                    for yi=1:model.nvals
                        n(yi,i) = 1;
                        for d=model.N1(i,:)
                            if d==-1, continue; end
                            n(yi,i) = n(yi,i)*m1(yi,d)^rho;
                        end
                        if isnan(n(yi,i))
                            keyboard
                        end
                        for d=model.N2(i,:)
                            if d==-1, continue; end
                            n(yi,i) = n(yi,i)*m2(yi,d)^rho;
                        end
                        if isnan(n(yi,i))
                            keyboard
                        end
                    end
                    % compute m(y_j)
                    for yj=1:model.nvals
                        m0(yj) = 0;
                        for yi=1:model.nvals
                            index = yi + (yj-1)*model.nvals; % correct?
                            S(yi,yj) = psi_ij(index,c)*psi_i(yi,i)*n(yi,i)/m1(yi,c);
                            m0(yj) = m0(yj) + S(yi,yj);
                            if isnan(sum(S(:)))
                                keyboard
                            end
                        end
                    end
                    if dorec, for yj=1:model.nvals, w(tree)=w(tree)+1; mstor(w(tree),tree) = m2(yj,c); end; end
                    k = sum(S(:));
                    m2(:,c) = m0 / k;
                    
                    if isnan(sum(m2(:)))
                        keyboard
                    end
                else
                    % compute n(y_j)
                    for yj=1:model.nvals
                        n(yj,j) = 1;
                        for d=model.N1(j,:)
                            if d==-1, continue; end
                            n(yj,j) = n(yj,j)*m1(yj,d)^rho;
                        end
                        for d=model.N2(j,:)
                            if d==-1, continue; end
                            n(yj,j) = n(yj,j)*m2(yj,d)^rho;
                        end
                    end
                    % compute m(y_i)
                    for yi=1:model.nvals
                        m0(yi) = 0;
                        for yj=1:model.nvals
                            index = yi + (yj-1)*model.nvals; % correct?
                            S(yi,yj) = psi_ij(index,c)*psi_i(yj,j)*n(yj,j)/m2(yj,c);
                            m0(yi) = m0(yi) + S(yi,yj);
                        end
                    end
                    if dorec, for yi=1:model.nvals, w(tree)=w(tree)+1; mstor(w(tree),tree) = m1(yi,c); end; end
                    k = sum(S(:));
                    m1(:,c) = m0 / k;
                end
                
                if isnan(sum(n(:)))
                    keyboard
                end
            end
        end
    end
    
    b_i_old = b_i;
    b_i = n.*psi_i;
    b_i = b_i ./ repmat(sum(b_i,1),model.nvals,1);
    
    conv = max(abs(b_i_old(:)-b_i(:)));
    if conv < convthresh
        % do backprop by same number of iterations
        maxiter = reps;
        break
    end
end

% recompute all inwards messages
for i=1:size(b_i,2)
    for yi=1:model.nvals
        n(yi,i) = 1;
        for d=model.N1(i,:)
            if d==-1, continue; end
            n(yi,i) = n(yi,i)*m1(yi,d)^rho;
        end
        for d=model.N2(i,:)
            if d==-1, continue; end
            n(yi,i) = n(yi,i)*m2(yi,d)^rho;
        end
    end
end

b_i0 = n.*psi_i;
Z_i  = sum(b_i0,1);
Z_i  = repmat(Z_i,model.nvals,1);
b_i  = b_i0 ./ Z_i;

b_ij0 = psi_ij;
for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    for yi=1:model.nvals
        for yj=1:model.nvals
            index = yi + (yj-1)*model.nvals;
            b_ij0(index,c) = b_ij0(index,c)*psi_i(yi,i)*psi_i(yj,j)*n(yi,i)*n(yj,j)/m1(yi,c)/m2(yj,c);
        end
    end
end
Z_ij = sum(b_ij0,1);
Z_ij = repmat(Z_ij,model.nvals^2,1);
b_ij = b_ij0 ./ Z_ij;

% compute loss
[L db_i db_ij] = loss(b_i, b_ij);
% propagate back to messages
dm1  = 0*m1;
dm2  = 0*m2;

% get derivs w.r.t. unnormalized beliefs
db_i0  = db_i  ./ Z_i  - repmat(sum(db_i .*b_i0 ./Z_i .^2,1),model.nvals  ,1);
db_ij0 = db_ij ./ Z_ij - repmat(sum(db_ij.*b_ij0./Z_ij.^2,1),model.nvals^2,1);

% starter propagations (onto params and n)
dn      = db_i0 .* psi_i;
dpsi_i  = db_i0 .*n;
dpsi_ij = db_ij0.*b_ij0./psi_ij;

for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    for yi=1:model.nvals
        for yj=1:model.nvals
            index = yi + (yj-1)*model.nvals;
            %b_ij0(index,c) = b_ij0(index,c)*psi_i(yi,i)*psi_i(yj,j)*n(yi,i)*n(yj,j)/m1(yi,c)/m2(yj,c);
            dpsi_i(yi,i) = dpsi_i(yi,i) + db_ij0(index,c)*b_ij0(index,c)/psi_i(yi,i);
            dpsi_i(yj,j) = dpsi_i(yj,j) + db_ij0(index,c)*b_ij0(index,c)/psi_i(yj,j);
            dn(yi,i)     = dn(yi,i)     + db_ij0(index,c)*b_ij0(index,c)/n(yi,i);
            dn(yj,j)     = dn(yj,j)     + db_ij0(index,c)*b_ij0(index,c)/n(yj,j);
            dm1(yi,c)    = dm1(yi,c)    - db_ij0(index,c)*b_ij0(index,c)/m1(yi,c);
            dm2(yj,c)    = dm2(yj,c)    - db_ij0(index,c)*b_ij0(index,c)/m2(yj,c);
        end
    end
end

% push from n back to m1 and m2
for i=1:size(b_i,2)
    for yi=1:model.nvals
        %n(yi,i) = 1;
        for d=model.N1(i,:)
            if d==-1, continue; end
            %n(yi,i) = n(yi,i)*m1(yi,d)^rho;
            dm1(yi,d) = dm1(yi,d) + rho * dn(yi,i)*n(yi,i)/m1(yi,d);
        end
        for d=model.N2(i,:)
            if d==-1, continue; end
            %n(yi,i) = n(yi,i)*m2(yi,d)^rho;
            dm2(yi,d) = dm2(yi,d) + rho * dn(yi,i)*n(yi,i)/m2(yi,d);
        end
    end
end

if isnan(sum(dm1(:)) + sum(dm2(:)))
    keyboard
end

for reps=1:maxiter
    for block=size(model.treeschedule,2):-1:1
        for treenum=1:size(model.treeschedule,1) % reverse order, but parallel!
            tree = model.treeschedule(treenum,block);
            if tree==0, continue; end
            
            ncliques = tree_ncliques(tree);
            for c0=2*ncliques:-1:1;
                if c0 <= ncliques
                    c = model.tree2clique(c0,tree);
                    mode = 1;
                else
                    c = model.tree2clique(2*ncliques-c0+1,tree);
                    mode = 2;
                end
                
                i = model.pairs(c,1);
                j = model.pairs(c,2);
                
                if mode==1
                    % compute n(y_i)
                    for yi=1:model.nvals
                        n(yi,i) = 1;
                        for d=model.N1(i,:)
                            if d==-1, continue; end
                            n(yi,i) = n(yi,i)*m1(yi,d)^rho;
                        end
                        for d=model.N2(i,:)
                            if d==-1, continue; end
                            n(yi,i) = n(yi,i)*m2(yi,d)^rho;
                        end
                    end
                    % compute m(y_j)
                    for yj=1:model.nvals
                        m0(yj) = 0;
                        for yi=1:model.nvals
                            index = yi + (yj-1)*model.nvals; % correct?
                            S(yi,yj) = psi_ij(index,c)*psi_i(yi,i)*n(yi,i)/m1(yi,c);
                            Snon(yi,yj) = psi_ij(index,c)*psi_i(yi,i)/m1(yi,c); % NEW
                            m0(yj) = m0(yj) + S(yi,yj);
                        end
                    end
                    k = sum(S(:));
                    
                    % propagate them derivs
                    dm0 = dm2(:,c)/k - sum(dm2(:,c).*m0)/k^2;
                    dS  = repmat(dm0',model.nvals,1);
                    dn  = 0*n(:,i);
                    for yj=1:model.nvals
                        for yi=1:model.nvals
                            index = yi + (yj-1)*model.nvals; % correct?
                            %S(yi,yj) = psi_ij(index,c)*psi_i(yi,i)*n(yi,i)/m1(yi,c);
                            dpsi_ij(index,c) = dpsi_ij(index,c) + dm0(yj)*S(yi,yj)/psi_ij(index,c);
                            dpsi_i(yi,i)     = dpsi_i(yi,i)     + dm0(yj)*S(yi,yj)/psi_i(yi,i);
                            %dn(yi)           = dn(yi)           + dm0(yj)*S(yi,yj)/n(yi,i);
                            dn(yi)           = dn(yi)           + dm0(yj)*Snon(yi,yj); % NEW
                            dm1(yi,c)        = dm1(yi,c)        - dm0(yj)*S(yi,yj)/m1(yi,c);
                            if isnan(dn(yi))
                                keyboard
                            end
                        end
                    end
                    for yi=1:model.nvals
                        for d=model.N1(i,:)
                            if d==-1, continue; end
                            dm1(yi,d) = dm1(yi,d) + rho * dn(yi)*n(yi,i)/m1(yi,d);
                            if isnan(dm1(yi,d))
                                keyboard
                            end
                        end
                        for d=model.N2(i,:)
                            if d==-1, continue; end
                            dm2(yi,d) = dm2(yi,d) + rho * dn(yi)*n(yi,i)/m2(yi,d);
                            if isnan(dm2(yi,d))
                                keyboard
                            end
                        end
                    end
                    
                    dm2(:,c) = 0;
                    if dorec, for yj=model.nvals:-1:1, m2(yj,c) = mstor(w(tree),tree); w(tree)=w(tree)-1; end; end
                else
                    % compute n(y_j)
                    for yj=1:model.nvals
                        n(yj,j) = 1;
                        for d=model.N1(j,:)
                            if d==-1, continue; end
                            n(yj,j) = n(yj,j)*m1(yj,d)^rho;
                        end
                        for d=model.N2(j,:)
                            if d==-1, continue; end
                            n(yj,j) = n(yj,j)*m2(yj,d)^rho;
                        end
                    end
                    % compute m(y_i)
                    for yi=1:model.nvals
                        m0(yi) = 0;
                        for yj=1:model.nvals
                            index = yi + (yj-1)*model.nvals; % correct?
                            S(yi,yj)     = psi_ij(index,c)*psi_i(yj,j)*n(yj,j)/m2(yj,c);
                            Snon(yi,yj) = psi_ij(index,c)*psi_i(yj,j)/m2(yj,c); % NEW
                            m0(yi) = m0(yi) + S(yi,yj);
                        end
                    end
                    k = sum(S(:));
                    
                    % propagate them derivs
                    dm0 = dm1(:,c)/k - sum(dm1(:,c).*m0)/k^2;
                    dS  = repmat(dm0,1,model.nvals);
                    dn  = 0*n(:,j);
                    for yj=1:model.nvals
                        for yi=1:model.nvals
                            index = yi + (yj-1)*model.nvals; % correct?
                            %S(yi,yj) = psi_ij(index,c)*psi_i(yj,j)*n(yj,j)/m2(yj,c);
                            dpsi_ij(index,c) = dpsi_ij(index,c) + dm0(yi)*S(yi,yj)/psi_ij(index,c);
                            dpsi_i(yj,j)     = dpsi_i(yj,j)     + dm0(yi)*S(yi,yj)/psi_i(yj,j);
                            %dn(yj)           = dn(yj)           + dm0(yi)*S(yi,yj)/n(yj,j);
                            dn(yj)           = dn(yj)           + dm0(yi)*Snon(yi,yj); % NEW
                            dm2(yj,c)        = dm2(yj,c)        - dm0(yi)*S(yi,yj)/m2(yj,c);
                            
                            if isnan(dn(yj))
                                keyboard
                            end
                        end
                    end
                    
                    for yj=1:model.nvals
                        for d=model.N1(j,:)
                            if d==-1, continue; end
                            dm1(yj,d) = dm1(yj,d) + dn(yj)*n(yj,j) * rho / m1(yj,d);
                            if isnan(dm1(yj,d))
                                keyboard
                            end
                        end
                        for d=model.N2(j,:)
                            if d==-1, continue; end
                            dm2(yj,d) = dm2(yj,d) + dn(yj)*n(yj,j) * rho / m2(yj,d);
                            if isnan(dm2(yj,d))
                                keyboard
                            end
                        end
                    end
                    
                    dm1(:,c) = 0;
                    if dorec, for yi=model.nvals:-1:1, m1(yi,c) = mstor(w(tree),tree); w(tree)=w(tree)-1; end; end
                end
            end
        end
    end    
    if isnan(sum(m1(:))+sum(m2(:))+sum(dpsi_i(:))+sum(dpsi_ij(:)))
        keyboard
    end
end

% psi_ij = psi_ij.^(1/rho);
dpsi_ij = (1/rho)*psi_ij.^(1-rho).*dpsi_ij;

dtheta_i  = dpsi_i .*psi_i;
dtheta_ij = dpsi_ij.*psi_ij0;


% % test if db_ij0 is correct (yup)
% L1 = L;
% e = 1e-5;
% b_ij0(5) = b_ij0(5) + e;
% Z_ij = sum(b_ij0,1);
% Z_ij = repmat(Z_ij,model.nvals^2,1);
% b_ij = b_ij0 ./ Z_ij;
% L2     = 0*sum(b_i(:).^2) + sum(b_ij(:).^2);
% (1/e)*(L2-L1)
% db_ij0(5)
% keyboard