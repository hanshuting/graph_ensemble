function [b_ij b_i A Z_ij] = trw(model,theta_ij,theta_i,rho,maxiter,damp,convthresh)

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_ij should have size [model.nvals^2 size(model.pairs,1)]');
end

% initialize messages
n  = ones(model.nvals,model.nnodes  );
m1 = ones(model.nvals,model.ncliques); % messages to first variable in clique
m2 = ones(model.nvals,model.ncliques); % messages to second variable

if nargin < 6 || isempty(damp)
    damp = 0;
end
if nargin < 7
    %convthresh = .000002; % ~ 0.001 accuracy
    convthresh = .00002; % ~ 0.01 accuracy
end

if damp~=0
    error('damping not implemented in regular trw (use fast).');
end

b_i = n;

psi_ij = exp(theta_ij);
psi_i  = exp(theta_i);

psi_ij = psi_ij.^(1/rho);

S  = zeros(model.nvals,model.nvals);
m0 = zeros(model.nvals,1);

for reps=1:maxiter
    tic
    for c0=1:2*model.ncliques
        
        if c0 <= model.ncliques
            c    = c0;
            mode = 1;
        else
            c = 2*model.ncliques-c0+1;
            mode = 2;
        end
        
        i = model.pairs(c,1);
        j = model.pairs(c,2);
        % compute n(y_i) and n(y_j)
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
                %m2(yj,c)=0;
                m0(yj) = 0;
                for yi=1:model.nvals
                    index = yi + (yj-1)*model.nvals; % correct?
                    %m2(yj,c) = m2(yj,c) + psi_ij(index,c)*psi_i(yi,i)*n(yi,i)/m1(yi,c);
                    S(yi,yj) = psi_ij(index,c)*psi_i(yi,i)*n(yi,i)/m1(yi,c);
                    m0(yj) = m0(yj) + S(yi,yj);
                end
            end
            %m2(:,c) = m2(:,c) / sum(m2(:,c));
            k = sum(S(:));
            m2(:,c) = m0 / k;
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
                %m1(yi,c)=0;
                m0(yi) = 0;
                for yj=1:model.nvals
                    index = yi + (yj-1)*model.nvals; % correct?
                    %m1(yi,c) = m1(yi,c) + psi_ij(index,c)*psi_i(yj,j)*n(yj,j)/m2(yj,c);
                    S(yi,yj) = psi_ij(index,c)*psi_i(yj,j)*n(yj,j)/m2(yj,c);
                    m0(yi) = m0(yi) + S(yi,yj);
                end
            end
            %m1(:,c) = m1(:,c) / sum(m1(:,c));
            k = sum(S(:));
            m1(:,c) = m0 / k;
        end
    end
    
    b_i_old = b_i;
    b_i = n.*psi_i;
    b_i = b_i ./ repmat(sum(b_i,1),model.nvals,1);
    
    b_ij = psi_ij;
    for c=1:model.ncliques
        i = model.pairs(c,1);
        j = model.pairs(c,2);
        for yi=1:model.nvals
            for yj=1:model.nvals
                index = yi + (yj-1)*model.nvals;
                b_ij(index,c) = b_ij(index,c)*psi_i(yi,i)*psi_i(yj,j)*n(yi,i)*n(yj,j)/m1(yi,c)/m2(yj,c);
            end
        end
    end
    Z_ij = sum(b_ij,1);
    b_ij = b_ij ./ repmat(Z_ij,model.nvals^2,1);
    

    %b_iSAVE{reps}  = b_i;
    
    conv = max(abs(b_i_old(:)-b_i(:)));
    if conv <= convthresh
        break
    end
end

%keyboard

% calculate approx log-partition function
who_i  = b_i >0;
who_ij = b_ij>0;
% count number of times each variable occurs in a clique
noccur = sum(model.N1~=-1,2)+sum(model.N2~=-1,2);
noccur = repmat(noccur',model.nvals,1);
%
A = sum(theta_ij(who_ij).*b_ij(who_ij)) + ...
    sum(theta_i (who_i ).*b_i (who_i )) - ...
    sum((1-rho*noccur(who_i)).*b_i(who_i).*log(b_i(who_i))) - ...
    rho*sum(b_ij(who_ij).*log(b_ij(who_ij)));