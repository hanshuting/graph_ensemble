function [b_ij b_i A] = meanfield(model,theta_ij,theta_i,maxiter, convthresh)

if sum(size(theta_i) ~= [model.nvals model.nnodes])
    error('theta_i should have size [model.nvals model.nnodes]');
end
if sum(size(theta_ij) ~= [model.nvals^2 size(model.pairs,1)])
    error('theta_i should have size [model.nvals^2 size(model.pairs,1)]');
end

psi_ij = exp(theta_ij);
psi_i  = exp(theta_i);

b_i = ones(model.nvals,model.nnodes);

for iter=1:maxiter
    b_iOLD = b_i;
    
    for j=[1:model.nnodes model.nnodes-1:-1:1]
       for yj=1:model.nvals
           b_i(yj,j) = psi_i(yj,j);
           for d=model.N1(j,:)
                if d==-1, continue; end
                i = model.pairs(d,2);
                for yi=1:model.nvals
                    index = yj + (yi-1)*model.nvals;
                    b_i(yj,j) = b_i(yj,j)*psi_ij(index,d)^b_i(yi,i);
                end
            end
            for d=model.N2(j,:)
                if d==-1, continue; end
                i = model.pairs(d,1);
                for yi=1:model.nvals
                    index = yi + (yj-1)*model.nvals;
                    b_i(yj,j) = b_i(yj,j)*psi_ij(index,d)^b_i(yi,i);
                end
            end
            %b_i2 = b_i2*
       end
       b_i(:,j) = b_i(:,j) / sum( b_i(:,j) );
    end
    
    conv = max(abs( b_i(:) - b_iOLD(:) ));
    %fprintf('iter: %d  conv: %f \n', iter, conv);
    if conv < convthresh
        iter
        break;
    end
end

b_ij = zeros(model.nvals^2,model.ncliques);
for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    for yi=1:model.nvals
        for yj=1:model.nvals
            index = yi + (yj-1)*model.nvals;
            b_ij(index,c) = b_i(yi,i)*b_i(yj,j);
        end
    end
end

A = sum(b_ij(:).*log(psi_ij(:))) + sum(b_i(:).*log(psi_i(:))) - ...
    sum( b_i(:) .* log( b_i(:) ));


% code below brute-force attempts to further reduce the KL-divergece.
% (which should be impossible!)
% options = optimset('Display','iter','GradObj','off','LargeScale','off');
% logb_i2 = fminunc(@eval,log(b_i),options);
% [~, b_i] = eval(logb_i2);

function [KL b_i] = eval(logb_i)
    % first, normalize
    b_i = exp(logb_i);
    for i=1:model.nnodes
        b_i(:,i) = b_i(:,i) / sum(b_i(:,i));
    end
    
    KL = 0;
    for i=1:model.nnodes
        for yi=1:model.nvals
            KL = KL + b_i(yi,i)*log(b_i(yi,i));
        end
    end
    for i=1:model.nnodes
        for yi=1:model.nvals
            KL = KL - log(psi_i(yi,i))*b_i(yi,i);
        end
    end
    
    for c=1:model.ncliques
        i = model.pairs(c,1);
        j = model.pairs(c,2);
        for yi=1:model.nvals
            for yj=1:model.nvals
                index = yi + (yj-1)*model.nvals;
                KL = KL - log(psi_ij(index,c))*b_i(yi,i)*b_i(yj,j);
            end
        end
    end
end

end

