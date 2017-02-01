function test_trw_bprop_fast

randn('state',5);
rand ('state',5);

% % first, make a randomly connected graph
% model.nnodes   = 1000000;
% model.ncliques = 2000000;
% model.nvals    = 2;
% rho            = .723;
% maxiter        = 5;
% 
% model.pairs = zeros(model.ncliques,2);
% where=1;
% for k=1:model.nnodes-1
%     model.pairs(where,:) = [k k+1];
%     where=where+1;
% end
% for where=where:model.ncliques
%     k = ceil(rand*model.nnodes);
%     l = k;
%     while l==k
%         l = ceil(rand*model.nnodes);
%     end
%     model.pairs(where,:) = [k l];
%     where=where+1;
% end
% 
% %model.pairs
% 
% % % compute neighborhoods
% % model.N1 = cell(model.nnodes,1);
% % model.N2 = cell(model.nnodes,1);
% % for i=1:model.nnodes
% %     model.N1{i} = [];
% %     model.N2{i} = [];
% % end
% % for c=1:model.ncliques
% %     i = model.pairs(c,1);
% %     j = model.pairs(c,2);
% %     model.N1{i} = [model.N1{i} c];
% %     model.N2{j} = [model.N2{j} c];
% % end
% 
% maxN = 10;
% model.N1 = -1*ones(model.nnodes,maxN);
% model.N2 = -1*ones(model.nnodes,maxN);
% where1 = ones(model.nnodes,1);
% where2 = ones(model.nnodes,1);
% for c=1:model.ncliques
%     i = model.pairs(c,1);
%     j = model.pairs(c,2);
%     model.N1(i,where1(i))=c;
%     model.N2(j,where2(j))=c;
%     where1(i)=where1(i)+1;
%     where2(j)=where2(j)+1;
% end

model = gridmodel(300,300,3);
rho            = .792;
maxiter        = 5;
convthresh     = 0;

% now, make some random potentials
psi_ij = rand(model.nvals^2,model.ncliques);
psi_i  = rand(model.nvals  ,model.nnodes  );
% normalize (not necessary)
for i=1:size(psi_i,2) psi_i(:,i) = psi_i(:,i)/sum(psi_i(:,i)); end

    function [L db_i db_ij] = loss(b_i, b_ij)
        L     = sum(b_i(:).^2) + sum(b_ij(:).^2);
        db_i  = 2*b_i;
        db_ij = 2*b_ij;
    end

%tic
[L2 dpsi_ij2 dpsi_i2 b_ij2 b_i2]   = trw_bprop_fast(model,psi_ij,psi_i,rho,maxiter,@loss);
%toc

%tic
[L dpsi_ij dpsi_i b_ij b_i]   = trw_bprop(model,psi_ij,psi_i,rho,maxiter,@loss);
%toc

fprintf('differences: %f | %f | %f \n', norm(L-L2), norm(dpsi_i-dpsi_i2), norm(dpsi_ij-dpsi_ij2));

end