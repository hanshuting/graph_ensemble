function [N1 N2] = find_node2pair_arrays(model)

% first, just compute # neighbors
NN = zeros(model.nnodes,1);
for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    NN(i)=NN(i)+1;
    NN(j)=NN(j)+1;
end

%maxN = 4;
maxN = max(NN);
N1 = -1*ones(model.nnodes,maxN);
N2 = -1*ones(model.nnodes,maxN);
where1 = ones(model.nnodes,1);
where2 = ones(model.nnodes,1);
for c=1:model.ncliques
    i = model.pairs(c,1);
    j = model.pairs(c,2);
    N1(i,where1(i))=c;
    N2(j,where2(j))=c;
    where1(i)=where1(i)+1;
    where2(j)=where2(j)+1;
end
