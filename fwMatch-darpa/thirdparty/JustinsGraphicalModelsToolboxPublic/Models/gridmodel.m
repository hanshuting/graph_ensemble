function model = gridmodel(ly,lx,nvals)

N = reshape(1:ly*lx,ly,lx);

model.nnodes    = ly*lx;
model.nodetype  = ones(ly*lx,1);
model.nodewhere = 1:ly*lx;

model.ncliques = ly*(lx-1) + (ly-1)*lx;
model.nvals    = nvals;
model.pairtype = zeros(model.ncliques,1);

model.pairs     = zeros(model.ncliques,2);
model.modeltype = 'grid';
model.treenum   = zeros(model.ncliques,1);
model.ly = ly;
model.lx = lx;

%model.tree2clique = zeros(max(ly-1,lx-1),lx+ly);

tree2clique_vert = zeros(ly-1,lx);
tree2clique_hor  = zeros(lx-1,ly);

where=1;
for x=1:lx
    % first, vertical connections
    for y=1:ly-1
        i = N(y  ,x);
        j = N(y+1,x);
        model.pairs(where,:)  = [i j];
        model.pairtype(where) = 1;
        model.treenum(where)  = 1;
        tree2clique_vert(y,x) = where;
        %model.tree2clique(y,x) = where;
        where = where+1;
    end
    
    % now, horizontal
    if x < lx
        for y=1:ly
            i = N(y,x  );
            j = N(y,x+1);
            model.pairs(where,:)  = [i j];
            model.pairtype(where) = 2;
            model.treenum(where)  = 2;
            %model.tree2clique(x,lx+y) = where;
            tree2clique_hor(x,y) = where;
            where = where+1;
        end
    end
end

%padx = max(0,ly-lx);
%pady = max(0,lx-ly);
%model.treeschedule = [[(1:lx)'; zeros(padx,1)]  [lx+(1:ly)'; zeros(pady,1)]];

%model.tree2clique = [tree2clique_vert tree2clique_hor];

% [ly lx] = size(tree2clique_vert);
% tree2clique_vert = reshape(tree2clique_vert,ly*25,lx/25);
% [ly lx] = size(tree2clique_hor);
% tree2clique_hor  = reshape(tree2clique_hor ,ly*25,lx/25);

model.tree2clique = padcat_hor(tree2clique_vert, tree2clique_hor);
s_vert = size(tree2clique_vert,2);
s_hor  = size(tree2clique_hor, 2);
model.treeschedule = padcat_hor((1:s_vert)',s_vert+(1:s_hor)');

%model.pairs

% compute neighborhoods
% model.N1 = cell(model.nnodes,1);
% model.N2 = cell(model.nnodes,1);
% for i=1:model.nnodes
%     model.N1{i} = [];
%     model.N2{i} = [];
% end
% for c=1:model.ncliques
%     i = model.pairs(c,1);
%     j = model.pairs(c,2);
%     model.N1{i} = [model.N1{i} c];
%     model.N2{j} = [model.N2{j} c];
% end

[model.N1 model.N2] = find_node2pair_arrays(model);