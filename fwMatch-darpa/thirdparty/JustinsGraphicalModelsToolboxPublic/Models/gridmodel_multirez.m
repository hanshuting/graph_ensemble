function model = gridmodel_multirez(ly,lx,nvals,nrez)

nnodes   = 0;
ncliques = 0;
for rez=1:nrez
    ly2 = ceil(ly*.5^(rez-1));
    lx2 = ceil(lx*.5^(rez-1));
    N{rez} = reshape(nnodes+(1:ly2*lx2),ly2,lx2);
    nnodes   = nnodes   + ly2*lx2;
    ncliques = ncliques + ly2*(lx2-1) + (ly2-1)*lx2 + double(rez<nrez)*ly2*lx2;
end

model.nnodes    = nnodes;
%model.nodetype  = ones(nnodes,1);
%model.nodewhere = 1:nnodes;

model.ncliques = ncliques;
model.nvals    = nvals;
model.pairtype = zeros(ncliques,1);

model.pairs     = zeros(ncliques,2);
model.modeltype = 'grid_multirez';
model.treenum   = zeros(ncliques,1);
model.ly = ly;
model.lx = lx;

%model.tree2clique = zeros(max(ly-1,lx-1),lx+ly);

%tree2clique_vert = zeros(ly-1,lx);
%tree2clique_hor  = zeros(lx-1,ly);

where=1;
for rez=1:nrez
    ly2 = ceil(ly*.5^(rez-1));
    lx2 = ceil(lx*.5^(rez-1));
    for x=1:lx2
        % first, vertical connections
        for y=1:ly2-1
            i = N{rez}(y  ,x);
            j = N{rez}(y+1,x);
            model.pairs(where,:)  = [i j];
            model.pairtype(where) = 1;
            %model.treenum(where)  = 1;
            where = where+1;
        end
        
        % now, horizontal
        if x < lx2
            for y=1:ly2
                i = N{rez}(y,x  );
                j = N{rez}(y,x+1);
                model.pairs(where,:)  = [i j];
                model.pairtype(where) = 2;
                %model.treenum(where)  = 2;
                where = where+1;
            end
        end
        
        % now, inter-scale
        if rez<nrez
            for y=1:ly2
                x_next = ceil(x/2);
                y_next = ceil(y/2);
                i = N{rez  }(y     ,x     );
                j = N{rez+1}(y_next,x_next);
                model.pairs(where,:)  = [i j];
                model.pairtype(where) = 2+rez;
                %model.treenum(where)  = 3;
                where = where+1;
            end
        end
    end
end

[model.N1 model.N2] = find_node2pair_arrays(model);