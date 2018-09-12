function G=bp(G);
%iterate bp algorithm

damping = .5;

oldM = G.m;

for i=1:length(G.R)
    source = G.R(i);
    dest = G.C(i);
    
    mprod = G.phi(source,:)';
    tmp = find(G.A(:,source));
    for j = 1:length(tmp)
        if (tmp(j)~=dest)
            %disp(['adding message from ' num2str(j) ' to ' num2str(source)]);
            mprod = mprod.*oldM{tmp(j)}{source};
        end
    end

    for val=1:G.k
        G.m{source}{dest}(val) = sum(G.psi{source}{dest}(:,val).*mprod); 
    end
    
    G.m{source}{dest} = G.m{source}{dest}/sum(G.m{source}{dest}(:));
    
    %damping
    G.m{source}{dest} = exp(safeplus(safelog(oldM{source}{dest}),damping*(safeplus(safelog(G.m{source}{dest}),-safelog(oldM{source}{dest})))));
    
    %G.m{source}{dest} = damping*oldM{source}{dest}+(1-damping)*G.m{source}{dest};
end

for i=1:G.N
    mprod = ones(G.k,1);
    tmp = find(G.A(:,i));
    for j = 1:length(tmp)
        mprod = mprod.*G.m{tmp(j)}{i};
    end

    G.b(i,:) = (G.phi(i,:).*mprod');
    
    G.b(i,:) = G.b(i,:)/sum(G.b(i,:));
end

%compute pairwise beliefs

%if (pairwise)
    for i=1:length(G.R)
        source = G.R(i);
        dest = G.C(i);

        mprod = ones(G.k,1);
        tmp = find(G.A(:,source));
        for j = 1:length(tmp)
            if (tmp(j)~=dest)
                %disp(['adding message from ' num2str(j) ' to ' num2str(source)]);
                mprod = mprod.*oldM{tmp(j)}{source};
            end
        end

        G.b2{source}{dest} = G.phi2{source}{dest}.*repmat(mprod,1,G.k);

        mprod = ones(G.k,1);
        tmp = find(G.A(:,dest));
        for j = 1:length(tmp)
            if (tmp(j)~=source)
                %disp(['adding message from ' num2str(j) ' to ' num2str(source)]);
                mprod = mprod.*oldM{tmp(j)}{dest};
            end
        end
        G.b2{source}{dest} = G.b2{source}{dest}.*repmat(mprod',G.k,1);

        G.b2{source}{dest} = G.b2{source}{dest}/sum(G.b2{source}{dest}(:));
    end
%end


G.iter = G.iter+1;