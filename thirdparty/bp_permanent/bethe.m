function [F, U, H] = bethe(G);

F = 0;

U = 0;
H = 0;

for i=1:G.N
    %for j=i+1:G.N
    for j=i+1:G.N
        if (G.A(i,j))
            nonzero = find(G.b2{i}{j}(:)~=0);
            
            F = F + sum(G.b2{i}{j}(nonzero).*log(G.b2{i}{j}(nonzero)./G.phi2{i}{j}(nonzero)));

            U = U-sum(G.b2{i}{j}(nonzero).*log(G.phi2{i}{j}(nonzero)));
            H = H-sum(G.b2{i}{j}(nonzero).*log(G.b2{i}{j}(nonzero)));
%            disp(['added H(' num2str(i) ',' num2str(j) ')']);
        end
    end
    
    nonzero = find(G.b(i,:)~=0);
    
    F = F+(G.k-1)*sum(G.b(i,nonzero).*log(G.phi(i,nonzero)./G.b(i,nonzero)));
    
    U = U+(G.k-1)*sum(G.b(i,nonzero).*log(G.phi(i,nonzero)));
    H = H+(G.k-1)*sum(G.b(i,nonzero).*log(G.b(i,nonzero)));
%    disp(['subtracted H(' num2str(i) ')']);
end
