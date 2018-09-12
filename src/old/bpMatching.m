function [logZ, beliefs, msgs] = bpMatching(G, varargin)

[n, m] = size(G);
p = inputParser;

%G(i,j) should be the weight on edge (i,j) in the graph defined by the
%non-zero entries of G

%'reweight' is an optional symmetric matrix specifying the reweighting parameters
defaultC = ones(n,1);
addOptional(p,'reweight',defaultC,@isnumeric);

%'msgs' is an optional initialization of the msgs
defmitoj = rand(n, n, 2);
for i = 1:n
    for j = 1:n
        if G(i,j) == 0
            defmitoj(i,j,1) = 0;
            defmitoj(i,j,2) = 0;
        end
    end
end
addOptional(p,'msgs',defmitoj,@isnumeric);

%control damping parameter
addOptional(p,'damp',.99,@isnumeric);

%control convergence tolerance
addParamValue(p,'TolMsg',1e-12,@isnumeric);

parse(p,varargin{:});

C = p.Results.reweight;
mitoj = p.Results.msgs;
damp = p.Results.damp;
TolMsg = p.Results.TolMsg;

conv = 0;

%used to test for convergence
temp = mitoj(1:n, 1:n, 1:2);

%iterate until messages converge
while conv == 0

    for i = 1:n

        %compute the product of all incoming msgs to i evaluated at 0
        %this speed-up only works if all msgs are not -Inf
        %this is NEVER a problem for connected graphs
        intoi0 = sum( (G(:,i) ~= 0).*(C(:).*mitoj(:,i,1) + (C(i)-1)*mitoj(i, :, 1)')  );

        %compute the sum of the product of all incoming msgs to i with
        %exactly one edge equal to one        
        nbrs = find(G(:,i) ~= 0);
        v = (intoi0 + C(nbrs).*(mitoj(nbrs, i, 2)-mitoj(nbrs, i, 1)) + (C(i)-1)*(mitoj(i, nbrs, 2) - mitoj(i, nbrs, 1))' + G(i,nbrs)');
        lsev = logSumExp(v);
        
        %compute the msgs from i to j
        for j = find(G(i,:))
            %remove j from the sum of zero messages
            diff = intoi0 - (C(j)*mitoj(j, i, 1) + (C(i)-1)*mitoj(i, j, 1));
            %extra term in v that belongs to j
            extra = intoi0 + C(j)*(mitoj(j, i, 2) - mitoj(j, i, 1)) + (C(i)-1)*(mitoj(i, j, 2) - mitoj(j, i, 1)) + G(i,j);

            if size(find((G(:,i) ~= 0))) < 2
                %This happens when the vertex is a leaf
                zero = 0;
                one = 0;
            else
                one = diff;
                zero = - (C(j)*mitoj(j, i, 1) + (C(i)-1)*mitoj(i, j, 1)) + real(logSumExp([intoi0; lsev; pi*complex(0,1) + extra]));
            end

            %normalize the messages
            nm = logSumExp([zero,  one]);

            zero = zero - nm;
            one = one - nm;

            mitoj(i, j, 1) = (1-damp)*mitoj(i, j, 1) + (damp)*zero;
            mitoj(i, j, 2) = (1-damp)*mitoj(i, j, 2) + (damp)*one;
        end        
    end

    %check for convergence  
    if max(max(max(abs(temp - mitoj(1:n, 1:n, 1:2))))) < TolMsg
        conv = 1;
    else
        temp = mitoj(1:n, 1:n, 1:2);
    end
end

%calculate beliefs at convergence
msgs = mitoj(1:n, 1:n, 1:2);
logZ = 0;


%edge beliefs
b = zeros(n,n,2);
for i = 1:n
    for j = i+1:n
        if G(i,j) ~= 0
            b(i,j,1) = mitoj(i,j,1) + mitoj(j,i,1);
            b(i,j,2) = mitoj(i,j,2) + mitoj(j,i,2) + G(i,j);
            nm = logSumExp([b(i,j,1) ,  b(i,j,2)]);
            b(i,j,1) = exp(b(i,j,1) - nm);
            b(i,j,2) = exp(b(i,j,2) - nm);
            b(j,i,1) = b(i,j,1);
            b(j,i,2) = b(i,j,2);

            logZ = logZ + G(i,j)*b(i,j,2) + (C(i) + C(j) - 1)*xlogx(b(i,j,1)) - xlogx(b(i,j,2));
        end
    end
end


%vertex contributions
for i = 1:n
    s = 0; 
    for j = find(G(i,:))
        s = s + b(i,j,2);
    end
    logZ = logZ - C(i)*xlogx(1-s);
end

beliefs = b(:,:,2);

