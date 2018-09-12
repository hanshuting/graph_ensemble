function N = bipcover(M)
k = 2;
n = size(M);
N = zeros(n*k);

comp = 0;

while comp ~= 1

    for r = 1:n
        for c = 1:n
            if(r == c)
                N((r-1)*k+1:r*k, (c-1)*k+1:c*k) = M(r,c)*eye(k);
            else
                %Create a random permutation and weight it by the r,c position in M
                P = [0 1; 1 0];

                N((r-1)*k+1:r*k, (c-1)*k+1:c*k) = M(r,c)*P;
            end
        end
    end

    %Needed for symmetric matrices?
    for r=1:n*k
        for c=1:r
            N(r,c) = N(c,r);
        end
    end

    G = (N - diag(diag(N))) ~= 0;
    [comp, C] = graphconncomp(sparse(G));
    %overload connectivity requirement
    comp = 1;
end