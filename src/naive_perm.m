function [p] = naive_perm(A)

    [D, D2] = size(A);
    assert(D == D2, 'A must be square');

    Ys = perms(1:D);    
    M = size(Ys, 1);

    p = 0;
    for m = 1:M
        q = 1;
        for i = 1:D
            q = q * A(i, Ys(m,i));
        end
        p = p + q;
    end

end

