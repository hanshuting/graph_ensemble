function mx = minSoFar(x)

    m  = Inf;
    T  = length(x);
    mx = zeros(size(x));
    
    for t = 1:T
        if x(t) < m
            m = x(t);
        end
        mx(t) = m;
    end

end
