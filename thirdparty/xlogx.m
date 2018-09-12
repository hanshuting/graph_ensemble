function y = xlogx(x)

    y = x .* log(x);
    y(x == 0) = 0;

end

