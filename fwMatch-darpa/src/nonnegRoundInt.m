function [y, scale] = nonnegRoundInt(x, intTy)
%nonnegRoundInt Scale and round nonnegative input to integral range
%
%   [y, scale, shift] = nonnegRoundInt(x, intTy) shifts and scales x to use
%   the full range of intTy ('int32', 'int64', etc.) using the formula
%   y = scale * x.

    lo = min(x(:));
    hi = max(x(:));
    
    assert(lo >= 0, 'x must be nonnegative.');
    
    scale = double(intmax(intTy) - 1) / hi;
    y = scale * x;

end

