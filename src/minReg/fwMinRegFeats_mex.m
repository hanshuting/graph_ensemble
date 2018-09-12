function [theta, beliefs, objs] = fwMinRegFeats_mex(Y, features, lambda, varargin)
p = inputParser;
p.KeepUnmatched = true; % Enable extraneous options.

p.addRequired('Y', @isnumeric);
p.addRequired('features', @iscell);
p.addRequired('lambda', @(x) isnumeric(Y) && isscalar(x));

[M N] = size(Y);

p.addParamValue('beliefs', []);
p.addParamValue('reweight', ones(N, 1));
p.addParamValue('TolGap', 1e-6);
p.addParamValue('maxIter', double(intmax('int32')));
p.addParamValue('debug', false);

p.parse(Y, features, lambda, varargin{:});
o = p.Results;

% Construct the big G matrix. Since MATLAB matrices are column-major, this
% makes G'*tau computations efficient in the loop.
G = cell2mat(cellfun(@vec, features, 'UniformOutput', false));

% Inefficient but whatevs
Ymats = cell(M, 1);
for m = 1:M
    Ymats{m} = expandPerm(Y(m,:), @zeros);    
end

y = vertcat(cell2mat(cellfun(@vec, Ymats, 'UniformOutput', false)));
Gty = G'*y;

if o.debug
    disp('Debug on!\n');
    o.reweight
    [theta, tau, objs] = fwMinRegFeats_mex_private_debug(N, Gty, G, lambda, o.reweight, o);
else
    [theta, tau, objs] = fwMinRegFeats_mex_private(      N, Gty, G, lambda, o.reweight, o);    
end

beliefs = reshape(tau, N, N, M);

end
