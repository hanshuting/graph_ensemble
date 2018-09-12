function [errs, aucs] = evalRoommateErrorMatrix(thetas, datafile)
% evalRoommates  Compute error of each parameter at each sample.

    % This takes a while
    data = load(datafile);
    M = length(data.Ys);
    S = length(thetas);

    errs = zeros(M, S);
    aucs = zeros(M, S);

    for m = 1:M
        for s = 1:S
            [errs(m,s), aucs(m,s)] = evalThetaUnipartite(thetas{s}, data.features(m,:), data.Ys(m));
        end
    end

end
