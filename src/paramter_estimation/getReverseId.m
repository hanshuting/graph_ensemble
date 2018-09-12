function [I,J] = getReverseId(N)
    I = [];
    J = [];
    for i=1:N-1
        I = [I;repmat(i,N-i,1)];
        J = [J;[i+1:N]'];
    end
end