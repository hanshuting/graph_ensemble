function p = unpackstruct(pvec,p);

% pvec = packstruct(p)
% turns some structure arrays of parameters into a simple vector
% p2 = unpackstruct(pvec2,p);
% takes a simple vector and turns it into a structure of parameters
% with the same form as p
%
% these functions are a must have if you want to do nonlinear
% optimization in almost any nontrivial case, where it is
% hard to specify the parameters in one vector.
%
% author: justin domke  contact: (author's last name)@cs.umd.edu

fields = sort(fieldnames(p));
for i=1:length(fields)
    f = fields{i};
    if iscell(p.(f))
        for j=1:length(p.(f))
            siz = size(p.(f){j});
            len = prod(siz);
            p.(f){j} = reshape(pvec(1:len),siz);
            pvec = pvec(len+1:end);
        end
    else
        siz = size(p.(f));
        len = prod(siz);
        p.(f) = reshape(pvec(1:len),siz);
        pvec = pvec(len+1:end);
    end
end