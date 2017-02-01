function st = structWrapStr(st)
% structWrapStr  Wrap any string values in st with cell
%
%   Purpose is to facilitate inclusion into dataset.

    
for v = fieldnames(st)'
    v1 = v{1};
    cv = getfield(st, v1);
    if isstr(cv)
        st = setfield(st, v1, {cv});
    end
end

end
