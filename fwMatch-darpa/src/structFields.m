function s = structFields(st, fields)
% structFields  Return structure with only specified fields.
    s = struct();
    for f = fields
        f1 = f{1};
        s  = setfield(s, f1, getfield(st, f1));
    end
end
