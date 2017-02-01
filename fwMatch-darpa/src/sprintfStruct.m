function str = sprintfStruct(st)

    fields = fieldnames(st);
    
    parts = cellfun(@(x) sprintf('%s = %g, ', x, getfield(st, x)), fields, 'UniformOutput', false);

    % Chop off the ', '
    parts{end}(end) = [];
    parts{end}(end) = sprintf('\n');
    
    str = [parts{:}];

end

