function reportHelper(exptName, report)
%    p = inputParser;
%    p.addParamValue('saveMat', []);
%    p.addParamValue('saveTxt', [])
%    p.parse(varargin{:});

%    o = p.Result;

% TODO: Report file should include provenance/versioning information.
%       Should error if result files have different provenances?
% TODO: Integrate report writing with Makefiles.

    existOpt = @(x) isfield(report, x) && ~isempty(getfield(report, x));

    for id = report.ids
        % Config is in struct params. 
        run(sprintf('expt/%s/config%d.m', exptName, id));

        % Result is in struct result. 
        try
            result = load(sprintf('expt/%s/result%d.mat', exptName, id));
        catch err
            fprintf('Error:\n');
            err
            continue;
        end

        % Shortcut: Handle empty configVars.
        if ~existOpt('configVars')
            c = params;
        else
            c = structFields(params, report.configVars);
        end

        % User *must* specify resultVars. This is because we almost always
        % will have non-scalar outputs in results, which we won't include.
        r = structFields(result, report.resultVars);

        % Computed fields. The fields of report.computeVar give the name of
        % each computed variables, and the values of report.computeVar give
        % the functions to compute, which take a single argument of result.
        %
        % TODO: Consider computed variables that also take params.
        cr = struct();
        if existOpt('computeVars')
            if ~isstruct(report.computeVars)
                error('makeExptReportHelper:computeVars', 'computeVars must be a struct. Format is field = variable; value = lambda = value of struct.')
            end

            for k  = fieldnames(report.computeVars)'
                k1 = k{1};
                v  = getfield(report.computeVars, k1);
                cr = setfield(cr, k1, v(result)); 
            end
        end

        % Wrap any strings we have as cell arrays, in order to concatenate
        % into a dataset.
        c  = structWrapStr(c);
        r  = structWrapStr(r);
        cr = structWrapStr(cr);

        rows{id} = horzcat(dataset(id), struct2dataset(c), struct2dataset(r), struct2dataset(cr));
    end

    rows = rows.';
    format long;
    summary = vertcat(rows{:})
    if existOpt('saveFile')
        export(summary, 'file', report.saveFile);
    end
end

