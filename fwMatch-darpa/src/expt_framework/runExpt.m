function runExpt(name, id)
% runExpt  Launch experiment at a particular configuration.

    % TODO: Check dependencies of run and warn about uncommitted files.
    % TODO: Automatic *branching* if run with uncommitted files.
    % TODO: Data versioning too.
    %        - Idea: Search for all params with "file" string (falsely
    %          detects output files)
    %        - Annotate (not idiotproof).
    %
    %        - Goal: Make idiotproof. (Usable after all-nighters w/o making
    %          mistake.)
    % TODO: Save the git configuration a
    
    t = timetic();
    tic(t);
    exptDir = sprintf('expt/%s', name);

    diary off;

    addpath(exptDir);

    % Ingest the params struct by evaluating the config file.
    % Ensure evaluating correct config file by CD'ing to its folder
    rootDir = pwd;
    cd(sprintf('%s/',exptDir))
    config = sprintf('config%d', id);
    feval(config);
    cd(sprintf('%s/',rootDir))

    % Automatically log everything that comes onscreen henceforth.
    %diary(sprintf('%s/diary%d.txt', exptDir, id));

    % Inject variables into params.
    params.exptDir        = exptDir;
    params.saveFileName   = sprintf('result%d.mat', id);
    params.saveFile       = sprintf('%s/result%d.mat', exptDir, id);
    params.checkpointFile = sprintf('%s/%s_cp.mat',      exptDir, config);
    params.logFile        = sprintf('%s/%s_log.txt',     exptDir, config);

    % Run!
    % eval('run(params)');
    run(params);

    % Clean up path, in case we were called interactively.
    rmpath(exptDir);
    toc(t)
end

