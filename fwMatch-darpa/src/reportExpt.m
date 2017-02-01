function reportExpt(name)
    exptDir = ['expt/' name];
    reportConfigs = dir([exptDir '/report*.m']);

    for n = 1:length(reportConfigs)
        fullPath = [exptDir '/' reportConfigs(n).name];
        run(fullPath);
        report
        reportHelper(name, report); 
    end
end

