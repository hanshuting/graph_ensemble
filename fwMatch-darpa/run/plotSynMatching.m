%% Load data
% Change this to something fuller (or adaptive reg) as it comes in.

% THIS ONE IS THE HIGH SIGNAL PROBLEM.
st = load('expt/synMatching/config3_cp.mat');

%% Plot optimal curve
% TODO: Make less manual.
% Ms = [10 30 100 300 1000 3000 10000];
Ms   = [10 30 100 300 1000 3000];
Mixs = 1:length(Ms);

ylims = [-16, -4];
xlims = [10 3000];

rwTrainRwTest     = st.meanObjs(Mixs,1);
betheTrainBetheTest = st.meanObjs(Mixs,2);

paperFig(0.5, 0.2);
semilogx(Ms, st.meanExactObjs(Mixs), '-o', ...
     Ms, st.meanObjs(Mixs,1), '-x', ...
     Ms, st.meanObjs(Mixs,2), '-v');
% axis tight;
% title('Objective values at individual minimizers');

xlabel('Samples');
ylabel('Objectives');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);

print -dpdf fig/syn_matching_individual_opt.pdf;

%% Plot cross-curves
% Compute cross-curves
rwTrainBetheTest = st.meanTestObjs(Mixs,1,2);
rwTrainExactTest = st.meanExactTestObjs(Mixs,1);

betheTrainRwTest   = st.meanTestObjs(Mixs,2,1);
betheTrainExactTest = st.meanExactTestObjs(Mixs,2);

paperFig(0.5, 0.2);
semilogx(Ms, rwTrainExactTest, '-o', ...
         Ms, rwTrainRwTest, '-x', ...
         Ms, rwTrainBetheTest, '-v');
     
% title('Objective values at TRW minimizer');
xlabel('Samples');
ylabel('Objectives');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);

print -dpdf fig/syn_matching_rw_opt.pdf;

%% other one
paperFig(0.5, 0.2);
semilogx(Ms, betheTrainExactTest, '-o', ...
         Ms, betheTrainRwTest, '-x', ...
         Ms, betheTrainBetheTest, '-v');

% title('Objective values at Bethe minimizer');
xlabel('Samples');
ylabel('Objectives');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);
legend({'Exact', 'RW (\rho = 0.5)', 'Bethe (\rho = 1.0)'}, 'Box', 'off', 'Location', 'best');

print -dpdf fig/syn_matching_bethe_opt.pdf;

%% True likelihood comparisons
paperFig(0.5, 0.2);
semilogx(Ms, st.meanExactObjs(Mixs), '-o', ...
         Ms, rwTrainExactTest, '-x', ...
         Ms, betheTrainExactTest, '-v');
xlabel('Samples');
ylabel('Reg. Like.');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);

% legend({'Exact', 'RW (\rho = 0.5)', 'Bethe (\rho = 1.0)'}, 'Box', 'off', 'Location', 'best');
% legend({'Exact', 'TRW (\rho = 0.5)', 'Bethe (\rho = 1.0)'});

print -dpdf fig/syn_matching_exact_likes.pdf;
