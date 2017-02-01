%% Load data
% Change this to something fuller (or adaptive reg) as it comes in.

% THIS ONE IS THE HIGH SIGNAL PROBLEM.
st = load('expt/synMatching/config30_cp.mat');

%% Plot optimal curve
% TODO: Make less manual.
% Ms = [10 30 100 300 1000 3000 10000];
Ms   = [10 30 100 300 1000 3000];
Mixs = 1:length(Ms);

ylims = [-23, -11.5];
xlims = [10 3000];

trwTrainTrwTest     = st.meanObjs(Mixs,1);
betheTrainBetheTest = st.meanObjs(Mixs,2);

paperFig(0.5, 0.2);
semilogx(Ms, st.meanExactObjs(Mixs), '-o', ...
     Ms, st.meanObjs(Mixs,1), '-x', ...
     Ms, st.meanObjs(Mixs,2), '-v');
% axis tight;
% title('Objective values at individual minimizers');
legend({'Exact', 'RW (\rho = 0.5)', 'Bethe (\rho = 1.0)'}, 'Box', 'off', 'Location', 'best');

xlabel('Samples');
ylabel('Objectives');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);

print -dpdf fig/synlo_matching_individual_opt.pdf;

%% Plot cross-curves
% Compute cross-curves
trwTrainBetheTest = st.meanTestObjs(Mixs,1,2);
trwTrainExactTest = st.meanExactTestObjs(Mixs,1);

betheTrainTrwTest   = st.meanTestObjs(Mixs,2,1);
betheTrainExactTest = st.meanExactTestObjs(Mixs,2);

paperFig(0.5, 0.2);
semilogx(Ms, trwTrainExactTest, '-o', ...
         Ms, trwTrainTrwTest, '-x', ...
         Ms, trwTrainBetheTest, '-v');
     
% title('Objective values at TRW minimizer');
xlabel('Samples');
ylabel('Objectives');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);

print -dpdf fig/synlo_matching_trw_opt.pdf;

paperFig(0.5, 0.2);
semilogx(Ms, betheTrainExactTest, '-o', ...
         Ms, betheTrainTrwTest, '-x', ...
         Ms, betheTrainBetheTest, '-v');

% title('Objective values at Bethe minimizer');
xlabel('Samples');
ylabel('Objectives');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);
legend({'Exact', 'RW (\rho = 0.5)', 'Bethe (\rho = 1.0)'}, 'Box', 'off', 'Location', 'best');

print -dpdf fig/synlo_matching_bethe_opt.pdf;

%% True likelihood comparisons
paperFig(0.5, 0.2);
semilogx(Ms, st.meanExactObjs(Mixs), '-o', ...
         Ms, trwTrainExactTest, '-x', ...
         Ms, betheTrainExactTest, '-v');
xlabel('Samples');
ylabel('Reg. Like.');

ylim(ylims);
axis([xlims ylims]);
set(gca, 'XTick', Ms, 'XTickLabel', Ms);

legend({'Exact', 'RW (\rho = 0.5)', 'Bethe (\rho = 1.0)'}, 'Box', 'off', 'Location', 'best');
% legend({'Exact', 'TRW (\rho = 0.5)', 'Bethe (\rho = 1.0)'});

print -dpdf fig/synlo_matching_exact_likes.pdf;
