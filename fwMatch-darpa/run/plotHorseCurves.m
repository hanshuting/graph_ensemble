%% Data files
% Objective traces
% mRes    = load('expt/horseReproduceCorrect/result2.mat');
j40Res  = load('expt/justinHorseCorrect/result1.mat');
j200Res = load('expt/justinHorseCorrect/result2.mat');
j20Res  = load('expt/justinHorseCorrect/result3.mat');
j10Res  = load('expt/justinHorseCorrect/result4.mat');
% Read from cp/justinHorseCache/train.mat
N_TRAIN = 200;

% Test error traces
gatherHorseTestErr;
load('cp/horseTestErr.mat');

%% Look at the beating results
beat = load('expt/horseReproduceBeatTiming/result1.mat');
beatFval = [beat.objHist.fval];
beatTime = [beat.objHist.totIterTime];


%% Extract series we want to plot

% m stands for me.
% mIter     = [mRes.objHist.iter];
% mTime     = [mRes.objHist.totIterTime];
% mTrainErr = [mRes.objHist.trainErr];
% mFval     = [mRes.objHist.fval];

j40Iter     = [j40Res.history.fval];
j40Time     = [j40Res.history.totTime];
j40Fval     = -[j40Res.history.fval] / j40Res.newScale;

j20Iter     = [j20Res.history.fval];
j20Time     = [j20Res.history.totTime];
j20Fval     = -[j20Res.history.fval] / j20Res.newScale;

j10Iter     = [j10Res.history.fval];
j10Time     = [j10Res.history.totTime];
j10Fval     = -[j10Res.history.fval] / j10Res.newScale;

% Time filter

% Don't show so much flat space on the plots.
MAX_TIME = 2e4;

findLastIx = @(times) find(times < MAX_TIME, 1, 'last');
% mIx    = 1:findLastIx(mTime);
j40Ix  = 1:findLastIx(j40Time);
% j200Ix = 1:findLastIx(j200Time);
j20Ix = 1:findLastIx(j20Time);
j10Ix = 1:findLastIx(j10Time);


beatIx = 1:findLastIx(beatTime);

%% Timing finder for 
time = beatTime(100);
j40MidIx = find(j40Time < time, 1, 'last');
j40MidTime = j40Time(j40MidIx);


%% Timing finder
lateTime    = beatTime(250);
j40LateIx   = find(j40Time < lateTime, 1, 'last');
j40LateTime = j40Time(j40LateIx);

%% Colors (see http://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html)
colors = [
    0.0000    0.4470    0.7410    0.4000 % MLE-Struct
    0.8500    0.3250    0.0980    1.0000 % MLE-Struct-wavg
    0.9290    0.6940    0.1250    1.0000 % domke40
    0.4940    0.1840    0.5560    1.0000 % domke20
    0.4660    0.6740    0.1880    0.5000 % domke10
    0.3010    0.7450    0.9330    1.0000
    0.6350    0.0780    0.1840    1.0000
    ];

lines = { '-', '-.', '--', ':', '-' };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% <IN PAPER: fig/horseObjCurve40+200.pdf>
h = paperFig(0.5, 0.25);

% Compute the primal-dual midpoint from the 200 values (to get a tighter gap)
% midPd = 0.5 * (min(mFval(mIx)) + max(j200Fval(j200Ix)));
% midPd = 0.5 * (min(mFval(mIx)) + max(j40Fval(j40Ix)));
midPd = 0.5 * (min(beatFval(beatIx)) + max([j40Fval(j40Ix) j20Fval(j20Ix)]));
% midPd = 0.5 * (min(mFval) + max(j40Fval));


% scaleFactor = (startDg - endDg) / (N_TRAIN * startDg);

% mCurve    = (mFval(mIx) - midPd)       / abs(midPd);
beatCurve = (beatFval(beatIx) - midPd)    / abs(midPd);
j40Curve  = (midPd - j40Fval(j40Ix))   / abs(midPd);
% j200Curve = (midPd - j200Fval(j200Ix)) / abs(midPd);
j20Curve = (midPd - j20Fval(j20Ix)) / abs(midPd);
j10Curve = abs(midPd - j10Fval(j10Ix)) / abs(midPd);

% % Manually specify colors to skip the mWavg line
% semilogy(beatTime(beatIx), beatCurve(beatIx), '-',  'Color', [  0    0.4470    0.7410]);
% hold on;
% semilogy(j200Time(j200Ix), j200Curve(j200Ix),  '-.',  'Color', [0.9290    0.6940    0.1250]);
% semilogy(j40Time(j40Ix), j40Curve(j40Ix),      ':', 'Color', [0.4940    0.1840    0.5560]);
% hold off;
% ylabel('log(Relative Duality Gap)');

% TODO: LINE DECORATIONS, CONSISTENT COLORS.

semilogy(beatTime(beatIx) / 60, beatCurve(beatIx),  lines{1}, 'Color', colors(1,:));
hold on;
% semilogy(j200Time(j200Ix), j200Curve(j200Ix));
semilogy(j40Time(j40Ix) / 60, j40Curve(j40Ix),      lines{3}, 'Color', colors(3,:));
semilogy(j20Time(j20Ix) / 60, j20Curve(j20Ix),      lines{4}, 'Color', colors(4,:));
semilogy(j10Time(j10Ix) / 60, j10Curve(j10Ix),      lines{5}, 'Color', colors(5,:));
hold off;
box off;
ylabel('Log Relative Duality Gap');
% legend({'MLE-Struct', 'domke40', 'domke20', 'domke10'}, 'box', 'off', 'Location', 'best')

% Linear instead of logarithmic plot
% plot(mTime(mIx), mCurve(mIx), '-',  'Color', [  0    0.4470    0.7410]);
% hold on;
% plot(j200Time(j200Ix), j200Curve(j200Ix),  '-.',  'Color', [0.9290    0.6940    0.1250]);
% plot(j40Time(j40Ix), j40Curve(j40Ix),      ':', 'Color', [0.4940    0.1840    0.5560]);
% hold off;
% ylabel('Relative Duality Gap');

% axis([-Inf Inf minVal*(1-epsilon) maxVal*(1+epsilon)]);
% ylim([0.07 0.5]);
ylim([1e-3, 1e2]);
xlim([0, 330]);

% title(sprintf('Relative difference to midpoint of primal and dual vs vs time, 200 TRW iters, final gap = %g', dg100));
xlabel('Time (min)');
     
print -dpdf fig/horseObjCurve40+200.pdf

%% <IN PAPER: figs/horseErrCurve40+200.pdf>
h = paperFig(0.5, 0.25);

% TODO: REPLACE BY BEAT...
plot(beatTime(beatIx) / 60, beatErr(beatIx),     lines{1}, 'Color', colors(1,:));
hold on;

plot(beatTime(beatIx) / 60, beatWavgErr(beatIx), lines{2},  'Color', colors(2,:));

% semilogy(j200Time(j200Ix), j200Curve(j200Ix));
plot(j40Time(j40Ix) / 60, j40Err(j40Ix),         lines{3},  'Color', colors(3,:));
plot(j20Time(j20Ix) / 60, j20Err(j20Ix),         lines{4}, 'Color', colors(4,:));
plot(j10Time(j10Ix) / 60, j10Err(j10Ix),         lines{5},  'Color', colors(5,:));
hold off;
box off;
legend({'MLE-Struct', 'MLE-Struct-wavg', 'domke40', 'domke20', 'domke10'}, 'box', 'off', 'Location', 'best');

ylim([0.08 0.45]);
xlim([0, 330]);
% title('Test error vs time, 200 TRW iters');
xlabel('Time (min)'); ylabel('Test Error');

print -dpdf fig/horseErrCurve40+200.pdf

%% Objective value plots
% h = paperFig(0.33, 0.33);
h = figure;
plot(mTime, mFval, '--k', ...
     j40Time, j40Fval, '.-b');

legend('BCFW', 'domke40');
dg40 = mFval(end) - j40Fval(end);
title(sprintf('Obj. vs time, 40 TRW iters, final gap = %g', dg40));
xlabel('Time (s)'); ylabel('Obj');

print -depsc fig/horseObjCurve40.eps

%% compare to j40 Runtime
gap = 0.05;
midpt = mean([mFval(end) j40Fval(end)]);
upper = midpt + gap*abs(midpt);
lower = midpt - gap*abs(midpt);

mSuboptTime   = mTime(find(mFval <= upper, 1))
j40SuboptTime = j40Time(find(j40Fval >= lower, 1))

speedup = j40SuboptTime / mSuboptTime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra junk we are not currently using.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% grrgh; relative duality gap
h4 = figure;     

semilogy(mTime,    mFval - min(mFval) + 1, '--k', ...
         j200Time, max(j200Fval) - j200Fval + 1, '.-b');
     
xlabel('Difference from Best Objective Value');
legend('BCFW to best primal', 'domke200 to best dual');
 
% yl = get(gca, 'ytick');
% set(gca, 'yticklabel', sign(yl).*10.^abs(yl));

% dg200 = mFval(end) - j200Fval(end);
% title(sprintf('Obj. vs time, 200 TRW iters, final gap = %g', dg200));
% xlabel('Time (s)'); ylabel('Obj');
% 
% print -depsc fig/horseObjCurve200.eps


%% compare to j200 Runtime
gap = 0.15;
midpt = mean([mFval(end) j200Fval(end)]);
upper = midpt + gap*abs(midpt);
lower = midpt - gap*abs(midpt);

mSuboptTime   = mTime(find(mFval <= upper, 1))
j200SuboptTime = j200Time(find(j200Fval >= lower, 1))

speedup = j200SuboptTime / mSuboptTime


%% ERR TIMES
% You must reproduce this statement:
% When comparing at same accuracy, on the horse dataset, our algorithm
% attained test error within 1 percentage point of the true error in just
% 11.2 min compared to 2.5h for [3], a 13.6 speedup.

gap = 0.01;

% bestErr  = min([mWavgErr ; j40Err]);
% upperErr = (1 + gap) * bestErr;

upperErr = (1 + gap) * min(j40Err);

mSuberrTime    = mTime(find(mWavgErr <= upperErr, 1)) / 3600
j40SuberrTime = j40Time(find(j40Err  <= upperErr, 1)) / 3600

speedup = j40SuberrTime / mSuberrTime


mSuberrTimeMin = mSuberrTime * 60
