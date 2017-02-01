%% Example for learnBethe.
% You should see obj converge to the negative maximum likelihood and dMu
% converge to zero.

Areasonable = rand(5,5);
Mreasonable = sample_perms(Areasonable, 10000);
[thetaReasonable, histReasonable] = learnBethe(Mreasonable, @minBethe, 'thetaTrue', Areasonable);
