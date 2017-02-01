% NOTE: We are using config2 for my results. (Better duality gap, and early
% behavior does not change.)
nBeatHorse = 2;

% Confusing stuff:
% jHorseTestErr:
%   result1, result2 correspond to justinHorseCorrect/result2.mat (200)
%   result3, result4 correspond to justinHorseCorrect/result1.mat (40)

resJ200Horse = [1 2]; 
resJ40Horse  = [3 4]; % result2

resJ20Horse  = [5 6];
resJ10Horse  = [7 8]; % eventually will need to add more

lenBeat    = 250;
% This is all just dirty hackery...
lenJ200 = 100;
lenJ40  = 96;
lenJ20  = 97;
lenJ10  = 200; % just for now

beatErr     = zeros(lenBeat, 1);
beatWavgErr = zeros(lenBeat, 1);

for n = 1:nBeatHorse
    st  = load(sprintf('expt/horseReproduceBeatTimingErr/config%d_cp.mat', n));
    who = st.errs ~= 0;
    beatErr(who)     = st.errs(who);
    beatWavgErr(who) = st.wavgErrs(who);
end

j200Err = zeros(lenJ200, 1);
for n = resJ200Horse
    st  = load(sprintf('expt/jHorseTestErr/result%d.mat', n));
    who = st.errs ~= 0;
    j200Err(who) = st.errs(who);
end

j40Err = zeros(lenJ40, 1);
for n = resJ40Horse
    st  = load(sprintf('expt/jHorseTestErr/result%d.mat', n));
    who = st.errs ~= 0;
    j40Err(who) = st.errs(who);
end

j20Err = zeros(lenJ20, 1);
for n = resJ20Horse
    st  = load(sprintf('expt/jHorseTestErr/result%d.mat', n));
    who = st.errs ~= 0;
    j20Err(who) = st.errs(who);
end

j10Err = zeros(lenJ10, 1);
for n = resJ10Horse
    st  = load(sprintf('expt/jHorseTestErr/result%d.mat', n));
    who = st.errs ~= 0;
    j10Err(who) = st.errs(who);
end

save('-v7.3', 'cp/horseTestErr.mat', 'beatErr', 'beatWavgErr', 'j200Err', 'j40Err', 'j20Err', 'j10Err');
% save('-v7.3', 'cp/horseTestErr.mat', 'mErr', 'j200Err', 'j40Err', 'j20Err');
