SD = enPruneChannels(d,SD,tIncMan,[0.003           10000000],1,[0  45],0);

dod = hmrIntensity2OD(d);

[tIncAuto,tIncChAuto] = hmrMotionArtifactByChannel(dod,t,SD,tIncMan,0.5,1,4,0.4);

[dod] = hmrMotionCorrectSpline(dod,t,SD,tIncChAuto,0.99,1);

[tIncAuto,tIncChAuto] = hmrMotionArtifactByChannel(dod,t,SD,tIncMan,0.5,1,4,0.4);

[s,tRangeStimReject] = enStimRejection(t,s,tIncAuto,tIncMan,[-1  1]);

dod = hmrBandpassFilt(dod,t,0.01,1);

dc = hmrOD2Conc(dod,SD,[5.1         5.1]);

[dcAvg,dcAvgStd,tHRF,nTrials,dcSum2,dcTrials] = hmrBlockAvg(dc,s,t,[-2  15]);

