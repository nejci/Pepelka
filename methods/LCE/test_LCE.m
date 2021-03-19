% Test script for DICLENS family

scoreExt = {'NMI'};

%% Load data
%[data,target] = pplk_loadData('artificial\D31');
[data,target] = pplk_loadData('real\iris_orig');
[N D] = size(data);
%% Generate ensemble
params = pplk_setParamsDefault();
%params.KM_nRuns = 10;
ensSize = 50;
K = floor(sqrt(N));

labelsEns = pplk_genEns(data,{'KM',ensSize,K,'rand'},params);

Kcons = max(target);

%% LCE
%(original implementation)
ticID1 = tic();
[labelsCons,numClust] = pplk_consEns(labelsEns,Kcons,'LCE_ORIG-SRS-SL',params);
t1 = toc(ticID1);

%pplk_scatterPlot(data,labelsCons);
[~,score1] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'LCE-ORIG | score: %f, time: %f\n', score1, t1);

%% LCE-fast
% fresh & faster implementation
ticID2 = tic();
[labelsCons,numClust]= pplk_consEns(labelsEns,Kcons,'LCE-SRS-SL',params);
t2 = toc(ticID2);
%pplk_scatterPlot(data,labelsCons);
[~,score2] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'LCE-fast | score: %f, time: %f\n', score2, t2);
