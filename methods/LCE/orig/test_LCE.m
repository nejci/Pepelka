% Test script for DICLENS family

scoreExt = {'NMI'};

%% Load data
[data,target] = pplk_loadData('artificial\halfring');
%[data,target] = pplk_loadData('real\iris_orig');
[N D] = size(data);
%% Generate ensemble
params = pplk_setParamsDefault();
%params.KM_nRuns = 10;
ensSize = 10;
K = floor(sqrt(N));

labelsEns = pplk_genEns(data,{'KM',ensSize,K,'rand'},params);

Kcons = max(target);

%% LCE
%(original implementation)
ticID = tic();
[labelsCons,numClust] = pplk_consEns(labelsEns,Kcons,'LCE',params);
t1 = toc(ticID);

pplk_scatterPlot(data,labelsCons);
[~,score1] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'LCE-ORIG | score: %f, time: %f\n', score1, t1);

%% LCE-fast
% fresh & faster implementation
ticID = tic();
[labelsCons,numClust]= pplk_consEns(labelsEns,Kcons,'LCE2',params);
t2 = toc(ticID);
pplk_scatterPlot(data,labelsCons);
[~,score2] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'DICLENS-fast | score: %f, time: %f\n', score2, t2);
