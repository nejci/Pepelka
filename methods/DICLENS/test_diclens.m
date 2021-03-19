% Test script for DICLENS family

scoreExt = {'ARI'};
SHOW = 0;
%% Load data
[data,target] = pplk_loadData('face320');
data = pplk_normalize(data,'zscore');
%[data,target] = pplk_loadData('real\iris_orig');
[N D] = size(data);
%% Generate ensemble
params = pplk_setParamsDefault();
%params.KM_nRuns = 10;
ensSize = 10;
K = floor(sqrt(N));

labelsEns = pplk_genEns(data,{'KM',ensSize,K,'rand'},params);

% ensemble accuracy
scoreEns = zeros(1,ensSize);
for i = 1:ensSize
   [~,scoreEns(i)] = pplk_validExt(target,labelsEns(:,i),scoreExt); 
end

fprintf(1,'Ensemble: min=%f, max=%f, mean=%f\n',min(scoreEns),max(scoreEns),mean(scoreEns));

Kcons = max(target);
Kcons = []; % To automatically determine

%% DICLENS-ORIG
%(original implementation with errors)
ticID = tic();
[labelsCons,numClust] = pplk_consEns(labelsEns,Kcons,'DICLENS-ORIG',params);
t1 = toc(ticID);

[~,score1] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'DICLENS-ORIG | score: %f, time: %f\n', score1, t1);

if SHOW
    pplk_scatterPlot(data,labelsCons);
end

%% DICLENS-fast
% fresh & faster implementation
ticID = tic();
[labelsCons,numClust]= pplk_consEns(labelsEns,Kcons,'DICLENS',params);
t2 = toc(ticID);
[~,score2] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'DICLENS-fast | score: %f, time: %f\n', score2, t2);
if SHOW
    pplk_scatterPlot(data,labelsCons);
end

%% DICLENS-W
% using weights

% pre-compute CVI if you want to
ticID = tic();
CVI = {'DN','DNG','DNS','DB','SDBW','CH','SIL','SSI','I','XB'};
[validStrct, list] = pplk_validInt(data, labelsEns, CVI, []);
PRM_raw = list';

%
params = [];
params.WEAC_CVI = CVI;
params.WEAC_data = data;
params.WEAC_CVImat = PRM_raw;
params.WEAC_unifyMeth = 'minmax';
params.WEAC_reduceMeth = 'SPEC';
params.WEAC_reduceDim = 'DANCo';
params.WEAC_weightMeth = 'wRankAggreg2';
params.WEAC_weightMode = 'RRA';

%% DICLENS-W slow (alternative, different from fast)
ticID = tic();
[labelsCons,numClust]= pplk_consEns(labelsEns,Kcons,'DICLENS-WALT',params);
t4 = toc(ticID);

[~,score4] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'DICLENS-Walt    | score: %f, time: %f\n', score4, t4);
if SHOW
    pplk_scatterPlot(data,labelsCons);
end


%% DICLENS-W flow (based on fast)
ticID = tic();
[labelsCons,numClust]= pplk_consEns(labelsEns,Kcons,'DICLENS-W',params);
t4f = toc(ticID);

[~,score4f] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'DICLENS-W    | score: %f, time: %f\n', score4f, t4f);
if SHOW
    pplk_scatterPlot(data,labelsCons);
end


%% OTHERS
KconsT = max(target);
ticID = tic();
[labelsCons,numClust]= pplk_consEns(labelsEns,KconsT,'STREHL-MCLA',params);
tMCLA = toc(ticID);
pplk_scatterPlot(data,labelsCons);
[~,scoreMCLA] = pplk_validExt(target,labelsCons,scoreExt);
fprintf(1,'STREHL_MCLA | score: %f, time: %f\n', scoreMCLA, tMCLA);
if SHOW
    pplk_scatterPlot(data,labelsCons);
end
