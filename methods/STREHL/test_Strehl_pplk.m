% test Strehl
eCVI = {'ARI'};

[data,target] = pplk_loadData('iris'); % try also iris_orig
data = pplk_normalize(data,'zscore');

Kcons = max(target);


%%
params = pplk_setParamsDefault();
% custom ensemble
labelsEns = pplk_genEns(data,{'KM',5,[12,32],'random'; 'random', 5, [12,32], 'random'},params);
%labelsEns = pplk_genEns(data,{'AL',2,Kcons,'fixed'; 'random',2,Kcons,'fixed'},params);

%% CSPA
tic();
cls = labelsEns';
clbs = clstoclbs(cls);
s = clbs' * clbs;
s = checks(s./size(cls,1));
toc();

tic();
CO = computeCO(labelsEns,[],'labels','full',[]);
toc();
assert(isequal(triu(s,1),triu(CO,1)));

labelsConsOrig = metis(s,Kcons);
labelsConsCO = metis(CO,Kcons);

labelsConsWrap = clusterensemble(labelsEns',Kcons,'cspa')';

[~,score_orig] = pplk_validExt(target,labelsConsOrig',eCVI)
[~,score_CO] = pplk_validExt(target,labelsConsCO',eCVI);
[~,score_wrap] = pplk_validExt(target,labelsConsWrap',eCVI);

assert((score_orig == score_CO) & (score_orig == score_wrap));

%% weights

% compute CVI
CVI = {'DN','DB','SDBW','CH','SIL','I','XB'};
%CVI = {'DN','DB','SDBW','CH','SIL','I','XB','XB*','DB*','DNg','DNs','CON','SD','SF','CI','GAMMA'};

[validStrct, list] = pplk_validInt(data, labelsEns, CVI, []);
PRM_raw = list';

%%
unifyMethod = 'minmax';
reduceMethod = 'FEKM';
reduceDim = 'DANCoFit';
weightMethod = 'wMean2';
weightMode = 'RRA';

PRMopt = [];
PRMopt.CVImat = PRM_raw;
PRMopt.reduceDim = reduceDim;

[weights,PRM,featInd] = ...
    pplk_partitionRelevance(data,labelsEns,CVI,...
    unifyMethod,reduceMethod,weightMethod,weightMode,PRMopt);

[N,numEns] = size(labelsEns);
%weights = ones(numEns,1)*0;

tic();
COw = computeCO(labelsEns,[],'labels','full',weights);
toc();
COw(1:N+1:N^2) = 1;
%spy(COw);
s2 = checks(COw);

labelsConsW = metis(COw,Kcons);

score_w = pplk_validExt(target,labelsConsW',eCVI)
%pplk_scatterPlot(data,labelsCons);

%% pepelka wrapper

params = [];
params.WEAC_CVI = CVI;
params.WEAC_CVImat = PRM_raw;
params.WEAC_unifyMeth = unifyMethod;
params.WEAC_reduceMeth = reduceMethod;
params.WEAC_reduceDim = reduceDim;
params.WEAC_weightMeth = weightMethod;
params.WEAC_weightMode = weightMode;
labelsConsPPLK = pplk_consEns(labelsEns,Kcons,'STREHL-MCLA',params);
score_wPplk = pplk_validExt(target,labelsConsPPLK',eCVI)

%% test easy
labelsEns2 = [1 1 2 2; 1 2 1 2]';
K = 2;
labelsConsHGPA = pplk_consEns(labelsEns2,K,'STREHL-HGPA',params);