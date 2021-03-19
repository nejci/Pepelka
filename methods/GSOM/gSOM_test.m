

[dataOrig, target, info]= pplk_loadData('carcinomas');
data = pplk_normalize(dataOrig,'zscore');

K = length(unique(target));

params=[];
params.GSOM_G = 0.0008; %0.0008;
params.GSOM_dG = 0.0; %0.045;
params.GSOM_alfa = 0.01;
params.GSOM_maxIter = 0;
params.GSOM_pOut = 0.1;
params.GSOM_msize = -2;%-1;
params.GSOM_shape = 'hexa';%'rect'; 
params.GSOM_showSOM = 0;
params.GSOM_showGrav = 0;
params.GSOM_distance='sqEuclidean';
params.GSOM_advanced={'eps',1e-3};

labels = pplk_runClusterer('gSOM',data,K,1,params);
pplk_validExt(target,labels,{'BCA'})

pplk_scatterPlot(data, labels);
pplk_scatterPlot(data, target);

