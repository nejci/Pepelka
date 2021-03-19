function [labels,moreInfo] = pplk_clustererSPECLS3SUB(data,K,params)

% [labels, moreInfo]=pplk_clustererSPECLS(data,K,params)
% Cluster data with Spectral local scaling algorithm (Zelnik-Perona, 2004)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modification: 28. October 2013 
% (C) Pepelka Package, Nejc Ilc (nejc.ilc@fri.uni-lj.si)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldPath=chdir('SPECLS');

% defaults
Knn = [];
mode = [];
KMruns = [];
KMiters = [];

if exist('params','var') && isstruct(params)
    if isfield(params,'SPECLS_Knn')
        Knn = params.SPECLS_Knn;        
    end
    if isfield(params,'SPECLS_mode')
        mode = params.SPECLS_mode;
    end
    if isfield(params,'SPECLS_KMruns')
        KMruns = params.SPECLS_KMruns;
    end
    if isfield(params,'SPECLS_KMiters')
        KMiters = params.SPECLS_KMiters;
    end
end

ticID = tic();
V = SPECLS_eigv(data,Knn);
[labels, numClust] = SPECLS3(V,K,mode,KMruns,KMiters);
time=toc(ticID);

moreInfo.time = time;
moreInfo.nClust = numClust;
moreInfo.options.Knn = Knn;
moreInfo.options.mode = mode;
moreInfo.options.KMruns = KMruns;
moreInfo.options.KMiters = KMiters;

chdir(oldPath);