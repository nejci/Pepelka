function [labels,moreInfo] = pplk_clustererSPECLS2(eigv,K,params)

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

if exist('params','var') && isstruct(params)
    if isfield(params,'SPECLS_Knn')
        Knn = params.SPECLS_Knn;        
    end
    if isfield(params,'SPECLS_mode')
        mode = params.SPECLS_mode;
    end
end

ticID = tic();
[labels, numClust] = SPECLS2(eigv,K,mode);
time=toc(ticID);

moreInfo.time = time;
moreInfo.nClust = numClust;
moreInfo.options.Knn = Knn;
moreInfo.options.mode = mode;

chdir(oldPath);