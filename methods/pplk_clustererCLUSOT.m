function [labels, moreInfo]=pplk_clustererCLUSOT(data,K,params)

% [labels, moreInfo]=pplk_clustererCLUSOT(data,K,params)
% Cluster data with CLUSOT algorithm (Brugger et al., 2008)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modification: 6. November 2013 
% (C) Pepelka Package, Nejc Ilc (nejc.ilc@fri.uni-lj.si)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldPath=chdir('CLUSOT');
% K is not relevant, number of clusters are determined automatically

% defaults
clustMethod = 'rec_flood';
somSize = -1;
somGrid = 'hexa';
somAlg = 'batch';
somTrainLen = 'default';
res = 0.05;
g = 0.5;
theta0 = 0.3;
theta  = 0.7;

if exist('params','var') && isstruct(params)
    if isfield(params,'CLUSOT_clustMethod'), clustMethod = params.CLUSOT_clustMethod; end
    if isfield(params,'CLUSOT_somSize'), somSize = params.CLUSOT_somSize; end
    if isfield(params,'CLUSOT_somGrid'), somGrid = params.CLUSOT_somGrid; end
    if isfield(params,'CLUSOT_somAlg'), somAlg = params.CLUSOT_somAlg; end
    if isfield(params,'CLUSOT_somTrainLen'), somTrainLen = params.CLUSOT_somTrainLen; end
    if isfield(params,'CLUSOT_res'), res = params.CLUSOT_res; end
    if isfield(params,'CLUSOT_g'), g = params.CLUSOT_g; end
    if isfield(params,'CLUSOT_theta0'), theta0 = params.CLUSOT_theta0; end
    if isfield(params,'CLUSOT_theta'), theta = params.CLUSOT_theta; end
end

ticID = tic();
% WARNING: labels may contain zeros - those points are considered outliers
[labels, mI] = clusot(data,'cluster_method',clustMethod,'som_size',somSize,...
    'som_grid',somGrid,'som_alg',somAlg,'som_trainlen',somTrainLen,...
    'res',res,'g',g,'theta',theta,'theta0',theta0);
time=toc(ticID);

moreInfo.time = time;
moreInfo.nClust = mI.numClust;
moreInfo.SOMstruct = mI.SOMstruct;
moreInfo.options = mI.params;

chdir(oldPath);