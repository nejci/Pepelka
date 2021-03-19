function [labels, moreInfo]=pplk_clustererGSOM(data,K,params)

if isempty(K)
    % stop the algorithm when at least 2 clusters are found
    K=2;
end

% defaults
G=0.0008;
dG=0.045;
alfa=0.01;
maxIter=100;
pOut=0.1;
msize=[];
shape='rect'; 
showSOM = 0;
showGrav = 0;
distance='sqEuclidean';
advanced={};

if exist('params','var') && isstruct(params)
    if isfield(params,'GSOM_G')
        G = params.GSOM_G;        
    end
    if isfield(params,'GSOM_dG')
        dG = params.GSOM_dG;        
    end
    if isfield(params,'GSOM_alfa')
        alfa = params.GSOM_alfa;        
    end
    if isfield(params,'GSOM_maxIter')
        maxIter = params.GSOM_maxIter;        
    end
    if isfield(params,'GSOM_pOut')
        pOut = params.GSOM_pOut;        
    end
    if isfield(params,'GSOM_msize')
        msize = params.GSOM_msize;        
    end
    if isfield(params,'GSOM_shape')
        shape = params.GSOM_shape;        
    end
    if isfield(params,'GSOM_showSOM')
        showSOM = params.GSOM_showSOM;        
    end
    if isfield(params,'GSOM_showGrav')
        showGrav = params.GSOM_showGrav;        
    end
    if isfield(params,'GSOM_distance')
        distance = params.GSOM_distance;        
    end
    if isfield(params,'GSOM_advanced')
        advanced = params.GSOM_advanced;        
    end    
end


if ~ismember(distance, {'sqEuclidean', 'correlation', 'dotprod','dotprodc'})
	error(['gSOM - Wrong distance:', distance]);
end

oldPath=chdir('GSOM');

retValGSOM = gSOM(data,G, dG, alfa,maxIter, pOut, msize, shape,...
				showSOM, showGrav,0, K, distance, advanced);

chdir(oldPath);
					
labels = retValGSOM.target;
moreInfo.time = retValGSOM.time(3);
moreInfo.nClust = retValGSOM.nClust;
moreInfo.nIters = retValGSOM.iter;
moreInfo.options.SOMprop = retValGSOM.SOMprop;
moreInfo.options.SOMtrain = retValGSOM.SOMtrain;
moreInfo.options.params = params;
% moreInfo.param.SOMprop=retValGSOM.SOMprop;
% moreInfo.param.SOMtrain=retValGSOM.SOMtrain;
% moreInfo.param.options=params;