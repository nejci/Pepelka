function [labels, moreInfo]=pplk_clustererSOMKM(data,K,params)
%Clustering of SOM with K-means (Vesanto & Alhoniemi 2000)
%
% data - matrix with an observation in each row
% K - number of clusters to evaluate
% params.SOMKM_
%    run2toK - logical: 0-run KM for K clusters, 1-run KM for range 2..K
%         cluster and return solution which minimizes Davies-Bouldin index
%    nRuns - number of runs for each cluster number
%    msize - size of SOM grid; could be [], [xdim,ydim] or scaling factor -f
%    shape - 'rect' or 'hexa'
%
% Writen by Nejc Ilc
% Last modification 20-December-2013

% Add path to modified SOMtoolbox (faster PCA, bug fix for rect/hexa size)
path_somtoolbox = [pplk_homeDir('libs'),'somtoolbox_pplk'];
addpath(path_somtoolbox);

N = size(data,1);

% defaults
if ~exist('K','var') || isempty(K)
    K = ceil(sqrt(N));
    run2toK = 1;
else
    run2toK = 0;
end

nRuns = 5;
msize = -1; % default size
shape = 'rect';
show = 0;
if exist('params','var') && isstruct(params)
    if isfield(params,'SOMKM_run2toK')
        run2toK = params.SOMKM_run2toK;
    end
    if isfield(params,'SOMKM_nRuns')
        nRuns = params.SOMKM_nRuns;
    end
    if isfield(params,'SOMKM_msize')
        msize = params.SOMKM_msize;
    end
    if isfield(params,'SOMKM_shape')
        shape = params.SOMKM_shape;
    end
    if isfield(params,'SOMKM_show')
        show = params.SOMKM_show;
    end
end

D = som_data_struct(data);
% normalize data?
%D = som_normalize(D, 'range');

% create SOM
ticID = tic();       
sMap = som_make(D,...
    'msize',msize,'lattice',shape,'algorithm','batch','shape','sheet',...
    'neigh','gaussian','init','lininit','training','default','tracking',0);
tSOM = toc(ticID);

% Compute BMU for every data point
data_bmus = som_bmus(sMap,D);

% BMUs indices
bmus = unique(data_bmus);

% Extract codebook vectors of BMUs - those are data points to cluster with
% K-means
coords = sMap.codebook(bmus,:);

% Limit the desired number of clusters if it exceeds the number of neurons.
K = min(size(coords,1),K);
    
ticID = tic();
if run2toK    
    %K-means: for each K run K-means nRuns times and compute 
    [~,clusters,~,best_ind] = kmeans_clusters(coords, K, nRuns, show);
    [~,best_i] = min(best_ind); % select the one with smallest DB index
    
    c=clusters{best_i};
    nClust=best_i;
else   
    % Run K-means nRuns times and return labels of best solution, which is
    % determined by minimal error value (K-means internal cost function).
    c = kmeans(coords,K, 'distance','sqEuclidean', ...
        'start','sample','replicates', nRuns, 'emptyaction','singleton',...
        'display','off');
    nClust=K;
end

% Projection of the clustered BMUs back to the original data points
assert(length(c)==length(bmus), 'Oops.');
labBin = bsxfun(@eq,data_bmus,bmus');
dataLabels = labBin*c;

tKM=toc(ticID);

if(show)
    plotOpt.title = ['K-means of SOM for k=',num2str(nClust)];
    pplk_scatterPlot(data,dataLabels,nClust,plotOpt);
    
    dim = size(data,2);
    if dim > 2
        %PCA projection of neurons
        neur=[sMap.codebook ; coords];
        new_CB=(D.pca_vec(:,1:2)' * neur')';
        numCB=size(sMap.codebook,1);        
        CB=new_CB(1:numCB,:);
        suns_new=new_CB(numCB+1:end,:);
        new_data=(D.pca_vec(:,1:2)' * D.data')';        
    else
        CB=sMap.codebook;
        suns_new=coords;
        new_data=D.data;        
    end
    % plot neurons, interconnections, data points and mark BMUs with red 
    figure();
    som_grid(sMap, 'Coord', [CB(:,1),CB(:,2)]);
    axis('square')
    hold on;
    plot(new_data(:,1),new_data(:,2),'bo');    
    for ind=1:length(bmus)
        plot(suns_new(ind,1),suns_new(ind,2),'r.', 'markersize', 5);
    end
    hold off;    
end

% Gather results
labels = dataLabels;
moreInfo.nClust = nClust;
moreInfo.time = tSOM+tKM;
moreInfo.options.SOMprop = sMap.topol;
moreInfo.options.SOMtrain = sMap.trainhist;
moreInfo.options.nRuns = nRuns;
moreInfo.options.run2toK = run2toK;

% Remove SOMtoolbox from the path
rmpath(path_somtoolbox);