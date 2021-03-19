function [labels,moreInfo] = SOMSTAR(data,somsize,smooth,show)
% Clustering of SOM with Starburst method (Hamel & Brown, 2011)
%--------------------------------------------------------------------------
% data - matrix with an observation in each row
% K - dummy; number of clusters is determined automatically
% params.SOMSTAR_
%    somsize  - size of SOM grid; could be [], [xdim,ydim] or scaling factor -f
%    smooth - smoothing parameter, bandwidth of the Gaussian kernel (>=0)
%    show   - if 1, plot the starburst map and clusters of data
%--------------------------------------------------------------------------
% Reference: 
% L. Hamel and C. Brown, "Improved Interpretability of the Unified Distance
% Matrix with Connected Components," in 7th International Conference on
% Data Mining, 2011, pp. 338–343.
%--------------------------------------------------------------------------
% Writen by Nejc Ilc
% Last modification 3-January-2014
% Requires Pepelka package

% Add path to SOMtoolbox
path_somtoolbox = [pplk_homeDir('misc'),'somtoolbox_pplk'];
addpath(path_somtoolbox);

% Input arguments check
if ~exist('somsize','var') || isempty(somsize)
    somsize = -1;
end
if ~exist('smooth','var') || isempty(smooth)
    smooth = 1;
end
if ~exist('show','var') || isempty(show)
    show = 0;
end
% defaults
somShape = 'rect'; % only rectangular grid is supported
somAlg = 'batch';
somTrainlen = 'default';

% Train SOM
D = som_data_struct(data);
ticID = tic();
sM = som_make(D,...
    'msize',somsize,'lattice',somShape,'algorithm',somAlg,'shape','sheet',...
    'neigh','gaussian','init','lininit','training',somTrainlen,'tracking',0);
T1 = toc(ticID);

ticID = tic();
% Compute smoothened U-matrix
umat = computeUmat(sM,smooth);
[xdim,ydim] = size(umat);
% Find coordinates of neurons that are centres of stars
coords = computeInternalNodes(umat);
sb = [coords.xcoords(:),coords.ycoords(:)];
% identify cluster centers (which neurons are centers) and neurons labels
[clusterCenters,~,labelsNeurons] = unique(sb,'rows');

% map data points to their bmus neurons
bmus = som_bmus(sM,D);
labels = labelsNeurons(bmus);
[uLabels,~,labels] = unique(labels);

T2 = toc(ticID);

% remove path and cd back to parent dir
rmpath(path_somtoolbox);

% return additional info
moreInfo.nClust = length(uLabels);
moreInfo.time = T1+T2;
moreInfo.options.SOMprop = sM.topol;
moreInfo.options.SOMtrain = sM.trainhist;
moreInfo.options.smooth = smooth;

if show
    % plot the starburst (figures are similar to those obtained by R code,
    % i.e., X dimension is horizontal and Y is vertical)
    load('mycmap.mat');
    fig1 = figure();    
    hm = heatmap(umat',[],[],[],'Colormap',mycmap);
    hold on;
    title('Enhanced U-Matrix');
    X = repmat((1:xdim)',1,ydim);
    Y = repmat((1:ydim),xdim,1);
    Xm = reshape(sb(:,1),xdim,ydim);
    Ym = reshape(sb(:,2),xdim,ydim);
    for i=1:xdim
        for j=1:ydim
            plot([X(i,j),Xm(i,j)],[Y(i,j),Ym(i,j)],'k-');
        end
    end
    for c = 1:size(clusterCenters,1)
        plot(clusterCenters(c,1),clusterCenters(c,2),'.k','MarkerSize',20);
    end
    hold off;
    
    pplk_scatterPlot(data,labels);
    
end
