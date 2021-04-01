function [labels,moreInfo] = pplk_clustererSOMNC(data,K,params)
% Clustering of SOM with Normalized Cuts (Yang et al. 2012)
%--------------------------------------------------------------------------
% data - matrix with an observation in each row
% K - number of clusters to evaluate
% params.SOMNC_
%    nRuns - number of runs for each cluster number
%    msize - size of SOM grid; could be [], [xdim,ydim] or scaling factor -f
%    shape - 'rect' or 'hexa'
%--------------------------------------------------------------------------
% Reference:
% L. Yang, Z. Ouyang, and Y. Shi, “A Modified Clustering
% Method Based on Self-Organizing Maps and Its Applications,” Procedia
% Computer Science, vol. 9, pp. 1371–1379, Jan. 2012.
%--------------------------------------------------------------------------
% Writen by Nejc Ilc
% Last modification 3-January-2014

path_somtoolbox = [pplk_homeDir('libs'),'somtoolbox_pplk'];
addpath(path_somtoolbox);

somSize = -1;
somShape = 'rect';
somAlg = 'batch';
somTrainlen = 'default';

if exist('params','var') && isstruct(params)
    if isfield(params,'SOMNC_msize')
        somSize = params.SOMNC_msize;
    end
    if isfield(params,'SOMNC_shape')
        somShape = params.SOMNC_shape;
    end
end

% create SOM
D = som_data_struct(data);
ticID = tic();
sM = som_make(D,...
    'msize',somSize,'lattice',somShape,'algorithm',somAlg,'shape','sheet',...
    'neigh','gaussian','init','lininit','training',somTrainlen,'tracking',0);
T1 = toc(ticID);

% run the Ncut algorithm on a graph made from SOM
NcutPath = pplk_homeDir('methods\NC');
oldPath = chdir(NcutPath);

somSize = sM.topol.msize;
xdim = somSize(1);
ydim = somSize(2);

% neighborhood matrix of SOM
Ne = som_neighbors(sM,'N1');
% compute distances between neighbors
DIST = sqrt(sqdistance2(sM.codebook));
DIST = Ne .* DIST;

% turn distances into similarities
W = pplk_dist2sim(DIST,'lin2',1);

options.offset = 0.5; %offset in the diagonal of W, default 0.5
options.verbose = 1; %0 for verbose off mode, 1,2,3 for verbose on modes
options.maxiterations = 100; %max number of iterations in eigensolver
options.eigsErrorTolerance = 1e-8; %error tolerance in eigensolver
options.valeurMin=1e-6; %truncates any values in W less than valeurMin

ticID = tic();
NcutDiscrete = ncutW(W,K,options);
T2 = toc(ticID);

%transform n-outOf-1 presentation
labelsNeurons = zeros(xdim*ydim,1);
for labInd=1:size(NcutDiscrete,2)
    ind = logical(NcutDiscrete(:,labInd));
    labelsNeurons(ind)=labInd;
end

% Projection of the clustered neurons back to the original data points
% Compute BMU for every data point
bmus = som_bmus(sM,D);
labels = labelsNeurons(bmus);
[uLabels,~,labels] = unique(labels);

% remove path and cd back to parent dir
chdir(oldPath);
rmpath(path_somtoolbox);

% return additional info
moreInfo.nClust = length(uLabels);
moreInfo.time = T1+T2;
moreInfo.options.SOMprop = sM.topol;
moreInfo.options.SOMtrain = sM.trainhist;
