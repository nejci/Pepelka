function [labels, numClust, distBest, distHist] = SPECLS2(eigv,K,mode)
% SPECLS
% MODIFIED for some speed-up (precomputed eigen values instead of data)
%
% [labels] = SPECLS(data,K,Knn,mode)
% Spectral clustering with local scaling (Zelnik et al., 2004)
%--------------------------------------------------------------------------
% INPUTS
%   data        (matrix)	matrix [N X d] with N d-dimensional samples
%
%   K           (scalar)	number of clusters to find in data.
%               (vector)    multiple numbers of clusters - automatically
%                           search for the true one among them. It works
%                           only when mode is 'RLS' or 'R'.
%	                        Leave empty [] to set the interval on [2,sqrt(n)]
%
%   Knn         (scalar)    number of neighbors to consider in local scaling;
%                           default is 7.
%   mode        (string)    'LS'  Locally Scaled clustering
%                           'RLS' Rotation clustering with local scaling
%                           'R'
%                           if mode and K are empty, 'RLS' is selected
%
%
% OUTPUTS
%   labels      (vector)    data labels
%               (matrix)    if CVI has more than one element and wMode is
%                           'single', consensus is computed for every index
%   numClust    (scalar)    number of output clusters
%--------------------------------------------------------------------------
% EXAMPLE
%
% [data,target] = pplk_loadData('rings3');
% K_true = max(target);
% labels_true = SPECLS(data,K_true);
% labels_auto = SPECLS(data,[]);
% labels_custom = SPECLS(data,2,10,'RLS');
%
% options.title = 'Spectral local scaling clustering';
% options.subtitle = {'True', 'Auto', 'Custom'};
% pplk_scatterPlot(data,[labels_true,labels_auto,labels_custom],[],options)
%
%------- REFERENCE --------------------------------------------------------
% Zelnik-Manor, L., & Perona, P. (2004). Self-tuning spectral clustering.
% Advances in neural information processing systems, 2, 1601�1608.
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 27-September-2013 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

% Maximum number of clusters to automatically search to 
MAX_CLUSTERS_AUTO = 50; 

N = size(eigv,1);

if ~exist('K','var')
    K = [];
end
if isempty(K)
    % automatically search on default interval for true number of cluster
    % limit upper bound of clusters
    bound = min(MAX_CLUSTERS_AUTO, ceil(sqrt(N))); 
    K = 2:bound; 
end

if isscalar(K)
    K_isvec = 0;
else
    K_isvec = 1;
end


if ~exist('mode','var') || isempty(mode)
    if K_isvec
        mode = 'RLS';
    else
        mode = 'LS';
    end
end

if K_isvec && strcmpi(mode,'LS')
   warning('pplk:SPECLS','Cannot use mode LS when K is a vector or empty. Switching mode to RLS.'); 
   mode = 'RLS';
end

% % centralize and scale the data
% data = bsxfun(@minus,data,mean(data,1));
% data = data/max(abs(data(:)));
% 
% % Build affinity matrix A
% D = dist2(data,data);    % Euclidean distance
% [~,A_LS] = scale_dist(D,Knn); % Locally scaled affinity matrix
% 
% % Zero out diagonal
% ZERO_DIAG = ~eye(size(data,1));
% A_LS = A_LS.*ZERO_DIAG;

% Clustering
switch upper(mode)
    case 'LS'
        % Zelnik-Perona Locally Scaled clustering
        [labels,distBest,distHist] = gcut2(eigv,K);        
        numClust = max(labels);
        
    case 'RLS'
        error('RLS not supported yet by SPECLS2.');
        
    otherwise
        error('Wrong mode!');
end



