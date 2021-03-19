function SF = indexSF(data,labels)
% SF = INDEXSF(data,labels)
%--------------------------------------------------------------------------
% Cluster internal validity index Score Function (Saitta et al., 2008).
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [N X D] with N D-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%--------------------------------------------------------------------------
% OUTPUTS:
%   SF          (scalar)	value of the index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Saitta, S., Raphael, B., & Smith, I. (2008). A comprehensive validity
% index for clustering. Intelligent Data Analysis, 12(6), 529�548.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 18-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


[N,~] = size(data);
K = max(labels);

% normalize data
%data = pplk_normalize(data,'norm');
% data = pplk_normalize(data,'range');
% data = bsxfun(@rdivide,data, sqrt(sum(data.^2,2))); % normalization

% compute the center of the data
Ztot = mean(data,1);

% Main loop: 
% 1. compute centroids of clusters
% 2. compute distances between centroids and data mean
% 3. compute distances between centroids and data points in each cluster
% 4. compute bcs and wcs and add them to cumulative sum

bcd = 0;
wcd = 0;
for i = 1:K
    C = data(labels==i,:); % data points in the cluster k
    Ni = size(C,1);     % number of data points in cluster
    Zi = mean(C,1);     % cluster centroid
    tmp = Zi - Ztot;    % prepare for Euclid. distance calculation
    
    % compute bcd for cluster k and add to sum
    bcd = bcd + tmp*tmp' * Ni; 
    % compute wcd for cluster k and add to sum
    wcd = wcd + sqrt(sum(sum(bsxfun(@minus,C,Zi).^2,2))/Ni);
end

% weight both quantities
bcd = bcd / (N*K);
wcd = wcd / K;

% return Score Function
SF = 1 - 1/exp(exp(bcd-wcd));


