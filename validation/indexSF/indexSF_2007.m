function SF = indexSF_2007(data,labels)
% SF = INDEXSF(data,labels)
%--------------------------------------------------------------------------
% Cluster internal validity index Score Function (Saitta et al., 2007).
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
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
% Saitta, S., Raphael, B., & Smith, I. (2007). A bounded index for cluster
% validity. In P. Perner (Ed.), Machine Learning and Data Mining in Pattern
% Recognition (pp. 174�187). Leipzig, Germany: Springer Berlin Heidelberg.
% doi:10.1007/978-3-540-73499-4_14
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 16-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


[N,D] = size(data);
[clusterIDs, ~, lbl] = unique(labels);
K = length(clusterIDs);

% normalize data
%data = pplk_normalize(data,'norm');
% data = pplk_normalize(data,'range');
% data = bsxfun(@rdivide,data, sqrt(sum(data.^2,2))); % normalization

% compute the center of the data
Ztot = mean(data,1);

% compute centroids of clusters
bcd = 0;
wcd = 0;
for i = 1:K
    C = data(lbl==i,:);
    Ni = size(C,1); % number of data points in cluster
    Zi = mean(C,1); % cluster center
    tmp = Zi - Ztot;
    
    bcd = bcd + sqrt(tmp*tmp') * Ni;
    wcd = wcd + sum(sqrt(sum(bsxfun(@minus,C,Zi).^2,2)))/Ni;
end

bcd = bcd / (N*K);

SF = 1 - 1/exp(exp(bcd-wcd));


