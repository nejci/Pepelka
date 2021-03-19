function VAR = indexVAR(data,labels)
% VAR = indexVAR(data,labels,alpha)
%--------------------------------------------------------------------------
% Cluster internal validity index - Variance. 
% Index should be minimized.
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%--------------------------------------------------------------------------
% OUTPUTS:
%   VAR         (scalar)	value of the Variance index 
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Handl, J., Knowles, J., & Kell, D. B. (2005). Computational cluster
% validation in post-genomic data analysis. Bioinformatics, 21(15),
% 3201�3212.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 23-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

[N,D] = size(data);
K = max(labels);

% compute centers and variance of clusters
V = zeros(K,D);
var = zeros(K,1);
for i = 1:K
    dataK = data(labels==i,:);
    V(i,:) = mean(dataK,1);
    
    % compute squared Euclidean distance (when q=1)
    var(i) = sum(sum(bsxfun(@minus,dataK,V(i,:)).^2,2));
end

VAR = sqrt(sum(var)/N);