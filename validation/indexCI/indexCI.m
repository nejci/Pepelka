function [CI,CImod] = indexCI(D,labels,Dmode,dtype)
% [CI, CImod] = indexCI(data,labels)
%--------------------------------------------------------------------------
% Cluster internal validity index 'C-index' 
%--------------------------------------------------------------------------
% INPUTS
%   D   		(matrix)	data matrix [n X d] with n d-dimensional samples
%               (matrix)    dissimilarity matrix [n X n] or vector [(n*(n-1))/2,1]
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%   Dmode       (bool)      0 - data is data matrix
%                           1 - data is dissimilarity matrix/vector
%   dtype       (string)    distance identifier that is passed to pdist;
%                           only relevant when Dmode is 0
%--------------------------------------------------------------------------
% OUTPUTS:
%   CI          (scalar)	value of the C-index
%   CImod       (scalar)    value of alternative C-index
%                           (as implemented in R package NbClust)
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Dalrymple-Alford, E. C. (1970). Measurement of clustering in free recall.
% Psychological Bulletin, 74(1), 32�34. doi:10.1037/h0029393
%
% Hubert, L. J., & Levin, J. R. (1976). A general statistical framework for
% assessing categorical clustering in free recall. Psychological Bulletin,
% 83(6), 1072�1080. doi:10.1037//0033-2909.83.6.1072
%
%------- VERSION ----------------------------------------------------------
% Version: 1.1
% Last modified: 8-August-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if ~exist('Dmode','var') || isempty(Dmode)
   Dmode = 0; 
end

if ~exist('dtype','var') || isempty(dtype)
   dtype = 'euclidean'; 
end

if ~Dmode
    % D is a data matrix
    % construct vector D - pairwise distances between points
    D = pdist(D,dtype);
else
   % D is a similarity matrix
   if ~isvector(D)
       D = D(tril(true(size(D)),-1))';
   end    
end

N = length(labels);
Nt = (N*(N-1))/2; % number of all pairs of data points
K = max(labels);

% construct binary vector Q; Qij is 1 if pair of data points i and j fall
% into the same cluster, and 0 otherwise. Here we use linear indices of i
% and j.
Q = false(1,Nt);
Nk = zeros(1,K); % number of data points in each cluster

cum = cumsum(1:N);

vmax = zeros(1,K);
vmin = zeros(1,K);

for i = 1:K
    C = find(labels==i);
    Nk(i) = length(C);
    if length(C) ==1
        continue;
    end
    I = nchoose2(C);
    J = I(:,2)';
    I = I(:,1)';
    ind = (I-1)*N + J - cum(I);
    Q(ind) = 1;
    
    t = D(ind);
    vmax(i) = max(max(t));
    vmin(i) = min(t);
end

% number of all pairs that share the same cluster
Nw = 0.5 * (sum(Nk.^2)-N);

Sw = Q * D';
Dsort = sort(D,'ascend');
Smin = sum(Dsort(1:Nw)); % sum of Nw smallest point-to-point distances
Smax = sum(Dsort(end-Nw+1:end)); % sum of Nw largest point-to-point distances

CI = (Sw - Smin) / (Smax - Smin);
if nargout == 1
    return;
end
% Variant as in NbClust
Smin2 = min(vmin);
Smax2 = max(vmax);
CImod = (Sw - Nw * Smin2)/(Smax2 * Nw - Smin2 * Nw);