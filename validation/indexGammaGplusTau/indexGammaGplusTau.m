function [G, Gmod, Gplus, Tau] = indexGammaGplusTau(D,labels,Dmode,dtype)
% [G, Gmod, Gplus, Tau] = INDEXGAMMA(data,labels)
%--------------------------------------------------------------------------
% Cluster internal validity indices Gamma (+mod), G+ and Tau.
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
%   G           (scalar)	value of the index Gamma
%   Gmod        (scalar)	value of the modified index Gamma
%   Gplus       (scalar)	value of the index Gplus
%   Tau         (scalar)	value of the index Tau
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCES -------------------------------------------------------
% Baker, F. B., & Hubert, L. J. (1975). Measuring the Power of Hierarchical
% Analysis Cluster. Journal of the American Statistical Association,
% 70(349), 31�38.
%
% Arbelaitz, O., Gurrutxaga, I., Muguerza, J., P�rez, J. M., & Perona, I.
% (2013). An extensive comparative study of cluster validity indices.
% Pattern Recognition, 46(1), 243�256. doi:10.1016/j.patcog.2012.07.021
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
   % D is a dissimilarity matrix
   if ~isvector(D)
       D = D(tril(true(size(D)),-1))';
   end    
end

N = length(labels);
Nt = (N*(N-1))/2; % number of all pairs

K = max(labels);

% construct binary vector B; Bij is 1 if pair of data points i and j fall
% into the same cluster, and 0 otherwise.
B = false(1,Nt);

Nk = zeros(1,K); % number of data points in each cluster

cum = cumsum(1:N);
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
    B(ind) = 1;
end

% number of all pairs that share the same cluster
Nw = 0.5 * (sum(Nk.^2)-N);
% number of all pairs that do not share the same cluster
Nb = Nt - Nw;

% compute concordances s+ and s-
Bzero = find(~B);

sPlus = 0;
sMinus = 0;

for b=Bzero
    sPlus = sPlus + sum(D(b) < D(B));
    sMinus = sMinus + sum(D(b) > D(B));
end

G = (sPlus - sMinus) / (sPlus + sMinus); % max is better
Gplus = 2*sMinus / (Nt*(Nt-1)); % min is better
Tau = (sPlus - sMinus) / sqrt(Nb*Nw*Nt*(Nt-1)*0.5); % max is better

% Alternative index as described in Arbelaitz et al., 2013
Gmod = sPlus / (Nw * (Nt-Nw)); % min is better


