function [HOM,HOMmin,SEP,SEPmax,WR,WRp] = indexHOM_SEP_WR(S,labels,Smode,dtype)
% [HOM,HOMmin,SEP,SEPmax,WR,WRp] = INDEXHOM_SEP_WR(S,labels,Smode,dtype)
%--------------------------------------------------------------------------
% Cluster internal validity indices: Homogeneity, Seperation, Weighted
% inter-intra ratio.
%--------------------------------------------------------------------------
% INPUTS
%   S   		(matrix)	data matrix [n X d] with n d-dimensional samples
%               (matrix)    dissimilarity/distance matrix [n X n]
%               (matrix)    similarity matrix [n X n]
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%   Smode       (string)    'data' - S is data matrix
%                           'dist' - S is distance matrix
%                           'sim'  - S is similarity matrix
%   dtype       (string)    distance identifier that is passed to pdist
%                           when Smode is 'data' or is used to compute
%                           similarities when Smode is 'dist'. Can be:
%                           'euclidean', 'cityblock', 'minkowski',
%                           'chebychev', 'correlation', 'spearman',
%                           'cosine', 'hamming','jaccard'.
%--------------------------------------------------------------------------
% OUTPUTS:
%   HOM          (scalar)	value of the average Homogeneity index
%   HOMMIN       (scalar)   value of the minumum Homogeneity index (worst case)
%   SEP          (scalar)	value of the average Separation index
%   SEPMAX       (scalar)   value of the maximum Separation index (worst case)
%   WR           (scalar)	value of the Weighted inter-intra ratio index
%   WRP          (scalar)   value of the Weighted inter-intra ratio index, 
%                           large number of clusters is penalized
%                           
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCEs -------------------------------------------------------
% Sharan, R., Maron-Katz, A., & Shamir, R. (2003). CLICK and EXPANDER: a
% system for clustering and visualizing gene expression data.
% Bioinformatics, 19(14), 1787�1799.
%
% Strehl, A. (2002). Relationship-based clustering and cluster ensembles
% for high-dimensional data mining. PhD thesis.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 10-August-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================
if ~exist('dtype','var') || isempty(dtype)
    dtype = 'euclidean';
end

if ~exist('Smode','var') || isempty(Smode)
    Smode = 'data';
end

if ~any(strcmpi(Smode,{'data','dist','sim'}))
    error('Wrong Smode value!');
end

% if S is data matrix compute distances
if strcmpi(Smode,'data')   
    S = squareform(pdist(S,dtype));
    Smode = 'dist';
end

% if S is a distance matrix compute similarities
if strcmpi(Smode,'dist')
    switch dtype
        case {'euclidean', 'cityblock', 'minkowski', 'chebychev'}
            S = 1./(1+S);
        case {'correlation', 'spearman'}
            S = 2-S; %correlation coeficient on interval [0,2] instead of [-1,1]
        case {'cosine', 'hamming','jaccard'}
            S = 1-S; %cosine coeficient on interval [0,1]
        otherwise
            error(['Do not know how to calculate similarities from distance ', dtype])
    end
end

N = size(S,1);
K = max(labels);

% % make diagonal elements 0
S(1:N+1:N^2)=0;

E = eye(K);
U = logical(E(:,labels));

% number of points in a cluster
clsize = sum(U,2);
% sum all similarities between points in the same cluster
sumIntra = sum(U.*(U*triu(S,1)),2);

M = (clsize.^2-clsize)./2;

% Homogeneity index
HOM = sum(sumIntra) / sum(M);
HOMmin = min(sumIntra ./ M);

% Separation index
if nargout > 2
    sumInter = zeros(1,(K^2-K)/2);
    idx = 1;

    for ci = 1:K-1
        for cj = ci+1:K
            sumInter(idx) = sum(sum(S(U(ci,:),U(cj,:))))/2;
            idx = idx+1;
        end
    end
    
    Cmat = clsize * clsize';
    Cmat = Cmat(tril(true(K),-1))';
    
    SEP = sum(sumInter) / ( (N^2-N)/2 - sum(M));
    SEPmax = max(sumInter ./ Cmat);    
end

% Weighted inter-intra ratio index
if nargout > 4
    nom = sum(sum(squareform(sumInter),2)./(N-clsize));
    singleton = clsize==1;
    if any(singleton)
        warning('pplk:validInt','WR: singleton clusters detected.');
    end    
    clsize(singleton) = 2; % singleton clusters will also be included
    denom = sum(2./(clsize-1).*sumIntra);
    WR = 1- nom/denom;
    % Penalized index to favorize low K.
    WRp = (1-2*K/N) * WR;    
end

