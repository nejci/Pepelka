function [weights,PRM,featInd]=pplk_partitionRelevance(data,labelsEns,indicesList,unifyMethod,reduceMethod,weightMethod,weightMode,options)
% [weights,PRM,featInd]=pplk_partitionRelevance(data,labelsEns,indicesList,
% unifyMethod,reduceMethod,weightMethod,weightMode,options)
% Partition Relevance step in clustering ensembles.
%
% INPUTS
%   data
%       A N-by-D matrix of data, where N is number of data samples and D is
%       number of dimensions. Can be empty if options.CVImat is provided
%       (values of cluster validation indices have been pre-computed).
%
%   labelsEns   
%       Cluster ensemble - labels are stored in the columns of an N-by-E
%       matrix, where N is number of data points and E is number of
%       clusterings.
%
%   indicesList
%       A cell containing one or more cluster validity indices:
%
%       Legend:
%       ID (Method name)[Best] Parameters
%
%       Optionally, you can use multiple methods, i.e., {'DN','SIL','WR'}.
%       If empty, none of them is calculated.
%
%       - 'APN' (Avg. proportion of non-overlap)[min] clMethod | labelsDel
%       - 'AD' (Average distance)[min] clMethod | labelsDel
%       - 'ADM' (Average distance between means)[min] clMethod | labelsDel
%       - 'BHI' (Biological homogeneity index)[max] genenames, annotations,
%         GO_aspect, GO_evidence
%       - 'BSI' (Biological stability index)[max] same as BHI
%       - 'CH' (Calinski-Harabasz) [max]
%       - 'CI' (C-index) [min]
%       - 'CON' (Connectivity index) [min] CON_L
%       - 'DB' (Davies-Bouldin index) [min] DB_p, DB_q
%       - 'DB*' (modified Davies-Bouldin index) [min]
%       - 'DN' (Dunn index) [max]
%       - 'DNG' (Dunn index using graphs) [max]
%       - 'DNS' (Modified Dunn index) [max]
%       - 'FOM' (Figure of merit)[min] clMethod | labelsDel
%       - 'GAMMA' (Gamma index)[min]
%       - 'GDI' (Generalized Dunn Indices (18))[max] GDI_interdist, 
%         GDI_intradist
%       - 'GPLUS' (G(+) index)[max]
%       - 'HOM' (Homogeneity (average))[max]
%       - 'HOMMIN' (Homogeneity (minimum))[max]
%       - 'I' (I-index or PBM)[max]
%       - 'SD' (SD index)[min] SD_alpha
%       - 'SDBW' (S_Dbw index)[min]
%       - 'SEP' (Separation (average))[min]
%       - 'SEPMAX' (Separation (maximum))[min]
%       - 'SF' (Score Function)[max]
%       - 'SIL' (Silhouette index (average))[max]
%       - 'SSI' (Simple Structure Index)[max]
%       - 'TAU' (Tau index)[min]
%       - 'VAR' (Variance index)[min]
%       - 'XB' (Xie-Beni index)[min]
%       - 'XB*' (Modified Xie-Beni index)[min]
%
%   unifyMethod
%       Which unification method to use:
%           'minmax'
%               max-like: value_i/max, min-like: min/value_i.
%           'range' 
%               Normalize on [0,1].
%           'range0'
%               Move range to start with 0.
%           'range1'
%               Move range to end with 1.
%           'pos'   
%               Move indeces with negative values to start with 0.
%           'prob'
%               Transform to probabilities, sum is 1.
%           'rank'
%               Rank values, 1 for the poorest value,..., ascending.
%           'rank10'
%               Rank values, 1 for the poorest value, other are less 
%               than 1 (suitable for rank aggregation).
%           empty or 'none'
%               Do not unify.
%
%   reduceMethod
%       Which method for feature reduction to use (to eliminate 
%       non-informative indices):
%
%       See pplk_featureReduce for more details.
%
%           'NONE'
%               Do not use any (default).
%
%           feature selection
%           'FSFS' 
%               Mitra et al., 2002.
%           'LS'   
%               Laplacian Score by He et al., 2005.
%           'SPEC' 
%               Spectral feature selection, Zhao & Liu 2007.
%           'FSKM'
%               Feature selection with k-medoids, Pepelka 2014.
%
%           feature extraction/transformation
%           'PCA'
%           'KPCA' 
%               Scholkopf et al., 1998.
%           'ICA'  
%               FastICA by Hyravinen et al., 2000.
%           any of supported methods in the Dimensionality reduction
%           toolbox (by van der Maaten).
%           Unsupervised (prefered):
%               'ProbPCA', 'MDS', 'FactorAnalysis',  'Isomap', 'Laplacian',
%               'HessianLLE', 'LTSA','FastMVU', 'DiffusionMaps', 'SNE',
%               'SymSNE', 'tSNE', 'SPE', 'Autoencoder'
%           Unsupervised (occasional singularity problem on PRM data)
%               'GPLVM', 'Sammon', 'LandmarkIsomap', 'LLE', 'MVU', 'CCA',
%               'LandmarkMVU', 'LPP', 'NPE', 'LLTSA', 'LLC',
%               'ManifoldChart', 'CFA'
%           Supervised / labeled
%               'LDA', 'GDA', 'NCA', 'MCML', 'LMNN'
%           'FEKM' 
%               Feature extraction with k-means, Pepelka 2014.
%
%   weightMethod, weigthMode
%       Which weighting method to use. You can choose among these:
%       
%       Legend:
%       Methods are in (), their possible modes in [].
%
%       - mean ('wMean')
%       - minimum ('wMin')
%       - maximum ('wMax')
%       - difference ('wDiff')
%       - Vega-Pons et al. 2010 ('wVegaPons')['CLK'|'CBK']
%       - rank aggregation ('wRankAggreg')
%         ['min'|'mean'|'median'|'geom.mean'|'stuart'|'RRA']
%
%   options.reduceDim 
%       Number of reduced dimensions or name of estimator (passed to 
%       pplk_featureReduce()).
%
%   options.CVImat
%       An ensembleSize-by-numIndices matrix with pre-computed indices 
%       values.  
%
%   options.dtype
%       Distance measure type:
%           1 - Euclidean distance
%           2 - Pearson correlation
%
%   options.NC 
%       Vector that contains numbers of clusters, e.g., [2 3 4 5]
%       length(NC) and number of labels sets in labels must match.
%
%	For other options, please see the help of pplk_validInt function.
%
%
% OUTPUTS
%   weigths
%       Vector of weights for each partition in ensemble.
%
%   PRM
%       An ensembleSize-by-numOfIndices Partition Relevance Matrix.
%
%   featInd
%       Indices of the remaining features, if selection step occured.
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka

callDir=chdir(pplk_homeDir());

reduceDim = [];

if ~exist('options','var')
    options = [];  
elseif isstruct(options)
    if isfield(options,'reduceDim')
        reduceDim = options.reduceDim;
        options = rmfield(options,'reduceDim');
    end
    if isfield(options,'CVImat')
        PRM_raw = options.CVImat;
        options = rmfield(options,'CVImat');
    end
else
    error('Variable options has to be struct.');
end

if ~exist('reduceMethod','var') || isempty(reduceMethod)
    reduceMethod = 'NONE';
end

if ~exist('weightMode','var') || isempty(weightMode)
    weightMode = '';
end

if ~exist('weightMethod','var') || isempty(weightMethod)
    %warning('pplk:PRM','No weight method specified, using none.');
    weightMethod = 'none';
end

if ~exist('unifyMethod','var') || isempty(unifyMethod)
    warning('pplk:PRM','No unification method specified, using range.');
    unifyMethod = 'range';
end

[N,ensembleSize] = size(labelsEns);

if ensembleSize < 2
    error('Ensemble size has to be at least 2!');
end

if ~exist('PRM_raw','var')
    nPoints = size(data,1);    
    if N ~= nPoints
        error('Number of data points inconsistent!');
    end
    % Run internal validity indices over ensemble of labels and obtain PRM matrix.
    [~, list] = pplk_validInt(data, labelsEns, indicesList, options);
    PRM_raw = list';
end

% Unification of values: scale on interval [0,1], where 0 represents the poorest
% and 1 the best quality; transforming into probabilities; ranking.

%	1. Determine, which indices needed to be reversed
listOfMinLike = {'APN', 'AD', 'ADM', 'CI', 'CON', 'DB', 'DB*','DBMOD', 'FOM', ...
    'GAMMA','SD','SDBW','SEP','SEPMAX','TAU','VAR','XB','XB*','XBMOD'};
mask = ismember(indicesList,listOfMinLike);

%	2. Unify indices
%   transform: min value means best result -> max value means best result
%   (note: this is not true for unify method 'rank01')
PRM = pplk_unifyPRM(PRM_raw, mask, unifyMethod);

%   3. Reduce features
%[PRM,featInd] = featureReduce(PRM, reduceMethod, 5);
% call: featureReduce(PRM, reduceMethod, k, varargin)
% if k is empty, MLE dimensionality estimator is used, otherwise k is a
% target dimensionality
[PRM_r,featInd] = pplk_featureReduce(PRM, reduceMethod, reduceDim);

% Check for monstrum columns with nans or infs
mask_nan_inf = any(isnan(PRM_r),1) | any(isinf(PRM_r),1);
if any(mask_nan_inf)
   PRM_r(:,mask_nan_inf) = [];
   featInd(mask_nan_inf) = [];
end
if size(PRM_r,2) == 0
   PRM_r = PRM; 
end

% if features were transformed (not selected), additional unification is
% necessary
if ~any(ismember(upper(reduceMethod),{'NONE','FSFS','SPEC','LS','FSKM'}))
    PRM_r = pplk_unifyPRM(PRM_r,[],unifyMethod);
end

%   4. Calculate weights
weights = pplk_weightPRM(PRM_r,weightMethod,weightMode);
% If all weights elements are zero, set them to 1.
% All-zeros may cause error in weighted consensus function.
if all(weights==0)
    weights = ones(size(weights));
end

chdir(callDir);
end
