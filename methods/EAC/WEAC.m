function [labelsCons,numClust,weights] = WEAC(data,labelsEns,K,consFun,CVI,wMethod,options)
% WEAC
% [labelsCons,numClust,weights] = WEAC(data,labelsEns,K,consFun,CVI,wMethod,CVIoptions)
% Ensemble clustering by Evidence accumulation clustering (Fred & Jain, 2005)
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns   (matrix)	clustering ensemble; each COLUMN corresponds
%							to one clustering (standard in Pepelka).
%
%   K           (scalar)	number of clusters in output consensus
%	                        clustering; leave empty [] to automatically
%	                        determine this number by lifetime criterion.
%
%   consFun     (string)    which clustering metod to use to cluster
%                           co-occurence matrix:
%                             'single'    - nearest distance (default)
%                             'complete'  - furthest distance
%                             'average'   - unweighted average distance (UPGMA)
%                                             (also known as group average)
%                             'weighted'  - weighted average distance (WPGMA)
%                             'centroid'  - unweighted center of mass distance (UPGMC)
%                             'median'    - weighted center of mass distance (WPGMC)
%                             'ward'      - inner squared distance (min variance algorithm)
%               (cell)      combination of methods in cell string
%
%   CVI         (string)    cluster validity index identifier:
%   wMethod     (string)    'single': (SWEAC) 
%                           'joint':  (JWEAC)
%                           others supported by pplk_partitionRelevance
%   options     (struct)    structure of options that is passed to
%                           the function pplk_partitionRelevance()
%
% OUTPUTS
%   labelsCons  (vector)    data labels - consensus clustering
%               (matrix)    if CVI has more than one element and wMode is 
%                           'single', consensus is computed for every index
%                           or if consFun is a cell, consensus is computed
%                           for every method specified
%   numClust    (vector)    number of clusters in consensus clustering
%   weights     (vector)    weights assigned to each ensemble partition -
%                             value of normalized CVI(s)
%               (matrix)    if CVI has more than one element and wMode is 
%                           'single', weigths are computed for every index 
%--------------------------------------------------------------------------
% EXAMPLE
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2014  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 23-June-2014 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if ~exist('K','var')
    K=[]; % use lifetime cirterion to determine number of clusters
end

if ~exist('consFun','var') || isempty(consFun)
    consFun = 'single';
end

% If no CVI is specified this method becomes (unweighted) EAC
if ~exist('CVI','var')
    CVI = [];
end

if ~exist('wMethod','var') || isempty(wMethod)
    wMethod = 'single';
end

unifyMethod = 'minmax';
reduceMethod = 'none';
weightMode = [];

%--------------------------------------------------------------------------
[N,ensSize] = size(labelsEns);

weights = [];

if ~exist('options','var')
    options = [];
elseif isstruct(options)
    if isfield(options,'unifyMethod')
        unifyMethod = options.unifyMethod;
        options = rmfield(options,'unifyMethod');
    end
    if isfield(options,'reduceMethod')
        reduceMethod = options.reduceMethod;
        options = rmfield(options,'reduceMethod');
    end
    if isfield(options,'weightMode')
        weightMode = options.weightMode;
        options = rmfield(options,'weightMode');
    end
    if isfield(options,'weights')
        weights = options.weights;
        options = rmfield(options,'weights');
    end
end


if isempty(weights)
    if isempty(CVI)
        weights = ones(ensSize,1);
    else
        % Compute CVI values for ensemble labels        
        if strcmpi(wMethod,'single')
            [~,PRM] = ...
                pplk_partitionRelevance(data,labelsEns,CVI,unifyMethod,reduceMethod,'none',weightMode,options);
            weights = PRM;
            
        elseif strcmpi(wMethod,'joint')
            [~,PRM] = ...
                pplk_partitionRelevance(data,labelsEns,CVI,unifyMethod,reduceMethod,'none',weightMode,options);
            weights = mean(PRM,2);
        else
            [weights,~] = ...
                pplk_partitionRelevance(data,labelsEns,CVI,unifyMethod,reduceMethod,wMethod,weightMode,options);
        end
    end
end

% If wMethod is 'single' and CVI contains more than 1 index, compute
% consensus for every index
numWeightsSingle = size(weights,2);

% If there are more consensus functions specified
numConsFun = 1;
if iscell(consFun)
    numConsFun = length(consFun);
else
    consFun = {consFun};
end
labelsCons = zeros(N,numConsFun,numWeightsSingle);
numClust = zeros(1,numConsFun,numWeightsSingle);


% some parameters for clusterer
params = pplk_setParamsDefault();
params.HCL_distance = [];

for w = 1:numWeightsSingle
    % Compute co-occurence matrix, return vector of pairs
    CO = computeCO(labelsEns,[],'labels','vec',weights(:,w));
    
    for c = 1:numConsFun
        params.HCL_clustMethod = consFun{c};
        % Run linkage algorithm on CO vector/matrix to obtain consensus parititon.
        % If K is empty, number ob final clusters is automatically determined.
        [labelsCons(:,c,w), moreInfo] = pplk_runClusterer('HCL',CO, K, 1, params);
        numClust(1,c,w) = moreInfo{1}.numClust;
    end
end
