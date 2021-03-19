function [weights] = pplk_weightPRM(PRM,weightMethod,weightMode)

callDir=chdir(pplk_homeDir());

if ~exist('weightMethod','var') || isempty(weightMethod)
    weightMethod = 'none';
end

if ~exist('weightMode','var') || isempty(weightMode)
    weightMode = '';
end


% According to unified PRM, calculate weigths of each ensemble member.
switch weightMethod
    case {'none','single'}
        % ensembleSize = size(PRM,1);
        % weights = ones(ensembleSize,1);        
        weights = PRM;        
        
    case 'wMin'
        weights = wMin(PRM);
    case 'wMin2'
        weights = min(PRM,[],2);
        
    case 'wMax'
        weights = wMax(PRM);
    case 'wMax2'
        weights = max(PRM,[],2);
        
    case 'wMean'
        weights = wMean(PRM);
    case {'wMean2','joint'}
        weights = mean(PRM,2);
        
    case 'wDiff'
        weights = wDiff(PRM);
    case 'wDiff2'
        weights = wDiff2(PRM);
        
    case 'wVegaPons'
        switch weightMode
            case {'','CLK'}
                weights = wVegaPonsCLK(PRM);
            case 'CBK'
                weights = wVegaPonsCBK(PRM);
            otherwise
                error('Wrong mode for weight method wVegaPons! Select CLK or CBK.');
        end
        
    case 'wVegaPons2'
        switch weightMode
            case {'','CLK'}
                weights = wVegaPonsCLK2(PRM);
            case 'CBK'
                weights = wVegaPonsCBK2(PRM);
            otherwise
                error('Wrong mode for weight method wVegaPons! Select CLK or CBK.');
        end
        
    case 'wRankAggreg'
        weights = wRankAggreg(PRM,weightMode);
        
    case 'wRankAggreg2'
        weights = wRankAggreg2(PRM,weightMode);
        
    otherwise
        error('Wrong weight method!');
end

chdir(callDir);
end

% =========================================================================

function weights = wMean(PRM)
% Mean index value for particular clustering
	weights = mean(PRM,2);
    
    % transform weights so that their sum is 1
    weights = weights ./ sum(weights);
end

function weights = wMin(PRM)
% ratio of min index value of a particular clustering to sum of min values 
% of all clusterings 
% NOTE: Be careful when PRM values are normalized to [0, X] (range, range0).
	minRow = min(PRM,[],2);
    minRowSum = sum(minRow);
    if(minRowSum == 0)
        warning('pplk:divBy0','Division by zero; sum of min indices is zero! Weights are set to 0.');
        weights = zeros(size(PRM,1),1);
    else
        weights = minRow ./ minRowSum;
    end
    
end

function weights = wMax(PRM)
% ratio of max index value to sum of max values for particular clustering
	weights = max(PRM,[],2)./ sum(max(PRM,[],2));
end

function weights = wDiff(PRM)
% ratio of difference of min and max index values to sum of differences - 
% prefer small ones, which mean some kind of stability across various indices.
% Mirror diff values to prefer smaller.

    minRow = min(PRM,[],2);
    maxRow = max(PRM,[],2);
    
    % mirror diff over mid-point
    tmp = maxRow-minRow;
    if all(tmp==0)
        warning('pplk:divBy0','Division by zero; sum of diff between min and max is zero! Weights are set to 1/numEl.');
        weights = ones(size(PRM,1),1)/size(PRM,1);
    else
         mid = (min(tmp,[],1) + max(tmp,[],1)) /2;
         tmp = tmp + 2 * bsxfun(@minus, mid, tmp);
         
         % calculate weights on mirrored data
         weights = tmp ./ sum(tmp);
    end
end

function weights = wDiff2(PRM)
% ratio of difference of min and max index values to sum of differences - 
% prefer small ones, which mean some kind of stability across various indices.
% Mirror diff values to prefer smaller.

    minRow = min(PRM,[],2);
    maxRow = max(PRM,[],2);
    
    % mirror diff over mid-point
    tmp = maxRow-minRow;
    if all(tmp==0)
        warning('pplk:divBy0','Division by zero; sum of diff between min and max is zero! Weights are set to 1.');
        weights = ones(size(PRM,1),1);
    else
         mid = (min(tmp,[],1) + max(tmp,[],1)) /2;
         tmp = tmp + 2 * bsxfun(@minus, mid, tmp);
         
         % calculate weights on mirrored data
         weights = tmp;
    end
end

% =========================================================================
function weights = wVegaPonsCLK(PRM)
% Clustering with lack of knowledge:
% "This situation is when the user does not know how to validate the results. 
% The user has a cluster analysis problem and decides to use a clustering 
% ensemble technique to face it, but he is not completely sure about the 
% properties that should be considered as good in its problem. 
% As we do not have any information about what partition is more appropriate 
% for the specific problem in the cluster ensemble, we will obey the decision 
% of the majority. It means that, we assign a small weight to a partition 
% which behaves very different to the rest of partitions, since it is probably 
% a noisy partition obtained by a not appropriate generation mechanism. 
% Otherwise, if a partition has an average behavior with respect to the PVIs, 
% it will have a high weight assigned. In other words, as we are trying to 
% find the median partition, we highlight the average behavior."
%
% PRM is a matrix [ensembleSize X numIndices] with normalized values on
% [0,1] (use unify with method range or prob)
%
% Reference: 
% Vega-Pons, S., Correa-Morris, J., & Ruiz-Shulcloper, J. (2010). 
% Weighted partition consensus via kernels. 
% Pattern Recognition, 43(8), 2712–2724

    if ~(all(max(PRM,[],1) <= 1) && all(min(PRM,[],1) >= 0))
        error('Index values are not within interval [0,1]!');
    end

    [m,l]=size(PRM);
    
    % 1. calculate probabilities
    A = sum(PRM,1);
    PRMp = bsxfun(@rdivide, PRM, A);

    % 2. calculate entropy
    H = PRMp .* log2(PRMp);
    H(isnan(H)) = 0;
    H = -sum(H,1);

    % 3. calculate column means
    mu = A./m;

    % 4. calculate weights
    weights = sum(bsxfun(@times,H,(1 - abs(bsxfun(@minus,PRM,mu)))),2);
    
    % 5. transform weights so that their sum is 1
    weights = weights ./ sum(weights);
        
end

function weights = wVegaPonsCLK2(PRM)
% DO NOT divide weights by sum of weights.
    if ~(all(max(PRM,[],1) <= 1) && all(min(PRM,[],1) >= 0))
        error('Index values are not within interval [0,1]!');
    end

    [m,l]=size(PRM);
    
    % 1. calculate probabilities
    A = sum(PRM,1);
    PRMp = bsxfun(@rdivide, PRM, A);

    % 2. calculate entropy
    H = PRMp .* log2(PRMp);
    H(isnan(H)) = 0;
    H = -sum(H,1);

    % 3. calculate column means
    mu = A./m;

    % 4. calculate weights
    weights = sum(bsxfun(@times,H,(1 - abs(bsxfun(@minus,PRM,mu)))),2);        
end

function weights = wVegaPonsCBK(PRM)
% Clustering with background knowledge:
% "This situation is when the user knows what characteristics in
% the results are considered good. The user has an index or a set of indexes
% that he wants to maximize.1 In this case, the weight assigned to each
% partition is a measure of how close the index values of a partition were
% to the maximum value of each index. We highlight the partitions that
% maximize the set of indexes by giving to them a high weight. That way,
% partitions with a behavior similar to the expected by the user, will have
% a high weight value assigned."
%
% PRM is a matrix [ensembleSize X numIndices] with normalized values on [0,1]
%
% Reference: 
% Vega-Pons, S., Correa-Morris, J., & Ruiz-Shulcloper, J. (2010). 
% Weighted partition consensus via kernels. 
% Pattern Recognition, 43(8), 2712–2724.

    if ~(all(max(PRM,[],1) <= 1) && all(min(PRM,[],1) >= 0))
        error('Index values are not within interval [0,1]!');
    end

    % 1. Calculate maximum index values in for each column (index)
    M=max(PRM,[],1);

    % 2. Calculate weights, which are sum of differences between index values
    % and maximum index value
    weights = sum(1 - abs(bsxfun(@minus,PRM,M)),2);
    
    % 3. transform weights so that their sum is 1
    weights = weights ./ sum(weights);

end

function weights = wVegaPonsCBK2(PRM)
% DO NOT divide weights by sum of weights.

    if ~(all(max(PRM,[],1) <= 1) && all(min(PRM,[],1) >= 0))
        error('Index values are not within interval [0,1]!');
    end

    % 1. Calculate maximum index values in for each column (index)
    M=max(PRM,[],1);

    % 2. Calculate weights, which are sum of differences between index values
    % and maximum index value
    weights = sum(1 - abs(bsxfun(@minus,PRM,M)),2);
end

function weights = wRankAggreg(PRM, weightMode)
    oldPath = chdir(['..',filesep,'methods',filesep,'PRM']);
    % always unify with rankRankAggreg (1=worst, min=best)
    % PRM should not be unified by rankRankAggreg before! Any other method is
    % allowed.
    PRM = pplk_unifyPRM(PRM,[],'rankRankAggreg');
    weights = aggregateRanks(PRM,[],weightMode,1);
    chdir(oldPath);
    
    % mirror over mid-point
    tmp=weights;
    mid = (min(tmp,[],1) + max(tmp,[],1)) /2;
    tmp = tmp + 2 * (mid-tmp);
        
    % transform weights so that their sum is 1
    weights = tmp ./ sum(tmp);
end

function weights = wRankAggreg2(PRM, weightMode)
% DO NOT divide weights by sum of weights.
    oldPath = chdir(['..',filesep,'methods',filesep,'PRM']);
    % always unify with rankRankAggreg (1=worst, min=best)
    % PRM should not be unified by rankRankAggreg before! Any other method is
    % allowed.
    PRM = pplk_unifyPRM(PRM,[],'rankRankAggreg');
    weights = aggregateRanks(PRM,[],weightMode,1);
    chdir(oldPath);
    
    % mirror over mid-point
    tmp=weights;
    mid = (min(tmp,[],1) + max(tmp,[],1)) /2;
    weights = tmp + 2 * (mid-tmp);
end

