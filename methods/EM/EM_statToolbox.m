function labels = EM_statToolbox(data,K,regularization,maxIter)
% Run MATLAB implementation of EM GMM algorithm.

options = statset();

if ~exist('regularization','var') || isempty(regularization)
    regularization = 0;
end
if exist('maxIter','var') && ~isempty(maxIter)
    options = statset('MaxIter',maxIter);
end

obj = gmdistribution.fit(data,K,'Regularize',regularization,'Options',options,'SharedCov',logical(1),'CovType','diagonal');
labels = cluster(obj,data);
