function [labels, moreInfo]=pplk_clustererEM(data,K,params)
%[labels, moreInfo]=pplk_clustererEM(data,K,params)
% Expectation-Maximization algoritem for Finite GAussian Mixtures
%
% INPUTS
%   data - input data matrix [N x d] of N samples, which are d-dimensional
%   K - desired number of clusters
%   params - parameters: maxIter, regularization
%   
% OUTPUTS
%   labels - vector of cluster assignments
%   moreInfo - running time
%

% Defaults
maxIter = 100;
regularization = 0;

if exist('params','var') && isstruct(params)
    if isfield(params,'EM_maxIter')
        maxIter = params.EM_maxIter;        
    end
    if isfield(params,'EM_regularization')
        regularization = params.EM_regularization;        
    end
end

runStatToolbox = 0;
if license('test','Statistics_Toolbox') == 1
    runStatToolbox = 1;
end

oldPath=chdir('EM');
ticID = tic();
if runStatToolbox
    labels = EM_statToolbox(data, K, regularization, maxIter);
else
    labels = emgm(data', K, regularization, maxIter);
end
time = toc(ticID);
chdir(oldPath);

labels=labels';
moreInfo.time=time;

end

