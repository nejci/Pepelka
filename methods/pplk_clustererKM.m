function [labels, moreInfo]=pplk_clustererKM(data,K,params)

% defaults
maxIter = 100;
replicates = 1;
distance = 'sqEuclidean';

if exist('params','var') && isstruct(params)
    if isfield(params,'KM_maxIter')
        maxIter = params.KM_maxIter;        
    end    
    if isfield(params,'KM_nRuns')
        replicates = params.KM_nRuns;
    end
    if isfield(params,'KM_distance')
        distance = params.KM_distance;
    end
end


if ~ismember(distance, {'sqEuclidean', 'correlation'})
	fprintf('KM - Wrong distance: ''%s'', using ''sqEuclidean'' instead.\n', distance);
end
    
tic;
[labels,c,sumd] = kmeans(data,K, 'distance',distance,'start','sample',...
    'replicates', replicates, 'emptyaction','singleton','maxIter',maxIter,...
    'display','notify');
time=toc;

moreInfo.time=time;
moreInfo.distance=distance;