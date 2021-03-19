function [C,match,cost] = getcmClust(target,predict)
% Compute confusion matrix for clustering - align clusters and classes with
% hungarian algorithm.

if ~exist('predict','var') || isempty(predict)
    isConfusion = 1;
else
    isConfusion = 0;    
end

% 1. compute confusion matrix
if ~isConfusion
    target = target(:);
    predict = predict(:);
    C1 = getcm(target,predict);
else
    C1 = target;
end

% 2. run hungarian algorithm to align clusters with classes
[match,cost] = hungarian(-C1);
C = C1*match';