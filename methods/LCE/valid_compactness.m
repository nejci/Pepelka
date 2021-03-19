function cp = valid_compactness(X, labels) 
%==========================================================================
% FUNCTION: cp = valid_compactness(X, labels) 
% DESCRIPTION: A function for computing Compactness validity index 
%              for a clustering result
%
% INPUTS:  X = a dataset, rows of X correspond to observations; columns
%              correspond to variables (exclude class labels!!)
%     labels = cluster labels from a clustering result (N-by-1 vector)
%
% OUTPUT: cp = compactness score
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

n = length(labels); % number of data points
C = unique(labels); %all clusters
k = length(C); %number of clusters
cp = 0; %initialize compactness
for i=1:k %for each cluster
    ind = find(labels==C(i)); %find data point members for the i-th cluster
    nk = length(ind);
    if nk <= 1 %singleton cluster
        cp = cp + 0;
    else
        sum_d = 0;
        sum_d = sum(pdist(X(ind,:)));
        cp = cp + (nk*(sum_d/(nk*(nk-1)/2)));
    end
end
cp = cp/n;
