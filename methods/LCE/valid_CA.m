function ca = valid_CA(labels, truelabels)
%==========================================================================
% FUNCTION: ca = valid_CA(labels, truelabels)
% DESCRIPTION: A function for computing Classification Accuracy for a clustering result
%
% INPUTS: labels = cluster labels from a clustering result (N-by-1 vector)
%     truelabels = known cluster labels for each data points (N-by-1 vector)
%
% OUTPUT:     ca = Classification Accuracy score
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

nrow = length(truelabels);
C = unique(labels);
k = length(C); % the number of clusters
clusters = cell(1,k);
ca = 0; %classification accuracy
for i = 1:k
    ind = find(labels==C(i));
    clusters{i} = ind; %find data point members for the i-th cluster
    n = length(ind);
    for j = 1:n
        clusters{1,i}(j,2) = truelabels(ind(j)); %find the true labels for the j-th data point
    end
    TC = unique(clusters{1,i}(:,2)); 
    kTC = length(TC); %the number true classes in this cluster
    ind = [];
    for l = 1:kTC 
        ind(l) = length(find(clusters{1,i}(:,2)==TC(l))); %find #points in each cluster
    end
    [M,I] = max(ind);
    clusters{1,i}(:,3) = TC(I); %re-labeling with the magority class
    ca = ca + M; %sum the number of data points in majority class
end
ca = ca/nrow;
