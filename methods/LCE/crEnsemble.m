function E = crEnsemble(X, M, k, scheme)
%==========================================================================
% FUNCTION: E = crEnsemble(X, M, k, scheme)
% DESCRIPTION: This function generates a cluster ensemble from the dataset X
%              using K-means algorithm
%
% INPUTS:   X = a dataset, rows of X correspond to observations; columns
%               correspond to variables (exclude class labels!!)
%           M = the prefered number of base clusterings in the ensemble
%           k = the prefered number of clusters in the base clusterings
%      scheme = cluster ensemble generating scheme (1 = Fixed k, 2 = Random k)
%
% OUTPUT:   E = matrix of cluster ensemble
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

E = [];
if scheme == 1 %Fixed k
    for i = 1: M
        Cl = kmeans(X,k,'emptyaction','singleton');
        while length(unique(Cl)) ~= k;
            Cl = kmeans(X,k,'emptyaction','singleton');
        end
        E = [E Cl];
    end
else %Random K
    for i = 1: M
        K = floor((k-2+1)*rand+2); % create randomized number of clusters in [2,k]
        Cl = kmeans(X,K,'emptyaction','singleton');
        while length(unique(Cl)) ~= K;
            Cl = kmeans(X,K,'emptyaction','singleton');
        end
        E = [E Cl];      
    end
end