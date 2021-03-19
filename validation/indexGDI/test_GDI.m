% test of functions

%% Load data
[dataOrig, target] = pplk_loadData('real\iris_orig');
data = dataOrig;

%% cluster data
res = [];
Kmax = 12;
for K = 2:Kmax
    labels = pplk_runClusterer('AL',data,K,1);

    centers = zeros(K,size(data,2));
    for i = 1:K
        Dtmp = data(labels==i,:);
        centers(i,:) = mean(Dtmp,1);
    end
    
    distM = squareform(pdist(data));
    
    C2C = squareform(pdist(centers));
    
    % compute indices
    %GDI = indexGDI(data,labels,[3,2,1])';
    GDI = indexGDI(data,labels,[],[],'euclidean',[],'euclidean')';
    res = [res, GDI(:)];
    
end