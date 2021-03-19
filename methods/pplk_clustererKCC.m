function [labels, moreInfo]=pplk_clustererKCC(data,K,params)
%K-centers clustering aka. K-medoids aka. PAM (Partitioning Around Medoids)

tic
[labels, energy, medoid] = kmedoids(data',K);
time=toc;

moreInfo.time=time;



function [label, energy, medoid] = kmedoids(X,m)
% X: d x n data matrix
% m: k (1 x 1) or label (1 x n, 1<=label(i)<=k) or medoid index(1 x k)
% Written by Mo Chen (mochen@ie.cuhk.edu.hk). March 2009.

% initialization
n = size(X,2);
s = length(m);
if s == 1
    k = m;
    medoid = randsample(n,k);
elseif s < n
    k = s;
    medoid = m;
elseif s == n
    k = max(m);
    label = m;
    medoid = zeros(k,1);
else
    error('ERROR: m is not valid.');
end

D = sqdistance(X);
if s ~= n
    [val,label] = min(D(medoid,:));
end

% main algorithm
last = 0;
while any(label ~= last)
    for i = 1:k
        idx = (label==i);
        [~,tmp] = min(sum(D(idx,idx),1));
        idx = find(idx);
        medoid(i) = idx(tmp);
    end  
    last = label;
    [val,label] = min(D(medoid,:));
end
energy = sum(val);



function D = sqdistance(A, B, M)
% Square Euclidean or Mahalanobis distances between all sample pairs
% A: d x n1 data matrix
% B: d x n2 data matrix
% M: d x d  Mahalanobis matrix
% D: n1 x n2 pairwise square distance matrix
% Written by Michael Chen (sth4nth@gmail.com). July 2009.

if nargin == 1
    A = bsxfun(@minus,A,mean(A,2));
    S = full(sum(A.^2,1));
    D = full((-2)*(A'*A));
    D = bsxfun(@plus,D,S);
    D = bsxfun(@plus,D,S');
elseif nargin == 2
    assert(size(A,1)==size(B,1));
    
    m = (sum(A,2)+sum(B,2))/(size(A,2)+size(B,2));
    A = bsxfun(@minus,A,m);
    B = bsxfun(@minus,B,m);
    D = full((-2)*(A'*B));
    D = bsxfun(@plus,D,full(sum(B.^2,1)));
    D = bsxfun(@plus,D,full(sum(A.^2,1)'));
elseif nargin == 3
    assert(size(A,1)==size(B,1));
    
    m = (sum(A,2)+sum(B,2))/(size(A,2)+size(B,2));
    A = bsxfun(@minus,A,m);
    B = bsxfun(@minus,B,m);
    D = full((-2)*(A'*M*B));
    D = bsxfun(@plus,D,full(sum(B.*(M*B),1)));
    D = bsxfun(@plus,D,full(sum(A.*(M*A),1)'));
end
