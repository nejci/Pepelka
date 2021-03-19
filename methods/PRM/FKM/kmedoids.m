function [label, energy, index] = kmedoids(X,k,maxIters)
% X: d x n data matrix
% k: number of cluster
% maxIters: maximum number of iterations, default 100
% Written by Mo Chen (sth4nth@gamil.com)
if ~exist('maxIters','var') || isempty(maxIters)
   maxIters = 100; 
end

v = dot(X,X,1);
D = bsxfun(@plus,v,v')-2*(X'*X);
n = size(X,2);
[~, label] = min(D(randsample(n,k),:),[],1);
last = 0;
iters = 1;
while any(label ~= last) && (iters <= maxIters)
    [~, index] = min(D*sparse(1:n,label,1,n,k,n),[],1);
    last = label;
    [val, label] = min(D(index,:),[],1);
    iters = iters + 1;
end
energy = sum(val);
