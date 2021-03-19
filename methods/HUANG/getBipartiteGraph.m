function B = getBipartiteGraph(bcs, bcsSegs, I_ncai, alpha)
% Huang Dong. 14 May, 2014.
% Optimized for speed by Nejc Ilc, 5 June 2014

% n:    the number of data points.
% nBase:    the number of base clusterings.
% nCls:     the total number of clusters (in all base clusterings).

[nCls,n] = size(bcsSegs);

W_Y = makeWeightsBtwCls_SACT(bcs, bcsSegs, I_ncai);

W_XY = zeros(nCls, n);
sim = I_ncai'.*alpha;
bcsIdx = bsxfun(@plus, bcs, (0:n-1)'*nCls);
bcsIdx = reshape(bcsIdx',numel(bcsIdx),1);
W_XY(bcsIdx) = repmat(sim,1,n);
B = [W_XY';W_Y];