function SACT = makeWeightsBtwCls_SACT(bcs, bcsSegs, I_ncai)
% Huang Dong. 19 April, 2013.
% Optimized for speed by Nejc Ilc, 10 June 2014

nCls = size(bcsSegs,1);

maxBcs = max(bcs);
minBcs = min(bcs);

bcsSegs_numel = sum(bcsSegs,2);

W_Y = eye(nCls);
nClsIdx = find(bcsSegs_numel);
for i = nClsIdx(:)'
    j = i+1:nCls;
    j = j(bcsSegs_numel(j)~=0);
    interCnt = sum(bsxfun(@and, bcsSegs(i,:), bcsSegs(j,:)),2);
    W_Y(i,j) = interCnt ./ (bcsSegs_numel(i)+bcsSegs_numel(j)-interCnt);
end


W_Y = W_Y + W_Y';
W_Y(1:nCls+1:nCls^2) = W_Y(1:nCls+1:nCls^2) - 2 + (bcsSegs_numel ~= 0)';
W_Y_binary = (W_Y>0);

nClsVec = (1:nCls)';
allWhichBcs = bsxfun(@le,nClsVec,maxBcs) & bsxfun(@ge,nClsVec,minBcs); 
I_ncaiBcs = allWhichBcs*I_ncai;

SACT_sum = zeros(size(W_Y));
for i = 1:nCls
    j = i+1:nCls;
    commonNeighbors = bsxfun(@and, W_Y_binary(:,i),W_Y_binary(:,j));
    weightedTerm = bsxfun(@times,commonNeighbors,I_ncaiBcs);    
    SACT_sum(i,j) = sum(weightedTerm.* bsxfun(@min,W_Y(:,i), W_Y(:,j)),1);
end
diagWeight = sum(bsxfun(@times,W_Y_binary,I_ncaiBcs) .* W_Y);
SACT_sum = SACT_sum + SACT_sum';
SACT_sum(1:nCls+1:nCls^2) = SACT_sum(1:nCls+1:nCls^2) + diagWeight;
SACT = SACT_sum / max(SACT_sum(:));
SACT(1:nCls+1:nCls^2) = 1;