function S = getWeightedMatrix(baseCls, I_ncai)
% Huang Dong. 14 May, 2014.
% Optimized for speed by Nejc Ilc, 10 June 2014

I_ncai = I_ncai(:);
w_sum = sum(I_ncai);

[n, nBC] = size(baseCls);
cntCol = max(baseCls);

S = zeros(n,n);
for k = 1:nBC
    tmp = double(bsxfun(@eq, 1:cntCol(k),baseCls(:,k)));
    idx = tmp * tmp';
    S = S + idx.*I_ncai(k);
end
S = S/w_sum;