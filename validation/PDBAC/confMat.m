function [C,gn] = confMat(target,predict)

target = target(:);
predict = predict(:);

gLen = size(target,1);
[idx,gn] =grp2idx([target;predict]);
gidx = idx(1:gLen);
ghatidx = idx(gLen+1:gLen*2);

%ignore NaN values in GIDX and GHATIDX
nanrows = isnan(gidx) | isnan(ghatidx);
if any(nanrows)
    gidx(nanrows) = [];
    ghatidx(nanrows) = [];
end

gLen = size(gidx,1);
gnLen =length(gn);

C = accumarray([gidx,ghatidx], ones(gLen,1),[gnLen, gnLen]);