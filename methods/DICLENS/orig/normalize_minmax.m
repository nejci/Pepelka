function norm_list = normalize_minmax(list)
% min-max normalization on interval [0,1]
minV = min(list);
maxV = max(list);

norm_list = (list - minV) / (maxV - minV);