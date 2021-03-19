function norm_list = normalize_minmax(list)
% min-max normalization
minV = min(list);
maxV = max(list);

% In a case minV equals maxV, there will occur division by zero. To prevent
% this set normalized values to zero.
if minV == maxV
    norm_list = zeros(size(list));
else
    norm_list = (list - minV) / (maxV - minV);
end