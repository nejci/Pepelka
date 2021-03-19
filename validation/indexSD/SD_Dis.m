function [Dis] = SD_Dis(data,labels)

[~,D] = size(data);
[clusterIDs, ~, lbl] = unique(labels);
c = length(clusterIDs);

% compute centroids of clusters
v = zeros(c,D);
for i = 1:c
    C = data(lbl==i,:);
    v(i,:) = mean(C,1);    
end

% compute total separation between clusters Dis
dist_v = sqrt(sqdistance2(v,v));

Dmax = max(dist_v(:));
Dmin = min(dist_v(dist_v>0));
if isempty(Dmin)
    Dmin = 0;
end

Dis = 0;
for k=1:c    
    z=1:c;
    z(k)=[];
    Dis = Dis + 1 / sum(dist_v(k,z));
end
Dis = Dmax/Dmin * Dis;
