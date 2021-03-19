function ICS = computeICS(c, clustersAll)
% ICS - Intra Cluster Similarity for a cluster c

[nClust, N] = size(clustersAll);
total = 0;

for ind = 1:nClust
	cluster = clustersAll(ind,:);
		
	dummy = c & cluster;
	s = sum(dummy);	
	total = total + s*(s-1)/2;
end

s = sum(c);
denom = s*(s-1)/2;

if denom == 0
	ICS = 0;
else
	ICS = total/denom;
end