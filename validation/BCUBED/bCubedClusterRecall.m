function B3CR = bCubedClusterRecall(refrPart, respPart)
% B-Cubed Recall cluster validity measure - averaged over clusters

score = 0;
clustRefrIDs = unique(refrPart);
nClustRefr = length(clustRefrIDs);

for clustID = 1:nClustRefr % for all clusters in reference partition
	cluster = (refrPart==clustRefrIDs(clustID));
	clusterInd = find(cluster);
	nElem = length(clusterInd);
	weight = (1/(nElem*nClustRefr));
	for elementID = 1:nElem
		element = clusterInd(elementID);
		score = score + weight * bCubedRecall(element,cluster,respPart);
	end
end

B3CR = score;

end
	