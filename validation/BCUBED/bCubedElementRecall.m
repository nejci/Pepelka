function B3CR = bCubedElementRecall(refrPart, respPart)
% B-Cubed Recall cluster validity measure - averaged over elements

score = 0;
clustRefrIDs = unique(refrPart);
nClustRefr = length(clustRefrIDs);
weight = (1/length(refrPart));

for clustID = 1:nClustRefr % for all clusters in reference partition
	cluster = (refrPart==clustRefrIDs(clustID));
	clusterInd = find(cluster);
	nElem = length(clusterInd);
	
	for elementID = 1:nElem
		element = clusterInd(elementID);
		score = score + weight * bCubedRecall(element,cluster,respPart);
	end
end

B3CR = score;

end
	