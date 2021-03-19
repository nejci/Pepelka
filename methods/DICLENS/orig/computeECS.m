function [sumECS] = computeECS(c1,c2,origClusters,nOrigClusterings)

[nClust,N] = size(origClusters);
 
sumECS = 0;
if ((sum(c1) > 0) && (sum(c2) > 0))

	c1_OR_c2 = c1 | c2;
	c1_AND_c2 =  c1 & c2;
	c1_ANDCMPL_c2 = c1 & ~c2;
	CMPL_c1_AND_c2 = ~c1 & c2;

	ECS = 0;

	for i=1:nClust

		dummy = c1_OR_c2 & origClusters(i,:);
		count = sum(dummy);

		dummy = c1_ANDCMPL_c2 & origClusters(i,:);
		count2 = sum(dummy);

		dummy = CMPL_c1_AND_c2 & origClusters(i,:);
		count3 = sum(dummy);

		dummy = c1_AND_c2 & origClusters(i,:);
		count4 = sum(dummy);

		ECS = ECS + (count * (count - 1) / 2 - count2 * (count2 - 1) / 2 - ...
			        count3 * (count3 - 1) / 2 + count4 * (count4 - 1) / 2);
	end
	
	ECS = ECS + sum(c1_AND_c2) * nOrigClusterings;
	ECS = ECS / (sum(c1) * sum(c2));

	sumECS = sumECS + ECS;
end