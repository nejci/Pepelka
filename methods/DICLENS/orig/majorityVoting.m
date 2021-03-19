function finalClusters = majorityVoting(comp, numComp, clusters)

[nClust, N] = size(clusters);

% calculate votes
votes = zeros(numComp, N);

for compInd = 1:numComp
	metaCluster = (comp==compInd);
	votes(compInd,:)=sum(clusters(metaCluster,:),1);
end

% It depends on the order of rows in matrix votes when searching for max
% value. If there is a tie, we choose one at random!
p = randperm(numComp);
votes = votes(p,:); % to simulate behaviour of Java Mimaroglu's implementation

% assign data samples to the meta cluster with the highest vote rate
finalClusters = false(numComp,N);
[maxVotes, maxVotesInd]=max(votes,[],1);
finalClusters( sub2ind([numComp,N],maxVotesInd,1:N) ) = 1;

% remove empty clusters
count = sum(finalClusters,2);
finalClusters(count==0,:)=[];
