function show_two_clusters(ca1,ca2,noclf)
% function show_two_clusters(ca1,ca2,noclf)
% visualizes the two cluseters side by side after "matching" the
% corresponding clusters. 

if size(ca1,1)==1
	ca1=ca1';
  end
  if size(ca2,1)==1
	ca2=ca2';
  end
  if nargin > 2 && noclf
  else
      clf;
  end
  
      
  % if they are the same size do the matching
  k=max(ca2);
  if max(ca1)==k
	[min_perm_error,permutation] = clustering_error( ca1,ca2)
	rev_perm(permutation)=1:k; 
	for i=1:length(ca2) 
	  ca2(i)=rev_perm(ca2(i)); 
	end
  end
  
  both_clusters=[ca1 ca2]; 
  imagesc(both_clusters); 
  colorbar; 
