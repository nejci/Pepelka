function wi=wallas_index(trueC,C) 
% function wi=wallas_index(trueC,C) 
% Computes the one-side Wallas index of C w.r.t to true C 
% is 1 if clustering is same and zero if disagrees on all the points.   
  kmax=max(trueC); 
  N11=0; 
  denom=0; 
  
  for k=1:kmax
	cluster_index=find(trueC==k); 
	clen=length(cluster_index);
	for i=1:(clen-1)
	  for j=(i+1):clen
		denom=denom+1; 
		if C(cluster_index(i))==C(cluster_index(j))
		  N11=N11+1; 
		end
	  end
	end
  end
  
  wi=N11/denom; 
  
  
	
  
