function mncut=compute_mncut(P,k,pi,ca) 
% function mncut=compute_mncut(P,k,pi,ca) 
%
% Computes the generalized MNCut of the cluster assignment ca given the Prob matrix
% P, number of clusters k and the stationary distribution pi. 
  
  ptotal_self=0; 
  for i=1:k
	ci=find( ca == i) ;
    if (length(ci) > 0) 
	  pic=pi(ci); 
	  den=sum(pic); 
	  pc=P(ci,ci); 
	  pcRowSum=sum(pc,2); 
	  num=sum(pcRowSum.*pic); 
	  if (num > den) 
		warning('compute_mncut: num > den');
	  end
	  
	  if (den>0)
		ptotal_self= ptotal_self+(num/den); 
	  end
	end
  end
  mncut=k-ptotal_self; 
  
