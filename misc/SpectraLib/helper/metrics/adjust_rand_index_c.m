function ari=adjust_rand_index_c(C) 
% function ari=adjusted_rand_index(C) 
%
% Calculates the adjusted rand index given the confusion matrix or otherwise
% called the contingency matrix. 
% NOTE : This expects the matrix to have the actual counts (integers and
% not the ratio).   
  
  
  leni=size(C,1); 
  lenj=size(C,2); 
  % sumi(i)= n_{i.}. 
  % sumi(j)= n_{.j}. 
  
  sumi=sum(C,2); 
  sumj=sum(C,1); 
  
  maxlen=max(max(sumi),max(sumj)); 
  
  nc2=zeros(1,maxlen); 
  for i=1:maxlen-1 
	nc2(i+1)=nc2(i)+i; 
  end
  
  t1=0;  % t1=sum_ij (n_ij C 2)
  t2=0; % t2=sum_i (n_i C 2)
  t3=0; % t3=sum_j (n_j C 2)
  
  for i=1:leni
	if sumi(i)>0
	  t2=t2+nc2(sumi(i)); 
	end
	
	for j=1:lenj
	  if C(i,j) > 0 
		t1=t1+nc2(C(i,j)); 
	  end
	end
  end
  
  for j=1:lenj
	if sumj(j)>0
	  t3=t3+nc2(sumj(j));
	end
  end
  n=sum(sum(C)); % sum_ij n_ij 
  t4=t2*t3 / (n*(n-1)/2);
  ari=(t1-t4)/( (0.5*(t2+t3)) - t4 );
  
	
  
  