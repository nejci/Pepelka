function [max_val,matching]=max_bipartite_matching(C)
% function [max_val,matching]=max_bipartite_matching(C)
% Given the Confusion matrix uses the mam bipartite matching to match
% cluster numbers in the two clusters assignment to minimize clustering
% error. 
  % C= m x n element matching. 
  % m <= n 
	
	m=size(C,1); 
	n=size(C,2); 
	
	
	Clin=zeros(m*n,1); 
	xtest=Clin; 
	
	for i=1:m
	  for j=1:n
		Clin((i-1)*n+j)=-C(i,j);
	  end
	  xtest((i-1)*n+i)=1;
	end
	
	A=zeros(m+n,m*n); 
	% for all i Sum_j A_ij <= 1  
	for i=1:m
	  A(i,(i-1)*n+[1:n])=1; 
	end

	% for all j Sum_i A_ij <= 1  
	for j=1:n
	  A(j+m,[0:m-1]*n+j)=1; 
	end

	b=ones(m+n,1);
	ub=ones(m*n,1);
	lb=ub-1;
 
	Ylin=linprog(Clin, A ,b,[],[],lb,ub,[],	optimset('Display','off'));
	max_val=-Clin'*Ylin;
	for i=1:m
	  temp=Ylin((i-1)*n+[1:n])';
	  [dummy matching(i)]=max(temp);
	end
	
