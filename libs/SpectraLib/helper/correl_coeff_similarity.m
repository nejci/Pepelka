function S=correl_coeff_similarity(xx) 
% function S=correl_coeff_similarity(xx) 
% Computes similarity between the points xx as the correlatino
% coefficient between them.  
  p=size(xx,1); 
  n=size(xx,2); 
	
  % first normalize all the vectors to have mean 0 and length 1
  for i=1:n 
	mean_i=mean(xx(:,i)); 
	xx(:,i)=xx(:,i)-mean_i; 
	length_i=norm(xx(:,i));
	xx(:,i)=xx(:,i)/length_i; 
  end
  
  % now the S_ij is nothing but dot(xx(:,i) ,xx(:,j)) 
  S=ones(n,n); % so S_ii is correct already.  
  for i=1:n-1 
	for j=(i+1):n
	  S(i,j) = dot(xx(:,i), xx(:,j)); 
	  S(j,i)=S(i,j); 
	end
  end
  
