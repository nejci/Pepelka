function viewSca(S,ca) 
% function viewSca(S,ca) 
% visualize the Similarity matrix and cluster assignment side by side. 

  figure(1); 
  clf; 
  imagesc(S); 
  colorbar; 
  figure(2); 
  clf; 
  S=S*0; 
  kmax=max(ca); 
  for i=1:kmax
	indices=find(ca==i); 
	S(indices,indices)=i; 
  end
  
  imagesc(S); 
  colorbar; 
