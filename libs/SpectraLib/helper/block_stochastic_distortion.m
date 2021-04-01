function [distortion, dist_matrix,mean_matrix]=block_stochastic_distortion(P,ca)
% function distortion=block_stochastic_distortion(P,ca)
%
% calcualtes  the distortion of P from being a block stochastic matrix
% w.r.t the cluster assignment ca. dist_matrix returns the distortion
% matrix for each subblock P_ij and distortion is the total . 
% NOTE : The distortion that is returned is absolute and the the
% distortion in the dist_matrix returned is relative. 
  
  kmax=max(ca); 
  for i=1:kmax
	cluster{i}=find(ca==i); 
  end
  
  A=zeros(kmax,kmax); 
  distortion=0;
  for i=1:kmax
	if length(cluster{i})
	  for j=1:kmax
		if length(cluster{j})  % non trivial sub block 
		  Pij=P(cluster{i},cluster{j});
		  sumPij=sum(Pij,2); 
		  varPij=var(sumPij);
		  A(i,j) = sqrt(varPij)/mean(sumPij); 
           mean2_matrix(i,j) = mean(sumPij); 
		  distortion = distortion + varPij; 
		end
	  end
	end
  end
  distortion=sqrt(distortion); 
  if nargout>1
	dist_matrix=A; 
  end 
  if nargout>2
      mean_matrix=mean2_matrix;
  end
    