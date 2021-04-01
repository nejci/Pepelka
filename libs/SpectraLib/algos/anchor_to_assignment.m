function cluster_assignment=anchor_to_assignment(Anchor, nanchpoints, nanchors, ianchpoints,nmin,k,  reassign_option)

% function cluster_assignment=anchor_to_assignment(Anchor, nanchpoints, nanchors, ianchpoints,reassign_option)
%
% To be used in conjuction with 'anchor' 
% Anchor = list of indices of the anchors
% nanchors = how many are real anchors (have >= nmin points)
% nanchpoints = number of points for each anchor
% ianchpoints = list of points for each anchor
% nmin = min number of points that are present in each "real" cluster.   
% k= number of "real" clusters. 
% reassign_option = the option which decides how the the cluster numbers
%                are assigned. 
% 0 :  assign as it is. 
% 1 :  make all the points which are assigned to "unreal" anchors as
%      outliers and the rest from one to k. 
% 2 :  renumber the main clusters from 1 to k and rest from k+1 onwards.   
% 3:  ****TODO*** reassign the "unreal" points to real points. 
  
  
kall=length(Anchor);
n=sum(nanchpoints)+nanchors;
cluster_assignment = zeros(n, 1 );
if reassign_option==0 
  % compute assignments 
  for ia = 1:kall;
	cluster_assignment( ianchpoints( 1:nanchpoints( ia ), ia )) = ia;
	cluster_assignment( Anchor( ia )) = ia;
  end;
end

if reassign_option==1
  k_counter=1; 
  for ia = 1:kall;
	if nanchpoints(ia)>=nmin
	  current_k=k_counter;
	  k_counter=k_counter+1;
	else
	  current_k=0; % outlier 
	end
	
	cluster_assignment( ianchpoints( 1:nanchpoints( ia ), ia )) = current_k;
	cluster_assignment( Anchor( ia )) = current_k;
  end;
  
end

if reassign_option==2
  real_k_counter=1;
  number_of_real_clusters= length( find( nanchpoints >= nmin ));
  unreal_k_counter=number_of_real_clusters+1; 
  for ia = 1:kall;
	if nanchpoints(ia)>=nmin
	  current_k=real_k_counter;
	  real_k_counter=real_k_counter+1;
	else
	  current_k=unreal_k_counter;
	  unreal_k_counter=unreal_k_counter+1;
	end
	
	cluster_assignment( ianchpoints( 1:nanchpoints( ia ), ia )) = current_k;
	cluster_assignment( Anchor( ia )) = current_k;
  end;
  
end
