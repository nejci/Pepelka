function [assignment_anchor]=cluster_point_anchor(xx, k,nmin) 
% function [assignment_anchor]=cluster_point_anchor(xx, k,nmin) 
% Clusters the given POINTS (Not similarity) using the anchor algorithm 
  
  p=size(xx,1); 
  n=size(xx,2);
  [ Anchor, nanchpoints, nanchors, ianchpoints ] = anchor(xx, k, nmin );
  kall = length( Anchor );
  if kall > k 
      fprintf('kall = %d (%s)\n',kall,sprintf('%d ',nanchpoints)); 
  end
 
  cluster_assignment=anchor_to_assignment(Anchor, nanchpoints, nanchors, ianchpoints,nmin,k, 2);
  
  % compute the centers of the first k clusters and then run kmeans 
  init_centers=zeros(p,k); 
  for i=1:k 
	ci=find(cluster_assignment==i); 
	if length(ci)==0 
	  warning('slight problem in anchor. Less than k clusters'); 
	  init_centers(:,i)=xx(:,i); 
	else
	  init_centers(:,i)=mean(xx(:,ci),2); 
	end
  end

  global KMEANS_THRESHOLD KMEANS_MAX_ITER;
  threshold = KMEANS_THRESHOLD;
  max_iter  = KMEANS_MAX_ITER; 
  [cluster_centers,assignment_anchor,distortion]=kmeans_nostat(max_iter,init_centers, threshold,xx);

  
