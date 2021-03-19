function plot2Dpoints_with_clusters(points2D, cluster_assignment,color_vector)  

% function plot2Dpoints_with_clusters(points2D, cluster_assignments,color_vector)  
% Plots the the 2D points with the colors chosen according the
% cluster_assignments and the color_vector. (The default is chosen if the
% color_vector are not given. 
% The 2Dpoints has first ROW as X and second ROW as Y 
% cluster assignment start from 1.

  n=size(points2D,2);
  if n>0
	if (nargin < 3)
	  color_vector=['kh'; 'ro'; 'bx'; 'g+'; 'ms'; 'c^'; 'yv'; 'kd'; 'k>'; 'k<' ];
	end 
	
	hold_state=ishold;
	hold on
	for i=1:n  
	  if cluster_assignment(i)<7
		fstring=color_vector(cluster_assignment(i)+1,:);
	  else
		fstring='k*';
	  end
	  plot(points2D(1,i),points2D(2,i),fstring);
	  
	end
	
	% restore the hold state as before
	if (~hold_state)
	  hold off 
	end
	
  end
  