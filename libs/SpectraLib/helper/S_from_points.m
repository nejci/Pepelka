function S=S_from_points(points,sigma,smoothing_constant,self_similarity_zero)
% function S=S_from_points(points,sigma,smoothing_constant,self_similarity_zero)
% Calculates the simlilarity as the affinity. 
% smoothing_constant : if nonzero the the similarity is smoothed. 
% self_similarity_zero : if 1 make the self similarity zero 
 
  % calculate the similarity in case sigma > 0 
  if sigma > 0 
  	similarity=AffinitySimilarityMatrix(points,sigma) ; 
  end

  % smooth the similarity in case the smoothing_constant > 0 
  if smoothing_constant > 0 
	similarity=smooth_similarity(similarity,smoothing_constant); 
  end 
  
  if nargin > 3
	if self_similarity_zero > 0 
	  similarity=makeDiagonalZero(similarity); 
	end
  end
  S=similarity; 
