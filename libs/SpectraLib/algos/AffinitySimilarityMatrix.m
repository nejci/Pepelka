function S=AffinitySimilarityMatrix(points,sigma) 

% function S=DistanceSimilarityMatrix(points,sigma) 
% Returns the similarity between two points based on the distance between
% the two points. NOTE: The vectors are assumed to be COLUMN VECTORS 
%
% S_ij= exp(distance_ij^2/(2*sigma*sigma)); 
  
  exp_mult_factor=-0.5/(sigma*sigma); 
  num_vectors=size(points,2); 

  %evklidska razdalja - default od pdist
  S=squareform(pdist(points')); % need to transpose it pdist expects row vectors
  S=exp(exp_mult_factor*S);   
  
  
