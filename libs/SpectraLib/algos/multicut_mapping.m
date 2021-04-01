function [vv]=multicut_mapping(S,k)
% function vv=multicut_mapping(S,k)
% Maps the similarity to the n points of dimension k-1 using the multicut
% lemma. 
  

  n=size(S,2);    

  D=diag(sum(S,2));
  % Compute the top k eigenvectors

  % use eigs instead of eig() 
  [ vv lambda] = myeigs_gen( S , D, k );

	vv = vv( :, 1:k ) ; %% Now vv is nxk
  vv=vv'; % makes it kxn i.e. n column vectors of dimension k
  
