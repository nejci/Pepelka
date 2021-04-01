function [ elems1, elems2, lambda1, lambda2, ww1, ww2 ] = split_a_cluster( S, ww, nmin,optimization_option )

% function [ elems1, elems2, lambda1, lambda2, ww1, ww2 ] = split_a_cluster( S, ww, nmin )
%
% internal function for the ShiMalik clustering alg. splits a set of points into 
% 2 clusters by the min NCut heuristic
%
% INPUT: 
% S = a similarity matrix
% ww = 2nd eigenvector of the corresponding P matrix
% nmin > 0 if the size of a cluster is <= nmin, then lambda <-- -2
%	   and ww <--- ones( 1, nmin ) signalling that the cluster
%	   will not be split again		  
%
% OUTPUT
% elems1,2 = elements of the clusters
% lambda1,2 = for each cluster, a new P is computed and its second
%	       eigenvalue stored in lambda; if the cluster has 1 
%	       element only, lambda1,2 <-- -2
% ww1,2     = the 2nd eigenvector for each resulting cluster 
%	      (=1 if cluster has 1 elem)

  % warning off MATLAB:eigs:SigmaChangedToLA;
  n = length( ww );
  tt = sum( S );

  [ww iperm ] = sort( ww );
  oneton = 1:n;
  [ idummy invperm ] = sort( oneton( iperm )); % get inverse permutation

  ncuts = zeros(1,n-1);    % find best cut
  opts.disp=0;

  for icut = 1:n-1;
    i1 = 1:icut;
    i2 = (icut+1):n;
    i1 = iperm( i1 );
    i2 = iperm( i2 );
	switch (optimization_option) 
	  case {'ncut'}
	   mult_factor=(1/sum( tt( i1 ))+1/sum( tt( i2)));
	   ncuts( icut ) = sum( sum( S(i1, i2 )))*mult_factor;
	  case {'conductance'}
	   mult_factor=(1/min( sum(tt(i1)), sum(tt(i2)) ));
	   ncuts( icut ) = sum( sum( S(i1, i2 )))*mult_factor;
	 case {'gap'} % cut at the maximum gap. So need to "minimize" -gap=ww(icut)-ww(icut+1) 
	  ncuts(icut)=ww(icut)-ww(icut+1);
	 otherwise
	  warning('******** WRONG OPTIMIZATION_OPTION IN CLUSTER_SHI_R ********'); 
%	  mult_factor=1; 
	end
%    ncuts( icut ) = sum( sum( S(i1, i2 )))*mult_factor;
  end;

  [dummy imin ] = min( ncuts );
  i1 = 1:imin;
  i2 = imin+1:n;
  elems1 = iperm( i1 );
  elems2 = iperm( i2 );



  if( length( i1 ) > nmin )
	S1 = S( elems1, elems1 );
	T1 = diag( sum( S1 ));
	P1 = T1 \ S1;
	
	%Changed to use myeigs 
	[ ww1, lambda1 ] = myeigs( P1, 2); 
	
	%    [ ww1 lambda1]=eigsort(P1); 
	%	lambda1=diag(lambda1); 
	
	ww1 = ww1( :, 2 )';
	lambda1 = lambda1( 2 );


  else
	lambda1 = -2;
	ww1 = ones( 1, length( i1 ) );
  end;

  if( length( i2 ) > nmin )
	S2 = S( elems2, elems2 );
	T2 = diag( sum( S2 ));
	P2 = T2 \ S2;
	%Changed to use the eigsort as eigs may not be stable enough 
%	opts.disp=0;
	[ ww2, lambda2 ] = myeigs( P2, 2);
	
    % [ ww2 lambda2]=eigsort(P2); 
	% lambda2=diag(lambda2); 
	
	ww2 = ww2( :, 2 )';
	lambda2 = lambda2( 2 );
  else
	lambda2 = -2;
	ww2 = ones( 1, length( i2 ));
  end;

