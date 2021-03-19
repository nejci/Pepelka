function [assignment_shi_r,param_string]=cluster_shi_r(S,k,optimization_option)

% clustering by the Shi & Malik recursive NCut algorithm,
% with two modifications:
% 1. the cluster with the largest lambda (2nd eigv of P) is split
% 2. the split is done by line search for the min ncut along the 
% sorted eigenvector 
% Note: modification 2. was suggested in the first
% version of the KVV paper where a similar alg was analyzed.
% (Therefore initially I called it the KVV algorithm.)
% The final version has a different algorithm, so I changed the name.
%
% Splits recursively until k clusters are obtained
%
% INPUTS
%    S = similarity matrix
%    k = final number of clusters < n PRECONDITION
%    nmin = min size cluster (a cluster <= nmin will not be split)
%   AS OF NOW the nmin is HARD WIRED TO HAVE THE VALUE 5
% LOCAL DATA 
%    KMAX size of internal tables
%    shi_lambda = 2nd eigenvalue for each [candidate] cut (equal to the NCut)
%    shi_n = number of elements in each cluster
%    shi_nclust = current number of clusters
%    shi_ww = weights in each [candidate] cluster
%    shi_elems = elements of each [candidate] cluster
%
% splits the cluster with the largest lambda. one of the resulting
% clusters ovewrites the parent clust, the other one becomes a new
% clust 

% compute the stochastic matrix. 
  
  T = sum( S );
  P = diag( T )\ S;
  n=length(S);
  nmin =1; % HARD CODED
  % Initialize
  KMAX=10; 
  shi_nclust = 1;
  shi_ww = zeros( KMAX, n );
  shi_elems = zeros( KMAX, n );
  shi_n = zeros( 1, KMAX );
  shi_lambda = zeros( 1, KMAX );
  assignment_shi_r = zeros( 1, n );

  % First split

  shi_elems( 1, : ) = 1:n;
  shi_n( 1 ) = n;

  isplit = 1; % index of split cluster
  %  [ vv lambda ] = eigsort( P );  % repeated from cluster_anchor
  [vv lambda]=myeigs(P,2); 
  shi_ww( 1, : ) = vv( :, 2 )';

  while ( shi_nclust < k )

	elems = shi_elems( isplit, 1:shi_n( isplit ));
	[ elems1, elems2, lambda1, lambda2, ww1, ww2 ] = split_a_cluster( S( elems, elems ), shi_ww( isplit, 1:shi_n( isplit)), nmin,optimization_option );

	nn1 = length( elems1 );
	nn2 = length( elems2 );
	shi_nclust = shi_nclust + 1;
	shi_elems( isplit, 1:nn1 ) = elems( elems1 );
	shi_elems( shi_nclust, 1:nn2 ) = elems( elems2 );
	shi_lambda( isplit ) = lambda1;
	shi_lambda( shi_nclust ) = lambda2;
	shi_ww( isplit, 1:nn1 ) = ww1;
	shi_ww( shi_nclust, 1:nn2 ) = ww2;
	shi_n( isplit ) = nn1;
	shi_n( shi_nclust ) = nn2;

	[ ldummy, isplit ] = max( shi_lambda( 1:shi_nclust ));
	if ( ldummy == -2 )
      disp( '** Error! cant find k clusters with these conditions' );
	end;

	%keyboard
  end; % while

  % Get assignments

  for iclust = 1:k;
    assignment_shi_r( shi_elems( iclust, 1:shi_n( iclust ))) = iclust;
  end;
  param_string=''; 
  
  