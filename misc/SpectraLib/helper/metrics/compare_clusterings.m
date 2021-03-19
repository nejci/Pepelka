function distance = compare_clusterings( assig1, assig2 )

%function distance = compare_clusterings( assig1, assig2 )
%
% Computes the Variation of Information distance between two clusterings
% given by assig1,2. The values in assig1 must range between 1 and k1, 
% the number of clusters in assig1. Similarly for assig2.
%
%
%assig1, assig2(1,n) = vectors of integers from 1:k that represent the
%                      assignment to clusters. must have the same length
%

  k1 = max( assig1 );         % number clusters
  k2 = max( assig2 );
%  keyboard
%  size(assig1)
%  size(assig2)
  n = length( assig1 );          % nuber data points
  if size(assig1,2)==1 % single column mode
	assig1=assig1';
  end
  if size(assig2,2)==1 % single column mode
	assig2=assig2';
  end
%  keyboard;
  confusion = zeros( k1, k2 );     % build confusion matrix
  for i1 = 1:k1;
	for i2 = 1:k2;
%	  keyboard
	  confusion( i1, i2 ) = length( find( (assig1 == i1) & (assig2 == i2 )));
	  
	end;
  end;

  confusion = confusion/n;
%  disp (confusion)

  p1 = sum( confusion, 1 );
  p2 = sum( confusion, 2 );

  idummy = find( confusion > 0 );
  pdummy = p2*p1; % the P1(x)p2(x) matrix
  cdummy = confusion;
  cdummy( idummy ) = confusion( idummy )./pdummy( idummy );
  distance = -sum( p1.*locallog( p1 ))-sum( p2.*locallog( p2 ))-2*sum( sum( confusion.* locallog( cdummy )));
  distance = distance/log(2);




