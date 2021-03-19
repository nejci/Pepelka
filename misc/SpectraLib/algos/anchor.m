function [ Anchor, nanchpoints, nanchors, ianchpoints ] =  anchor(xx, k, nmin )

% function [ Anchor, nanchpoints, nanchors, ianchpoints  ] =  anchor( p, n, xx, k, nmin )
%
% takes as input a p dimensional data set xx with n points
% and k the number of anchors. produces a list of anchors
% plus the assignments of points to them in anchor
%
% Anchor = list of indices of the anchors
% nanchors = how many are real anchors (have >= nmin points)
% nanchpoints = number of points for each anchor
% ianchpoints = list of points for each anchor
%
% MAXANCHORS = 50;  total nr of anchors it tests
%
% nanchors is either =k or less if the limit MAXANCHORS is reached 
% before k real anchors are found

MAXANCHORS = 25;  % total nr of points it tests

% xx(p,n)
p=size(xx,1);
n=size(xx,2); 


Anchor = zeros( 1, MAXANCHORS );  % the anchors
ianchpoints = zeros( n-1, MAXANCHORS );
                                  % the point list for each anchor
distanch = zeros( n-1, MAXANCHORS );
                                  % the distances to the anchor 
                                  % (sorted in increasing order)
nanchpoints = zeros( 1, MAXANCHORS ); 
                                  % nr points in each list

nanchors = 0;

for ia = 1:MAXANCHORS;

% First anchor

  if( ia == 1 )   % First anchor
     [ yy imax ] = max( xx(1,:) );  
     Anchor( 1 ) = imax;
     ianchpoints( :, 1 ) = [ 1:imax-1, imax+1:n ]';
     nanchpoints( 1 ) = n-1;
     xxtemp = xx( :, ianchpoints(:,1) );
     dtemp = sqrt( sum((xxtemp - repmat( xx( :, Anchor(1)), [ 1, n-1 ])).^2, 1));
     [ dtemp2 itemp ] = sort( dtemp );
     distanch( :, 1 ) = dtemp2';
     ianchpoints( :, 1 ) = ianchpoints( itemp, 1 ); 
  else

% Find new anchor

     [ dtemp imax ] = max( max( distanch ));
     Anchor( ia ) = ianchpoints( nanchpoints( imax ), imax );
     ilist = [];      % temporary storage for new anchor - points
     na = 0;          % number points
     da = [];          % distances
     xa = xx( :, Anchor( ia ));


%  delete new anchor from old list

     ianchpoints( nanchpoints( imax ), imax ) = 0; 
     distanch( nanchpoints( imax ), imax ) = 0;
     nanchpoints( imax ) = nanchpoints( imax )-1;

% Find points that belong to that anchor

     for iolda = 1:ia-1   

%  steal points from anchor iolda

       d = sqrt( sum( (xx( :, Anchor( iolda )) - xx( :, Anchor( ia ))).^2, 1));
       icandidates = ianchpoints( find( distanch( :, iolda ) > d/2 ), iolda );
       ncandidates = length( icandidates );
       xxtemp = xx( :, icandidates );
       ntrunc = nanchpoints( iolda ) - ncandidates;  % number points of iaa
                                                    % definitely not stolen
         
       dtemp = sqrt( sum( (xxtemp - repmat( xa, [1, ncandidates ] )).^2, 1));
       dtempold = distanch( ntrunc+1:nanchpoints( iolda ), iolda )';
       isteal = find( dtemp < dtempold );
       istay = find( dtemp >= dtempold );
       
% give stolen points to new anchor ia

       ilist = [ ilist; icandidates( isteal ) ];
       na = na + length( isteal );
       da = [ da dtemp( isteal ) ];

% give other points back to old anchor iolda

       nold = nanchpoints( iolda );
       nnew = ntrunc + length( istay );
       distanch( ntrunc+1:nnew, iolda ) = dtempold( istay )';
       distanch( nnew+1:nold, iolda ) = 0;
       ianchpoints( ntrunc+1:nnew, iolda ) = icandidates( istay );
       ianchpoints( nnew+1:nold, iolda ) = 0;
       nanchpoints( iolda ) = nnew;
  
     end; % for 

% put data of new anchor in place

     nanchpoints( ia ) = na;
     [ dtemp2 itemp ] = sort( da );
     ianchpoints( 1:na, ia ) = ilist( itemp );
     distanch( 1:na, ia ) = dtemp2';
  end % if new anchor

% see how many real anchors 

  nanchors = length( find( nanchpoints >= nmin ));
  if nanchors == k;
     break;
  end;
end; % for

Anchor = Anchor( 1: ia );
ianchpoints = ianchpoints( :, 1:ia );
nanchpoints = nanchpoints( 1:ia );
distanch = distanch( :, 1:ia );


