function [C]=vector_edge_fill(C,B,null,fill)
% [C]=vector_edge_fill(C,B,null,fill)
%
% vector_edge_fill fills a polygon given by border pixels in B
% in the matrix C using fill. Background pixels are assumed to
% have pixel value null.
%
% This algorithm is the one described in Dunlavey83, with the
% exception, that border pixels are additionally set after the
% algorithm has terminated. Thus is doesn't matter if points in B
% are given in clock- or counterclockwise order -- border pixels
% are always drawn. (TODO: Might want to generalize it - see notes
% in paper)
%
% Inputs: 
% C - (matrix) "screen" where pixels are set
% B - (matrix) border pixels (n x 2)
% null - value of background pixels
% fill - fill value for polygon
% 
% $Id: vector_edge_fill.m 474 2005-03-23 09:25:13Z dome $
% D. Brugger, 02 March 2005
% algo/vector_edge_fill.m

if(nargin == 0)
  test_vector_edge_fill();
  return;
end

%[n,m]=size(C);
n=size(C,1);

% init q, q(y)=NaN means, q(y) is empty
q=ones(n,1)*NaN;

% translate border pixels to first quadrant, so we can use them for
% indexing matrix C
[x_min y_min]=min(B);
if(x_min < 1)
  B(:,1)=B(:,1)+abs(x_min)+1;
end
if(y_min < 1)
  B(:,2)=B(:,2)+abs(y_min)+1;
end

% treat special case - just one point
if(size(B,1) == 1)
  C(B(1,2),B(1,1))=fill; return;
end

x=B(1,1); y=B(1,2); 
[firstclass,dx,dy]=mov_class(B(1,:),B(2,:));
x=x+dx; y=y+dy;
oldclass=firstclass;
k=2;
while(k < size(B,1))
  [newclass,dx,dy]=mov_class(B(k,:),B(k+1,:));
  [q,C]=take_action(x,y,oldclass,newclass,C,null,fill,q);
  x = x+dx; y = y+dy;
  oldclass = newclass;
  k = k + 1;
end

if(k <= size(B,1))
  [newclass,dx,dy]=mov_class(B(k,:),B(1,:));
  [q,C]=take_action(x,y,oldclass,newclass,C,null,fill,q);
  x = x+dx; y = y+dy;
  oldclass = newclass;
  [q,C]=take_action(x,y,oldclass,firstclass,C,null,fill,q);
end

% set border pixels
for k=1:size(B,1)
  C(B(k,2),B(k,1))=fill;
end

function [q,C]=take_action(x,y,oldclass,newclass,C,null,fill,q)
% movement table
persistent table;
% build table on first call
if(isempty(table)) 
  % Coding: 
  %  movement: NorthEast=1, SouthEast=2, East=3, West=4,
  %  action:  Clear=1, Clear1=2, Invert=3, Nop=4
  table = [ 1 3 4 1; ...
            3 2 2 4; ...
            1 4 4 1; ...
            4 2 2 4 ...
          ];
end
switch table(oldclass,newclass)
 case 1 % Clear
  [q,C]=clear(C,null,fill,x,y,q);
 case 2 % Clear1
  [q,C]=clear(C,null,fill,x-1,y,q);
 case 3 % Invert
  C=invert_pixel(x,y,C,null,fill);
 case 4
  ; % Nop
 otherwise
  error('!!! Unkown action code %d !!!', table(oldclass,newclass))
end

function [d,dx,dy]=mov_class(p1,p2)
dx=p2(1)-p1(1); dy=p2(2)-p1(2);
if(dy == 1)
  d=1; % North
  return; 
end
if(dy == -1)
  d=2; % South
  return; 
end
if(dx > 0)
  d=3; return; % East
end
if(dx < 0)
  d=4; return; % West
end

error(['!!! Why did we get here? dx=%g dy=%g p1=(%g,%g) p2=(%g,%g)!!' ...
       '!'],dx,dy, p1(1), p1(2), p2(1), p2(2))

function C=invert_pixel(x,y,C,null,fill)
if(C(y,x) == null)
  C(y,x) = fill;
else 
  if(C(y,x) == fill)
    C(y,x) = null;
  end
end

function C=invert_scan_segment(x1,x2,y,C,null,fill)
if(x1 < x2)
  ind_null=find(C(y,x1+1:x2) == null);
  ind_fill=find(C(y,x1+1:x2) == fill);
  C(y,x1+ind_null)=fill; C(y,x1+ind_fill)=null;
%  for k=x1+1:x2
%    C=invert_pixel(k,y,C,null,fill);
%  end
else
  ind_null=find(C(y,x2+1:x1) == null);
  ind_fill=find(C(y,x2+1:x1) == fill);
  C(y,x2+ind_null)=fill; C(y,x2+ind_fill)=null;
%  for k=x2+1:x1
%    C=invert_pixel(k,y,C,null,fill);
%  end
end

function [q,C]=clear(C,null,fill,x,y,q)
%C=invert_scan_segment(0,x,y,C,null,fill);
if(isnan(q(y)))
  q(y) = x;
else
  C=invert_scan_segment(q(y),x,y,C,null,fill);
  q(y) = NaN;
end

function test_vector_edge_fill()
null=0; fill=1; 
% Test case #1
C=zeros(8,8);
eC=[0 0 0 0 0 0 0 0; ...
    1 1 1 0 0 1 1 0; ...
    1 1 1 0 1 1 1 1; ...
    0 1 1 1 1 1 1 0; ...
    0 0 1 1 1 1 1 0; ...
    0 0 1 1 1 1 1 0; ...
    0 0 1 1 1 1 0 0; ...
    0 0 0 1 1 0 0 0];
B=[1 2; ...
   1 3; ...
   2 4; ...
   3 5; ...
   3 6; ...
   3 7; ...
   4 8; ...
   5 8; ...
   6 7; ...
   7 6; ...
   7 5; ...
   7 4; ...
   8 3; ...
   7 2; ...
   6 2; ...
   5 3; ...
   4 4; ...
   3 3; ...
   3 2; ...
   2 2];
%B=flipdim(B,1);
C=vector_edge_fill(C,B,null,fill)
check_equal(eC,C,'eC','C');

% Test case #2
C=zeros(9,8);
eC=[0 0 0 0 0 0 0 0; ...
   0 0 0 0 1 0 0 0; ...
   0 0 0 1 1 1 0 0; ...
   0 0 1 1 1 1 1 0; ...
   0 1 1 1 1 1 1 1; ...
   0 0 1 1 1 1 1 0; ...
   0 0 0 1 1 1 0 0; ...
   0 0 0 0 1 0 0 0; ...
   0 0 0 0 0 0 0 0];
B=[5 2; ...
   4 3; ...
   3 4; ...
   2 5; ...
   3 6; ...
   4 7; ...
   5 8; ...
   6 7; ...
   7 6; ...
   8 5; ...
   7 4; ...
   6 3];
%B=flipdim(B,1);
C=vector_edge_fill(C,B,null,fill)
check_equal(eC,C,'eC','C');

% Test case #3
eC=[0 2 2 2 0; ...
    1 0 2 0 1; ...
    1 1 1 1 1; ...
    1 1 1 1 0; ...
    0 1 1 0 0; ...
    0 0 1 1 0; ...
    0 0 0 0 0];
C=[0 2 2 2 0; ...
   0 0 2 0 0; ...
   0 0 0 0 0; ...
   0 0 0 0 0; ...
   0 0 0 0 0; ...
   0 0 0 0 0; ...
   0 0 0 0 0];
B=[1 2; ...
   1 3; ...
   1 4; ...
   2 5; ...
   3 6; ...
   4 6; ...
   3 5; ...
   4 4; ...
   5 3; ...
   5 2; ...
   4 3; ...
   3 3; ...
   2 3];
C=vector_edge_fill(C,B,null,fill)
check_equal(eC,C,'eC','C');

% Test case #4
C=zeros(7,7);
eC=[0 0 1 0 0 0 0; ...
    0 1 1 1 1 1 0; ...
    1 1 1 1 1 1 1; ...
    1 1 1 1 1 1 1; ...
    1 1 1 1 1 1 1; ...
    0 1 1 1 1 1 0; ...
    0 0 1 0 1 0 0]
B=[3 1; ...
   2 2; ...
   1 3
   1 4
   1 5
   2 6
   3 7
   4 6
   5 7
   6 6
   7 5
   7 4
   7 3
   6 2
   5 2
   4 2
   4 3
   4 4
   4 3
   4 2];
C=vector_edge_fill(C,B,null,fill)
check_equal(eC,C,'eC','C');
fprintf('test_vector_edge_fill succeded\n')

   
   
   