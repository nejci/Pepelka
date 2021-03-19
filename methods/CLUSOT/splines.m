function [a,b,c,d]=splines(p,t)
% [a,b,c,d]=splines(p,t)
%
% splines computes a cubic spline s(t) defined by the coefficients
% a,b,c,d. The dimension of s(t) is determined by the dimension of
% points p. 
%
% Input:
% p - (matrix) points in R^d, e.g. for R^2 it is a (n x 2)-matrix.
%      Note: s(t) is a closed curve such that p(1) should equal p(n).
% t - (vector) knot vector where elements should be monotonically
%      increasing, e.g. t(1) <= t(2) <= ... <= t(n)
%
% Output:
% a,b,c,d - (matrices) coefficients defining s(t), e.g. for R^d
%             each will be a (n x d) matrix. To evaluate s(t) at
%             a specific point tt use 'spline_ipol(t,a,b,c,d,tt)'.
% 
% $Id: splines.m 479 2005-03-24 08:56:22Z dome $
% D. Brugger, 08 March 2005
% algo/splines.m

if(nargin == 0)
  test_splines();
  return;
end

[n,m]=size(p);
% first and last point are equal
n=n-1;
% preallocate
delta=zeros(n,1); dp=zeros(n,m);
v=zeros(n,m); A=zeros(n,n);
a=zeros(n,m); b=zeros(n,m); d=zeros(n,m);

for k=1:n
  % compute points differences ...
  dp(k,:)=p(k+1,:)-p(k,:);
  % ...and parameter differences
  delta(k)=t(k+1)-t(k);
  % set coefficients a
  a(k,:)=p(k,:);
end

% build up right hand side v and matrix A
v(1,:)=3*((dp(1,:)/delta(1))-(dp(n,:)/delta(n)));
A(1,1)=2*(delta(1)+delta(n)); A(1,2)=delta(1); A(1,n)=delta(n);
for l=2:n-1
  v(l,:)=3*((dp(l,:)/delta(l))-(dp(l-1,:)/delta(l-1)));
  A(l,l)=2*(delta(l)+delta(l-1));
  A(l,l-1)=delta(l-1);
  A(l,l+1)=delta(l);
end
v(n,:)=3*((dp(n,:)/delta(n))-(dp(n-1,:)/delta(n-1)));
A(n,1)=delta(n); A(n,n-1)=delta(n-1); A(n,n)=2*(delta(n)+delta(n-1));

% compute coefficients c
c=A\v;

% compute coefficients b and d
for l=1:n-1
  b(l,:)= dp(l,:)/delta(l) - (delta(l)/3 * (c(l+1,:)+ 2*c(l,:)));
  d(l,:)=(c(l+1,:)-c(l,:))/(3*delta(l));
end
b(n,:)= dp(n,:)/delta(n) - (delta(n)/3 * (c(1,:)+ 2*c(n,:)));
d(n,:)=(c(1,:)-c(n,:))/(3*delta(n));

function test_splines()
% Test case #1 - a circle
p=[1 0; ...
   0 1; ...
   -1 0; ...
   0 -1; ...
   1 0];
t=[0 pi/2 pi (3/2)*pi 2*pi];
tt=[0:2*pi/255:2*pi];
%t=[0 1 2 3 4];
%tt=[0:0.1:4];
ea=[1 0; ...
    0 1; ...
    -1 0; ...
    0 -1];
eb=[0 1.5; ...
    -1.5 0;...
    0 -1.5; ...
    1.5 0];
ec=[-1.5 0; ...
    0 -1.5; ...
    1.5 0; ...
    0 1.5];
ed=[0.5 -0.5; ...
    0.5 0.5; ...
    -0.5 0.5; ...
    -0.5 -0.5];
[a,b,c,d]=test_splines_helper(p,t,tt,'Test case #1 - a circle')
check_equal(a,ea,'a','ea');
check_equal(a,ea,'b','eb');
check_equal(a,ea,'c','ec');
check_equal(a,ea,'d','ed');

% Test case #2 - ellipse
p=[2 0; ...
   0 1; ...
   -2 0; ...
   0 -1; ...
   2 0];
test_splines_helper(p,t,tt,'Test case #2 - ellipse')

% Test case #3 - irregular ellipse
p=[2 0; ...
   0 1; ...
   -1.5 0; ...
   0 -2.5; ...
   2 0];
test_splines_helper(p,t,tt,'Test case #3 - irregular ellipse')

% Test case #4 - more than 4 points
p=[2 0; ...
   1 1.5; ...
   0 2.5; ...
   -0.5 2; ...
   -1 3; ...
   -0.5 -0.5; ...
   0 -2; ...
   2 0];
t=[0 1 2 3 4 5 6 7];
tt=[0:0.001:7];
test_splines_helper(p,t,tt,'Test case #4 - more than 4 points')

% Test case #5 - just 2 points
p=[2 0; ...
   -2 0; ...
   2 0];
t=[0 1 2];
tt=[0:0.01:2];
test_splines_helper(p,t,tt,'Test case #5 - just 2 points')

% Test case #6 - 3 points but only in one halfplane
p=[2 0.5; ...
   1 1.5; ...
   -1.5 1; ...
   2 0.5];
t=[0 1 2 3];
tt=[0:0.01:3];
test_splines_helper(p,t,tt,'Test case #6 - 3 points but only in one halfplane')

% Test case #7 - curve in 3d space
p=[1 0 0; ...
   0 1 0; ...
   0 0 1; ...
   0 -1 0; ...
   0 -1 -2; ...
   0 0.5 0.5; ...
   1 0 0];
t=[0 1 2 3 4 5 6];
tt=[0:0.01:6];
test_splines_helper(p,t,tt,'Test case #7 - curve in 3d space',1);

% Test case #8 - corner neuron in grid
p=[0 0
   0 1; ...
   1 0; ...
   0 0];
t=[0 1 2 3];
tt=[0:0.01:3];
test_splines_helper(p,t,tt,'Test case #8 - corner neuron in grid');

% Test case #9 - edge neuron in grid
p=[0 0; ...
   1 1; ...
   0 2; ...
   0 0];
t=[0 1 2 3];
tt=[0:0.01:3];
test_splines_helper(p,t,tt,'Test case #9 - edge neuron in grid',0,[0 1]);

% Test case #10
p=[3 2.5; ...
   2 4; ...
   1 3.5; ...
   2 2; ...
   3 2.5];
t=[0.4636, 1.5708, 2.1588, 4.4528, 6.7468];
tt=[0.4636:0.01:6.7468];
%t=[0 1 2 3 4];
%tt=[0:0.01:4];
test_splines_helper(p,t,tt,'Test case #10',0,[2 2]);

% Test case #11
p=[3 2.5; ...
     2 4; ...
     2 0; ...
   3 2.5];
t=[0.4636 1.5708 4.7124 6.7468];
tt=[0.4636:0.01:6.7468];
%t=[0 1 2 3];
%tt=[0:0.01:3];
test_splines_helper(p,t,tt,'Test case #11',0,[2 ...
                    2]);

np=[0 0];
p=[2 0; ...
   0 1; ...
   -2 0; ...
   0 -1; ...
   2 0];
t=[0 pi/2 pi 3/2*pi 2*pi];
tt=[0:0.01:2*pi];
test_splines_helper(p,t,tt,'Test case #12',0,np);

function [a,b,c,d]=test_splines_helper(p,t,tt,str,flag,nij)
if(nargin < 5)
  flag=0; % flag will be set for 3d curves
end
if(nargin < 6)
  nij=[0 0];
end
[a,b,c,d]=splines(p,t)
val=spline_ipol(t,a,b,c,d,tt);
%for k=1:size(tt,1)
%  val(k,:)=spline_ipol(t,a,b,c,d,tt(k));
%end
figure;
if(flag)
  plot3(val(:,1),val(:,2),val(:,3));
else
  plot(val(:,1),val(:,2));
end
hold on;
for k=1:size(p,1)
  if(flag)
    plot3(p(k,1),p(k,2),p(k,3),'rx','MarkerSize',10);
  else
    plot(p(k,1),p(k,2),'rx','MarkerSize',10);
  end
end
if(flag)
  plot3(0,0,0,'go','MarkerSize',15);
  xlabel('x');ylabel('y');zlabel('z');
else
  plot(nij(1),nij(2),'go','MarkerSize',15);
  text(nij(1)-eps,nij(2)-eps,'N_{ij}');
  xlabel('x');ylabel('y');
end
grid on;
title(str);
