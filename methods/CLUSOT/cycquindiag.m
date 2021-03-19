function [x]=cycquindiag(a,b,c,d,e,S)
% [x]=cycquindiag(a,b,c,d,e,S)
%
% cycquindiag solves a system of linear equations with an periodic
% quindiagonal coefficient matrix A using the algorithm of Benson and
% Evans as described in "The Computer Journal", vol. 16, no. 3,
% 1973, p. 278/279.
%
% Inputs:
% a,b,c,d,e - (vectors) diagonal elements, e.g.
%             +-                                -+
%             | c1 d1 e1                    a1 b1|
%             | b2 c2 d2 e2                    a2|
%             | a3 b3 c3 d3 e3                   |
%             |    a2 ...                        |
%             |                             e_n-2|
%             | e_n-1     a_n-1 b_n-1 c_n-1 d_n-1|
%             | d_n e_n           a_n  b_n  c_n  |
%             +-                                -+
% S - (matrix) of right hand sides
%
% Output:
% X - (matrix) solution of A*X=S
%
% $Id: cycquindiag.m 601 2005-04-20 13:41:26Z dome $
% D. Brugger, 19 April 2005
% algo/cycquindiag.m

if(nargin == 0)
  test_cycquindiag();
  return;
end

% (1) reduce to upper triangular form
n=size(a,2);
g=a(2); a(2)=0; h=0;

for idx=2:n-4
  idxm1=idx-1;
  idxp1=idx+1;

  v=b(idx)/c(idxm1);
  c(idx)=c(idx)-v*d(idxm1);
  d(idx)=d(idx)-v*e(idxm1);
  b(idx)=g-v*b(idxm1);
  a(idx)=a(idx)-v*a(idxm1);
  S(idx,:)=S(idx,:)-v*S(idxm1,:);

  v=a(idxp1)/c(idxm1);
  b(idxp1)=b(idxp1)-v*d(idxm1);
  c(idxp1)=c(idxp1)-v*e(idxm1);
  a(idxp1)=-v*a(idxm1);
  g=-v*b(idxm1);
  S(idxp1,:)=S(idxp1,:)-v*S(idxm1,:);

  v=e(n-1)/c(idxm1);
  c(n-1)=c(n-1)-v*a(idxm1);
  d(n-1)=d(n-1)-v*b(idxm1);
  e(n-1)=h-v*d(idxm1);
  h=-v*e(idxm1);
  S(n-1,:)=S(n-1,:)-v*S(idxm1,:);

  v=d(n)/c(idxm1);
  d(n)=e(n)-v*d(idxm1);
  e(n)=-v*e(idxm1);
  b(n)=b(n)-v*a(idxm1);
  c(n)=c(n)-v*b(idxm1);
  S(n,:)=S(n,:)-v*S(idxm1,:);
end

e(n-3)=e(n-3)+a(n-3);
a(n-3)=g;
a(n-1)=a(n-1)+h;

idx=n-3;
idxm1=idx-1;
idxp1=idx+1;
idxp2=idx+2;
idxp3=idx+3;

v=b(idx)/c(idxm1);
c(idx)=c(idx)-v*d(idxm1);
S(idx,:)=S(idx,:)-v*S(idxm1,:);
d(idx)=d(idx)-v*e(idxm1);
e(idx)=e(idx)-v*a(idxm1);
a(idx)=a(idx)-v*b(idxm1);

v=a(idxp1)/c(idxm1);
b(idxp1)=b(idxp1)-v*d(idxm1);
c(idxp1)=c(idxp1)-v*e(idxm1);
S(idxp1,:)=S(idxp1,:)-v*S(idxm1,:);
d(idxp1)=d(idxp1)-v*a(idxm1);
e(idxp1)=e(idxp1)-v*b(idxm1);

v=e(idxp2)/c(idxm1);
a(idxp2)=a(idxp2)-v*d(idxm1);
b(idxp2)=b(idxp2)-v*e(idxm1);
S(idxp2,:)=S(idxp2,:)-v*S(idxm1,:);
c(idxp2)=c(idxp2)-v*a(idxm1);
d(idxp2)=d(idxp2)-v*b(idxm1);

v=d(idxp3)/c(idxm1);
e(idxp3)=e(idxp3)-v*d(idxm1);
a(idxp3)=a(idxp3)-v*e(idxm1);
b(idxp3)=b(idxp3)-v*a(idxm1);
c(idxp3)=c(idxp3)-v*b(idxm1);
S(idxp3,:)=S(idxp3,:)-v*S(idxm1,:);

idx=n-2;
idxm1=idx-1;
idxp1=idx+1;
idxp2=idx+2;

v=b(idx)/c(idxm1);
c(idx)=c(idx)-v*d(idxm1);
S(idx,:)=S(idx,:)-v*S(idxm1,:);
d(idx)=d(idx)-v*e(idxm1);
e(idx)=e(idx)-v*a(idxm1);

v=a(idxp1)/c(idxm1);
b(idxp1)=b(idxp1)-v*d(idxm1);
c(idxp1)=c(idxp1)-v*e(idxm1);
S(idxp1,:)=S(idxp1,:)-v*S(idxm1,:);
d(idxp1)=d(idxp1)-v*a(idxm1);

v=e(idxp2)/c(idxm1);
a(idxp2)=a(idxp2)-v*d(idxm1);
b(idxp2)=b(idxp2)-v*e(idxm1);
S(idxp2,:)=S(idxp2,:)-v*S(idxm1,:);
c(idxp2)=c(idxp2)-v*a(idxm1);

idx=n-1;
idxm1=idx-1;
idxp1=idx+1;

v=b(idx)/c(idxm1);
c(idx)=c(idx)-v*d(idxm1);
S(idx,:)=S(idx,:)-v*S(idxm1,:);
d(idx)=d(idx)-v*e(idxm1);

v=a(idxp1)/c(idxm1);
b(idxp1)=b(idxp1)-v*d(idxm1);
c(idxp1)=c(idxp1)-v*e(idxm1);
S(idxp1,:)=S(idxp1,:)-v*S(idxm1,:);

idx=n;
idxm1=n-1;
v=b(idx)/c(idxm1);
c(idx)=c(idx)-v*d(idxm1);
S(idx,:)=S(idx,:)-v*S(idxm1,:);

% (2) back substitution
x(n,:)=S(n,:)/c(n);
x(n-1,:)=(S(n-1,:)-d(n-1)*x(n,:))/c(n-1);
x(n-2,:)=(S(n-2,:)-d(n-2)*x(n-1,:)-e(n-2)*x(n,:))/c(n-2);
x(n-3,:)=(S(n-3,:)-d(n-3)*x(n-2,:)-e(n-3)*x(n-1,:)-a(n-3)*x(n,:))/c(n-3);
for idx=n-4:-1:1
  x(idx,:)=(S(idx,:)-d(idx)*x(idx+1,:) ...
            -e(idx)*x(idx+2,:)-a(idx)*x(n-1,:) ...
            -b(idx)*x(n,:))/c(idx);
end

function test_cycquindiag()

range=10:10:2500;
for k=1:size(range,2)
  [t1,t2]=mytimer(range(k));
  t_inv(k)=t1;
  t_cqd(k)=t2;
end

figure; plot(range,t_inv,'r',range,t_cqd,'g');
title('Lapack vs. Cqd');
legend('Lapack','Cqd');

fprintf('test_cycquindiag succeded\n');

function [t1,t2]=mytimer(N)
alpha=0.5; beta=1; gamma=0.2;
alpha = alpha* ones(1,N);
beta = beta*ones(1,N);

% produce the five diagonal vectors
alpham1 = [alpha(2:N) alpha(1)];
alphap1 = [alpha(N) alpha(1:N-1)];
betam1 = [beta(2:N) beta(1)];
betap1 = [beta(N) beta(1:N-1)];

a = betam1;
b = -alpha - 2*beta - 2*betam1;
c = alpha + alphap1 +betam1 + 4*beta + betap1;
d = -alphap1 - 2*beta - 2*betap1;
e = betap1;

% generate the parameters matrix
A = diag(a(1:N-2),-2) + diag(a(N-1:N),N-2);
A = A + diag(b(1:N-1),-1) + diag(b(N), N-1);
A = A + diag(c);
A = A + diag(d(1:N-1),1) + diag(d(N),-(N-1));
A = A + diag(e(1:N-2),2) + diag(e(N-1:N),-(N-2));
A = A + gamma * diag(ones(1,N));
tic;
eIA=inv(A);
t1=toc;

tic;
E=eye(N);
IA=cycquindiag([A(1,N-1) A(2,N) diag(A,-2)'], ...
               [A(1,N) diag(A,-1)'], ...
               diag(A), ...
               [diag(A,1)' A(N,1)], ...
               [diag(A,2)' A(N-1,1) A(N,2)], ...
               E);
t2=toc;
%IA
check_equal(eIA,IA,'eIA','IA',0.000000001);

