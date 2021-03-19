function [val]=spline_ipol(t,a,b,c,d,tt)
% [val]=spline_ipol(t,a,b,c,d,tt)
%
% spline_ipol computes the value of s(t) as defined
% by coefficients a,b,c,d. spline_ipol is very fast when tt
% is a vector or a matrix as opposed to calling spline_ipol in a loop.
%
% Inputs:
% t - (vector) knot vector where elements should be monotonically
%      increasing, e.g. t(1) <= t(2) <= ... <= t(n)
% a,b,c,d - (matrices) coefficients defining s(t), e.g. for R^d
%             each will be a (n x d) matrix. 
% tt - (vector|matrix) points at which s(t) will be evaluated
%
% Output:
% val - (vector|matrix) the values of s(t) at points tt
% 
% $Id: spline_ipol.m 447 2005-03-17 13:25:03Z dome $
% D. Brugger, 08 March 2005
% algo/spline_ipol.m

lo=NaN; hi=NaN;
n=size(tt,1);
m=size(tt,2);
u=zeros(n,m);
for k=1:n
  for l=1:m
    if(~isnan(lo) && ~isnan(hi) && ...
       t(lo) <= tt(k,l) && t(hi) > tt(k,l))
      u(k,l)=lo;
      % nothing to do, use cached values
    else
      % find interval using bisection
      lo=1; hi=size(t,2);
      while(hi-lo > 1)
        idx=floor((hi+lo)/2);
        if(t(idx) > tt(k,l))
          hi=idx;
        else
          lo=idx;
        end
      end
      % tt(k,l) is now in [t(lo),t(hi)),
      % store index value in u
      u(k,l)=lo;
    end
  end
end
% compute values s(t) at points tt
if(n == 1) 
  % tt is a vector
  tmp1=repmat(tt-t(u),size(a,2),n)';
  %tmp2=tmp1.*tmp1; tmp3=tmp2.*tmp1;
  size(tmp1)
  size(a(u,:))
  val=a(u,:)+tmp1.*(b(u,:)+tmp1.*(c(u,:)+tmp1.*d(u,:)));
  %au=a(u,:); bu=b(u,:); cu=c(u,:); du=d(u,:);
else
  % tt is a matrix
  tmp1=tt-t(u); %tmp2=tmp1.*tmp1; %tmp3=tmp2.*tmp1;
  a1=a(:,1); b1=b(:,1);  c1=c(:,1);  d1=d(:,1); 
  a2=a(:,2); b2=b(:,2);  c2=c(:,2);  d2=d(:,2); 
  val(:,:,1)=a1(u)+tmp1.*(b1(u)+tmp1.*(c1(u)+tmp1.*d1(u)));
  val(:,:,2)=a2(u)+tmp1.*(b2(u)+tmp1.*(c2(u)+tmp1.*d2(u)));
%  val(:,:,1)=a1(u)+b1(u).*tmp1+c1(u).*tmp2+d1(u).*tmp3;
%  val(:,:,2)=a2(u)+b2(u).*tmp1+c2(u).*tmp2+d2(u).*tmp3;
end


