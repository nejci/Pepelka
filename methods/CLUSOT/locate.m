function [u]=locate(t,tt)
% [u]=locate(t,tt)
%
% locate finds the index vector u such that t(u(k)) <= tt(k) <
% t(u(k)) using bisection. Worst case running time for |t|=n,
% |tt|=m is thus O(m log_2 n).
%     
% Inputs:
% t - (vector) vector representing intervals, thus 
%       elements in t should be monotonically
%       increasing, e.g. t(1) <= t(2) <= ... <= t(n)
% tt - (vector) points for which the interval needs to be determined.
%
% Output:
% u - (vector) index vector, see description above.
% 
% $Id: locate.m 417 2005-03-09 16:52:45Z dome $
% D. Brugger, 09 March 2005
% algo/locate.m

if(nargin == 0)
  test_locate();
  return
end

lo=NaN; hi=NaN;
n=size(tt,2);
u=zeros(1,n);
for k=1:n
  if(~isnan(lo) && ~isnan(hi) && ...
     t(lo) <= tt(k) && t(hi) > tt(k))
    u(k)=lo;
    % nothing to do, use cached values
  else
    % find interval using bisection
    lo=1; hi=size(t,2);
    while(hi-lo > 1)
      idx=floor((hi+lo)/2);
      if(t(idx) > tt(k))
        hi=idx;
      else
        lo=idx;
      end
    end
    % tt(k) is now in [t(lo),t(hi)),
    % store index value in u
    u(k)=lo;
  end
end

function test_locate()
%  1  2   3   4   5   6   7 
t=[0 0.2 0.4 0.5 0.7 0.9 1.0];
tt=[0.45 0.89 1.0 0.1 0.2 0.3 0.34 0.21 0.50001 0.82];
eu=[3 5 6 1 2 2 2 2 4 5]
u=locate(t,tt)
check_equal(u,eu,'u','eu');
fprintf('test_locate succeded\n')
