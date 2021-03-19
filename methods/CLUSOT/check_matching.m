function [b]=check_matching(d1,d2)
% [b]=check_matching(d1,d2)
%
% TODO: <short description>
%
% TODO: <explain vars>
% 
% $Id: check_matching.m 594 2005-04-18 16:52:24Z dome $
% D. Brugger, 20 March 2005
% test/check_matching.m

if(nargin == 0)
  test_check_matching();
  return;
end

if(ischar(d1) && ischar(d2))
  b=strcmp(d1,d2); 
elseif(isnumeric(d1) && isnumeric(d2))
  b=check_size(d1,d2);
  if(b)
    [n,m]=size(d1);
    % assume that NaN's at the same place are equal
    in1=isnan(d1); in2=isnan(d2);
    if(check_size(in1,in2) && sum(sum(in1 == in2)) == n*m)
      b=(sum(sum((d1 == d2) + isnan(d1))) == n*m);
    else
      b=(sum(sum(d1 == d2)) == n*m);
    end
  end
elseif(isstruct(d1) && isstruct(d2))
  %  fprintf('Both are structs!!!\n')
  b=check_size(d1,d2);
  if(b)
    fn1=fieldnames(d1); fn2=fieldnames(d2);
    b=check_matching(fn1,fn2);
    %  fprintf('Fieldnames equal!\n');
    if(b)
      b=check_size(fn1,fn2);
      if(b)
        [dim1,dim2]=size(d1);
        [n,m]=size(fn1);
        for kk=1:dim2
          d1_kk=d1(kk); d2_kk=d2(kk);
          for k=1:n
            for l=1:m
              fv1=getfield(d1_kk,fn1{k,l});
              fv2=getfield(d2_kk,fn2{k,l});
              b=check_matching(fv1,fv2);
              if(~b)
                return;
              end
            end
          end
        end
      end
    end
  end
elseif(iscell(d1) && iscell(d2))
  b=check_size(d1,d2);
  if(b)
    [n,m]=size(d1);
    for k=1:n
      for l=1:m
        b=check_matching(d1{k,l},d2{k,l});
        if(~b)
          return;
        end
      end
    end
  end
elseif(islogical(d1) && islogical(d2))
  b = d1 == d2;
else
  b=0;
end

function b=check_size(d1,d2)
[n1,m1]=size(d1); [n2,m2]=size(d2);
if(n1 ~= n2 || m1 ~= m2)
  b=0;
else
  b=1;
end

function test_check_matching()
% Test case #1 - two strings
s1='abc'; s2='def';
eb1=1; eb2=0;
b1=check_matching(s1,s1);
check_equal(b1,eb1,'b1','eb1');
b2=check_matching(s1,s2);
check_equal(b2,eb2,'b2','eb2');

% Test case #2 - two numbers
n1=3; n2=4;
eb1=1; eb2=0;
b1=check_matching(n1,n1);
check_equal(b1,eb1,'b1','eb1');
b2=check_matching(n1,n2);
check_equal(b2,eb2,'b2','eb2');

% Test case #3 - two matrices
m1=[1 2 0; 3 4 7]; m2=[1 1 0; 0 0 0];
eb1=1; eb2=0;
b1=check_matching(m1,m1);
check_equal(b1,eb1,'b1','eb1');
b2=check_matching(m1,m2);
check_equal(b2,eb2,'b2','eb2');

% Test case #4 - two structs
st1.a='abc'; st2.a='abc';
st1.b=3; st2.b=4;
eb1=1; eb2=0;
b1=check_matching(st1,st1);
check_equal(b1,eb1,'b1','eb1');
b2=check_matching(st1,st2);
check_equal(b2,eb2,'b2','eb2');

% Test case #5 - two cell arrays
ca1={1,'abc'}; ca2={'tds',2};
b1=check_matching(ca1,ca1);
check_equal(b1,eb1,'b1','eb1');
b2=check_matching(ca1,ca2);
check_equal(b2,eb2,'b2','eb2');

% Test case #6 - nested
nst1={st1, s1, m1};
nst2={st1, s1, m2};
b1=check_matching(nst1,nst1);
check_equal(b1,eb1,'b1','eb1');
b2=check_matching(nst1,nst2);
check_equal(b2,eb2,'b2','eb2');

% Test case #7 - two logicals
log1=true; log2=false;
b1=check_matching(log1,log1);
check_equal(b1,eb1,'b1','eb1');
b2=check_matching(log1,log2);
check_equal(b2,eb2,'b2','eb2');

fprintf('test_check_matching succeded!\n');
