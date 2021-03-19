function [d]=calc_d(inp,iv)
% [d]=calc_d(inp,iv)
%
% calc_d determines the euclidian distance between neighboring
% neurons of a som or between vectors in data space.
%
% Input:
% inp - (som_struct) trained som_map or training data som_data. 
%
% Output:
% d - (matrix) with distances (num_n x num_n)
% iv - (vector) interval used to normalize distance values, e.g.
%        [0, 0.99]. If not given values are left unchanged. [optional]
%
% TODO: Need more test cases
%
% $Id: calc_d.m 853 2005-07-25 20:44:00Z dome $
% D. Brugger, 01 February 2005
% som/calc_d.m

if(nargin == 0)
  test_calc_d();
  return;
end

if(nargin < 2)
  flag=0;
else
  flag=1;
end

if(isstruct(inp))
  switch lower(inp.type)
   case 'som_map'
    % number of neurons
    num_n = size(inp.codebook,1);
    % weight matrix of size (num_n,dim)
    w = inp.codebook;
    C = som_unit_neighs(inp);
    d = ones(num_n,num_n).*Inf;
    if(flag) % normalize distances to iv
      % check if distances have already been computed,
      % this is the case if som_kbatch was used for training
      if(isfield(inp,'cbv_d'))
        d = inp.cbv_d;
        for k=1:num_n
          % get index of neighbors at distance 1
          n_ind=find(C(k,:) == 1);
          for l=1:length(n_ind)
            if(k==1 && l==1)
              minimum=d(k,n_ind(l));
              maximum=d(k,n_ind(l));
            else
              if(d(k,n_ind(l)) > maximum)
                maximum=d(k,n_ind(l));
              end
              if(d(k,n_ind(l)) < minimum)
                minimum=d(k,n_ind(l));
              end
            end
          end
        end
      else
        for k=1:num_n
          % get index of neighbors at distance 1
          n_ind=find(C(k,:) == 1);
          for l=1:length(n_ind)
            d(k,n_ind(l))=sqrt(sum((w(k,:)-w(n_ind(l),:)).^2));
            if(k==1 && l==1)
              minimum=d(k,n_ind(l));
              maximum=d(k,n_ind(l));
            else
              if(d(k,n_ind(l)) > maximum)
                maximum=d(k,n_ind(l));
              end
              if(d(k,n_ind(l)) < minimum)
                minimum=d(k,n_ind(l));
              end
            end
          end
        end
      end
      d = lin_scale(d,[minimum maximum],iv);
    else
      for k=1:num_n
        % get index of neighbors at distance 1
        n_ind=find(C(k,:) == 1);
        for l=1:length(n_ind)
          d(k,n_ind(l))=sqrt(sum((w(k,:)-w(n_ind(l),:)).^2));
        end
      end
    end
    
   case 'som_data'
    num_n=size(inp.data,1);
    w = inp.data;
    d = zeros(num_n,num_n);
    if(flag) % normalize distances to iv
      for k=1:num_n
        for l=k:num_n
          d(k,l)=sum((w(k,:)-w(l,:)).^2);
          if(k==1 && l==1)
            minimum=d(k,l);
            maximum=d(k,l);
          else
            if(d(k,l) > maximum)
              maximum=d(k,l);
            end
            if(d(k,l) < minimum)
              minimum=d(k,l);
            end
          end
        end
      end
      d=sqrt(d);
      d=triu(d,1)+d';
      d = lin_scale(d,[minimum maximum],iv);
    else
      for k=1:num_n
        for l=k:num_n
          d(k,l)=sum((w(k,:)-w(l,:)).^2);
        end
      end
      d=sqrt(d);
      d=triu(d,1)+d';
    end
    
   otherwise
    error('!!! Unkown input type ''%s'' !!!', inp.type)
  end
else
  num_n = size(inp,1);
  w = inp;
  d = zeros(num_n,num_n);
  if(flag) % normalize distances to iv
    for k=1:num_n
      for l=k:num_n
        d(k,l)=sum((w(k,:)-w(l,:)).^2);
        if(k==1 && l==1)
          minimum=d(k,l);
          maximum=d(k,l);
        else
          if(d(k,l) > maximum)
            maximum=d(k,l);
          end
          if(d(k,l) < minimum)
            minimum=d(k,l);
          end
        end
      end
    end
    d=sqrt(d);
    d=triu(d,1)+d';
    d = lin_scale(d,[minimum maximum],iv);
  else
    for k=1:num_n
      for l=k:num_n
        d(k,l)=sum((w(k,:)-w(l,:)).^2);
      end
    end
    d=sqrt(d);
    d=triu(d,1)+d';
  end
end

% check if result is valid
%mn = max(sum(d > 0,2));
%for k=1:num_n
%  for l=1:num_n
%    if(d(k,l) < 0)
%      error('!!! Distances should be positive !!!');
%    end
%    if(d(k,l) ~= d(l,k))
%      error('!!! Matrix d should be symmetric !!!');
%    end
%  end
%end

function test_calc_d()
test.data=[1 0 0; ...
           0 1 0; ...
           0 0 0];
test.type='som_data';
ed = [0 sqrt(2) 1; ...
      sqrt(2) 0 1; ...
      1 1 0];
d=calc_d(test);

if(sum(sum(ed==d)) ~= 9)
  error('!!! Test failed !!!');
end
fprintf('test_calc_d succeded\n');
