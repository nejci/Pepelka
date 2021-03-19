function [f,F]=calc_f(map, td, iv)
% [f,F]=calc_f(map, td, iv)
%
% calc_f determines the frequencies of neurons in map.
% By definition the total number of training patterns
% mapped to a specific neuron.
%
% Inputs:
% map - trained som
% td  - training data
%
% Output:
% f - (vector) neuron frequencies, row vector (1,num_n).
% F - (matrix) where F(k,l)=|f(k)-f(l)|
% iv - (vector) interval used to normalize frequency values in F, e.g.
%        [0, 0.99]. If not given values are left unchanged. [optional]
% 
% $Id: calc_f.m 520 2005-03-30 11:54:44Z dome $
% D. Brugger, 01 February 2005
% som/calc_f.m


if(nargin < 3)
  flag=0;
else
  flag=1;
end

% determine best matching units (BMUs)
bmus=som_bmus(map,td);
% number of neurons
[num_n, dump] = size(map.codebook);
% number of training patterns num_p
[num_p, dump] = size(td.data);
f=zeros(1,num_n);
for k=1:num_n
  f(k)=sum(bmus == k);
end

C = som_unit_neighs(map);
F = ones(num_n,num_n).*Inf;
if(flag)
  for k=1:num_n
    % get index of neighbors at distance 1
    n_ind=find(C(k,:) == 1);
    for l=1:length(n_ind)
      F(k,n_ind(l))=abs(f(k)-f(n_ind(l)));
      if(k==1 && l==1)
        minimum=F(k,n_ind(l));
        maximum=F(k,n_ind(l));
      else
        if(F(k,n_ind(l)) > maximum)
          maximum=F(k,n_ind(l));
        end
        if(F(k,n_ind(l)) < minimum)
          minimum=F(k,n_ind(l));
        end
      end
    end
  end
  F = lin_scale(F,[minimum maximum],iv);
else
  for k=1:num_n
    % get index of neighbors at distance 1
    n_ind=find(C(k,:) == 1);
    for l=1:length(n_ind)
      F(k,n_ind(l))=abs(f(k)-f(n_ind(l)));
    end
  end
end

% check validity
if(sum(f) ~= num_p)
  error(['!!! Sum of neuron frequencies should match total number' ...
	 ' of training patterns  !!!'])
end
