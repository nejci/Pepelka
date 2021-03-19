function [index,centres] = km(data,k,varargin)
%
% implements the standard k-means clustering algorithm for clusterer_ensemble.m
%
% ATTN: This package is free for academic usage. The code was developed by Mr. W. Tang (wtang314@yahoo.com). You can run
% it at your own risk. For other purposes, please contact Prof. Zhi-Hua Zhou (zhouzh@nju.edu.cn)
%
% index     - the index, or order, of the clusters
%
% centres   - the centres of the clusters
%
% data      - a feature matrix data(n,m), where m is the number of instances and n is the number of attributes 
%
% k         - the number of data groups to be clustered
% 
% varargin  - an input varible with five arguments:
%       'metric'          - the similarity metric, the default value is Minkowski distance
%       'maxiteration'    - the maximum iteration to be executed, the default value is 100
%       'errorthreshold'  - the error threshold for terminating the iterative process, the default value is 1e-5
%       'power'           - the power of the Minkowski distance, the default value is 2, 
%                           that is, the default metric is Euclidean distance
%       'centres'         - the initial centres, the default value is randomly selected
%
%       The default values can be changed as shown by the following examples:
%           km(data,3,'metric','Func')      - set the metric as the function Func you specified
%           km(data,4,'maxiteration',200)   - set the maximum iteration to be executed as 200
%           km(data,5,'power',1)            - use Manhattan distance
%           km(data,2,'maxiteration',200,'errorthreshold',1e-4)
%                                           - set both the maximum iteration and the error threshold 
%
% ATTN2: This package was developed by Mr. W. Tang (wtang314@yahoo.com). For any problem concerning the code, please feel
% free to contact Mr. Tang.
%

if nargin < 2
    error('at least two arguments required.');
end

[data_dim,data_num] = size(data);
perm = randsample(data_num,k);      %randomly select k instances as the initial centres, refer the function randsample below
pnames = {'metric' 'maxiteration' 'errorthreshold','power','centres'};
defaults = {'minkowski' 100 1e-5 2 data(:,perm)};
[errmsg,metric,maxiteration,errorthreshold,power,centres] = getargs(pnames, defaults, varargin{:});
error(errmsg);

if k > data_num
    error('more cluster centres than data!');
end

id = eye(k);

disp('start k-means clustering:');
for i = 1:maxiteration
    old_centres = centres;
    distance = feval(metric, centres, data, power);
    [minval,index] = min(distance);
    position = id(index,:);
    points_num = sum(position);
    for j = 1:k
        if points_num(j) > 0
            centres(:,j) = sum(data(:,find(position(:,j))),2) / points_num(j);
        end
    end
    err = sum(minval);
    fprintf(1,'iteration:%4d  error:%11.6f\n',i,err);
    if i > 1
        if (max(max(abs(centres - old_centres))) < errorthreshold) & (abs(old_err - err) < errorthreshold)
            disp('error threshold reached!');
            return;
        end
    end
    old_err = err;
end
disp('maximum iteration reached!');
return;


%--------------------------------------------------------------------------
% This function implements the Minkowski distance
%--------------------------------------------------------------------------
function dis = minkowski(x,y,power)
dis = zeros(size(x,2), size(y,2));
for i = 1:size(x,2)
    for j = 1:size(y,2)
        dis(i,j) = sum(abs(x(:,i) - y(:,j)) .^ power).^(1 / power);
    end
end


%--------------------------------------------------------------------------
% This function randomly selects k instances to be the initial centres
%--------------------------------------------------------------------------
function y = randsample(n, k)
if 4 * k > n
    rp = randperm(n);
    y = rp(1:k);
else
    x = zeros(1,n);
    sumx = 0;
    while sumx < k
        x(ceil(n * rand(1,k - sumx))) = 1;
        sumx = sum(x);
    end
    y = find(x > 0);
    y = y(randperm(k));
end
