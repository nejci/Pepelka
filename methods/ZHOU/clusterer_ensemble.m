function micro_p = clusterer_ensemble(dataset,k,varargin)
%
% create an ensemble of clusterers
% 
% ATTN: This package is free for academic usage. The code was developed by Mr. W. Tang (wtang314@yahoo.com). You can run
% it at your own risk. For other purposes, please contact Prof. Zhi-Hua Zhou (zhouzh@nju.edu.cn)
%
% micro_p   - the micro-p metric for evaluating the quality of the clustering
%
% dataset   - file name for data set, e.g. 'c:\clusterer_ensemble\clusterdata.mat'. This file should contain two parts, 
%             i.e. inputs and labels. The data set is represented by a feature matrix data(n,m), where n is the number
%             of attributes while m is the number of instances. label is a class vector represents the prior knowledge 
%             about the distribution of the data.
%             for example, the data and label for parity problem are:
%                  data:                  label:      
%                  0  0  1  1        
%                  0  1  0  1              1 2 2 1
%                  1  1  1  1
%             Note that here the labels are required because micro-p is employed to evaluate the quality of the clustering
%
% k         - the number of data groups to be clustered
%
% varargin  - an input varible with three arguments: 
%       'populationsize'    - the ensemble size, the default value is 20
%       'selective'         - 1 means selective ensemble, 0 means non-selective ensemble, the default value is 0
%       'weighted'          - 1 means weighted voting, 0 means unweighted voting, the default value is 1
%
%       The default values can be changed as shown by the following examples:
%           sel_w_voting('c:\clusterer_ensemble\clusterdata.mat',4,'populationsize',30)     - set the ensemble size as 30
%           sel_w_voting('c:\clusterer_ensemble\clusterdata.mat',3,'selective',1)           - to build selective ensemble
%           sel_w_voting('c:\clusterer_ensemble\clusterdata.mat',3,'populationsize',50,'selective',1,'bootstrapsample',1)
%                                                                                           - set three arguments together
%
%
% Reference: Z.-H. Zhou and W. Tang. Clusterer ensemble. Knowledge-Based Systems, 2006, 19(1): 77-83. 
%
% ATTN2: This package was developed by Mr. W. Tang (wtang314@yahoo.com). For any problem concerning the code, please feel
% free to contact Mr. Tang.
%

tic;
%rand('state',sum(100 * clock));
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
 
if nargin < 2
    error('at least two arguments required!');
end

argnames = {'populationsize' 'selective' 'weighted'};
defaults = {20 0 1};
[errmsg,populationsize,selective,weighted] = getargs(argnames,defaults,varargin{:});
error(errmsg);

if selective ~= 0 && selective ~= 1 && weighted ~= 0 && weighted ~= 1
    error('wrong input arguments!');
end

if exist(dataset,'var')
    [filebase,filename] = fileparts(dataset);
    load(dataset);
    new_dir = strcat(filebase,'\',filename);
    if ~exist(new_dir,'dir')    % if the new directory does not exist, it is created
        mkdir(filebase,filename);
    end
    filebase = new_dir;
else
    error('file does not exist!');
end

[datadim,datanum] = size(data);

% normalize the original data into 0-1 scale
minimum = min(data,[],2);
maximum = max(data,[],2);
for i = 1:datanum
    data(:,i) = (data(:,i) - minimum) ./ (maximum - minimum);
end

filename = strcat(filebase,'\','cluster_components.mat');
if exist(filename,'file')
    load(filename);
else
    % generate the component clusterers
    cluster_components = zeros(populationsize,datanum);
    for i = 1:populationsize     % here the base learner is k-means implemented in km.m
        lamda = km(data,k);      % run k-means without bootstrap sample
        cluster_components(i,:) = lamda;
    end
        
    % train the clusterer ensemble
    % align all the component clusterers
    baseline = cluster_components(floor(rand * size(cluster_components,1)) + 1,:);   % randomly select a component as the baseline
    for i = 1:populationsize
        matrix_index = label2matrix(cluster_components(i,:),k);
        matrix_overlap = count(baseline,cluster_components(i,:));
        position = zeros(k);
        for j = 1:k
            tmpidx = find(max(matrix_overlap(:)) == matrix_overlap);                 % find the number of elements with the maximum value
            % if there exist more than one maximum, select the one that contains the minimum number of instances in the i-th cluster
            if length(tmpidx) > 1
                minval = datanum;
                for n = 1:length(tmpidx)
                    if length(find(cluster_components(i,:) == tmpidx(n))) < minval
                        maxidx = tmpidx(n);
                        minval = length(find(cluster_components(i,:) == tmpidx(n)));
                    end
                end
            else
                maxidx = tmpidx;
            end
            position(mod(maxidx - 1,k) + 1,ceil(maxidx / k)) = 1;
            matrix_overlap(mod(maxidx - 1,k) + 1,:) = 0;
            matrix_overlap(:,ceil(maxidx / k)) = 0;
        end
        [maxval,idx] = max(position);                                                % maxval is useless
        [maxval,cluster_components(i,:)] = max(matrix_index(:,idx)');
    end
    
    % compute the weight according to the mutual information
    weight = ones(1,populationsize);
    if populationsize ~= 1
        matrix_NMI = zeros(populationsize);
        for i = 1:populationsize
            for j = 1:i - 1
                matrix_NMI(i,j) = compute_NMI(cluster_components(i,:),cluster_components(j,:),k);
                matrix_NMI(j,i)  = matrix_NMI(i,j);
            end
        end
        reciprocal_belta = (populationsize - 1) ./ sum(matrix_NMI);
        weight = reciprocal_belta ./ sum(reciprocal_belta);
    end
    save(filename,'cluster_components','weight');
end


% test the clusterer ensemble

% combine the component clusterers based on the weight
if selective
    selected = ((weight > mean(weight)) | (abs(weight - mean(weight)) < 1e-10));
    if weighted                              % selective weighted voting 
        consensus_cluster = weighted_vote(cluster_components(selected,:),weight(selected),k);
        filename = strcat(filebase,'\','s_w_consensus.mat');
    else
        weight = ones(1,populationsize);     % selective voting
        consensus_cluster = weighted_vote(cluster_components(selected,:),weight(selected),k);        
        filename = strcat(filebase,'\','s_nw_consensus.mat');
    end
    save(filename,'consensus_cluster','selected');
else
    if weighted                              % weighted-voting   
        consensus_cluster = weighted_vote(cluster_components,weight,k);
        filename = strcat(filebase,'\','ns_w_consensus.mat');
    else                                     % voting                              
        weight = ones(1,populationsize);
        consensus_cluster = weighted_vote(cluster_components,weight,k);
        filename = strcat(filebase,'\','ns_nw_consensus.mat');
    end
    save(filename,'consensus_cluster');
end

% evaluate the quality of the clustering with micro-p
micro_p = micro_precision(consensus_cluster,label);
toc;


%--------------------------------------------------------------------------
%  This function generates cluster label based on the centres of k-means
%--------------------------------------------------------------------------
function index = generatelabel(centres,data)
distance = zeros(size(centres,2),size(data,2));    % the matrix of distance between cluster centres and the data
for i = 1:size(centres,2)
    for j = 1:size(data,2)
        distance(i,j) = sum(abs(centres(:,i) - data(:,j)).^2).^(1 / 2);
    end
end
[minval,index] = min(distance);                    % the minval is useless


%--------------------------------------------------------------------------
%  This function changes a label vector to a matrix presentation
%--------------------------------------------------------------------------   

%--------------------------------------------------------------------------   
%   For example,
%      if       input   = 1 2 2 3 1 3 2     k = 3
%
%      then     output  = 1 0 0 0 1 0 0
%                         0 1 1 0 0 0 1
%                         0 0 0 1 0 1 0 
%--------------------------------------------------------------------------
function output = label2matrix(input,k);
id = eye(k);
output = id(input,:);


%--------------------------------------------------------------------------
%  This function caculates the overlap matrix between two labels 
%--------------------------------------------------------------------------
function matrix_overlap = count(lamda_x,lamda_y)
k = length(unique(lamda_x)); % get the number of clusters
matrix_overlap = zeros(k);
for i = 1:length(lamda_x)
    for j = 1:length(lamda_y)
        matrix_overlap(lamda_x(i),lamda_y(j)) = matrix_overlap(lamda_x(i),lamda_y(j)) + 1;
    end
end


%--------------------------------------------------------------------------
%  This function computes the normalized mutual information 
%--------------------------------------------------------------------------
function NMI = compute_NMI(lamda_x,lamda_y,k)
NMI = 0;
for i = 1:k
    for j = 1:k
        if sum((lamda_x == i) & (lamda_y == j))   % ensure the number of instances in either cluster is zero
            NMI = NMI + sum((lamda_x == i)&(lamda_y == j)) * log(sum((lamda_x == i)...
                &(lamda_y == j)) * length(lamda_x) / (sum(lamda_x == i) * sum(lamda_y == j))) / log(k * k);
        end
    end
end
NMI = NMI * 2 / length(lamda_x);                  % the number of instances in lamda_x is equal to that in lamda_y


%--------------------------------------------------------------------------
%   This function votes the component clusterers according to given weight
%   vector. Note that when all the elements of the weight vector is 1, 
%   this function votes the component clusteres without weight
%--------------------------------------------------------------------------
function consensus = weighted_vote(components,weight,k)
id = eye(k);
[componentno,exampleno] = size(components);
consensus = zeros(1,exampleno);
for i = 1:exampleno
    tmp = label2matrix(components(:,i),k);
    for j = 1:componentno
        tmp(j,:) = tmp(j,:) * weight(j);
    end
    if length(find(sum(tmp,1) == max(sum(tmp,1)))) == 1
        [maxoutput,consensus(1,i)] = max(sum(tmp,1));       % maxoutput is useless
    end
end

% if tie appears
tag = unique(consensus(find(consensus ~= 0)));
label = 1:k;
for i = 1:length(tag)
    label = label(1,label ~= tag(i));
end
if ~isempty(find(consensus == 0))
    if isempty(label)   % if all the winners have got some instances, then current instance is randomly assigned to one of the winners
        for i = find(consensus == 0)
            consensus(1,i) = floor(rand * k) + 1;
        end
    else                % if some winners have not got instances, then current instance is randomly assigned to one of them
        for i = find(consensus == 0)
            index = floor(rand * length(label)) + 1;
            consensus(1,i) = label(index);
        end
    end
end


%--------------------------------------------------------------------------
%  This function computes micor-p
%--------------------------------------------------------------------------
function mp = micro_precision(T,D)
class = unique(D);
mp = 0;
for i = 1:length(class)
    tmp = label2matrix(T(find(D == class(i))),length(class));
    mp = mp + max(sum(tmp));
end
mp = mp / length(D);

