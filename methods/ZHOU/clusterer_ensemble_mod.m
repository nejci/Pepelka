function [labelsCons, K, time] = clusterer_ensemble_mod(labelsEns,k,varargin)
%
% [labelsCons, time] = clusterer_ensemble_mod(labelsEns,k,varargin)
%--------------------------------------------------------------------------
% Create a consensus from ensemble of clusterings implementing:
%	- voting
%	- weighted voting
%	- selective voting
%	- selective weighted voting
%
%	!!!  Every clustering in the ensemble has to have same number of
%	clusters k !!!
%
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns	(matrix)	clustering ensemble; each COLUMN corresponds
%							to one clustering (standard in Pepelka).
%
%	k			(scalar)	the number of data groups to be clustered
%
%	varargin	(list)		an input variable with two property-value pair: 
%							'selective'			1 : selective ensemble 
%												0 : non-selective ensemble [default]
%							'weighted'          1 : weighted voting [default] 
%												0 : unweighted voting
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   labelsCons		(vector)	ensemble clustering result - labels of
%								consensus clustering in a column vector
%	time			(scalar)	execution time in seconds
%
%------- LEGAL NOTICE -----------------------------------------------------
% This package is free for academic usage. 
% The code was developed by Mr. W. Tang (wtang314@yahoo.com). 
% You can run it at your own risk.
% Part of Pepelka package.
%
%------- REFERENCE --------------------------------------------------------
% Z.-H. Zhou and W. Tang. 
% Clusterer ensemble. Knowledge-Based Systems, 2006, 19(1): 77-83.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 11-May-2012 by Nejc Ilc
%
% Modification by Nejc to the original code written by W. Tang:
%	- accepts ensemble labels as input
%	- does not generate ensemble
%	- does not calculate micro-precision validity measure
%	- does not read/save files
%
%------- CONTACT ----------------------------------------------------------
% Please write to:	W. Tang (wtang314@yahoo.com) or 
%					Zhi-Hua Zhou (zhouzh@nju.edu.cn)
%==========================================================================

% Test
% [labelsEns, moreInfo]=pplk_genEns(data, {'KM',20,2,'fixed'}, []);
% [labelsCons, K,time] = clusterer_ensemble_mod(labelsEns,2,'selective',1,'weighted','1');
% CVI = pplk_validExt(target,labelsCons,{'CA','NMI','VM','B3E'});

if any(max(labelsEns) ~= k)
	error(['Number of clusters in ensemble does not match k (',num2str(k),').']);
end

% initialize random generator
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

tic;

if nargin < 2
    error('at least two arguments required!');
end

argnames = {'selective' 'weighted'};
defaults = {0 1};
[errmsg,selective,weighted] = getargs(argnames,defaults,varargin{:});
error(errmsg);

if selective ~= 0 && selective ~= 1 && weighted ~= 0 && weighted ~= 1
    error('wrong input arguments!');
end


cluster_components = labelsEns';

[populationsize,datanum] = size(cluster_components);

% train the clusterer ensemble
% align all the component clusterers

% 1. randomly select a component as the baseline
baseline = cluster_components(floor(rand * size(cluster_components,1)) + 1,:);
for i = 1:populationsize
	matrix_index = label2matrix(cluster_components(i,:),k);
	matrix_overlap = count(baseline,cluster_components(i,:));
	position = zeros(k);
	for j = 1:k
		% find the number of elements with the maximum value
		tmpidx = find(max(matrix_overlap(:)) == matrix_overlap);
		
		% if there exist more than one maximum, select the one that
		% contains the minimum number of instances in the i-th cluster
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
	[maxval,idx] = max(position);
	[maxval,cluster_components(i,:)] = max(matrix_index(:,idx)');
end

% compute the weight according to the mutual information
weight = ones(1,populationsize);
if populationsize ~= 1
	matrix_NMI = zeros(populationsize);
	for i = 1:populationsize
		for j = 1:i - 1
			matrix_NMI(i,j) = compute_NMI(cluster_components(i,:),cluster_components(j,:),k);
			matrix_NMI(j,i) = matrix_NMI(i,j);
		end
	end
	reciprocal_belta = (populationsize - 1) ./ sum(matrix_NMI);
	weight = reciprocal_belta ./ sum(reciprocal_belta);
end

% Consensus function
% combine the component clusterers based on the weight
if selective
    selected = ((weight > mean(weight)) | (abs(weight - mean(weight)) < 1e-10));
    
	% selective weighted voting 
	if weighted
        consensus_cluster = weighted_vote(cluster_components(selected,:),weight(selected),k);
    
	% selective voting 
	else
        weight = ones(1,populationsize);
        consensus_cluster = weighted_vote(cluster_components(selected,:),weight(selected),k);        
	end
    
else
	% weighted-voting 
	if weighted
        consensus_cluster = weighted_vote(cluster_components,weight,k);

	% voting
	else
        weight = ones(1,populationsize);
        consensus_cluster = weighted_vote(cluster_components,weight,k);
	end
end

time=toc;
labelsCons = consensus_cluster';
K = length(unique(labelsCons));

% if K ~= k
% 	warning('Pepelka:numberOfClusters',['Number of obtained clusters (',num2str(K),') different from specified (',num2str(k),')!']);
% end

end % function

%--------------------------------------------------------------------------
%  This function changes a label vector to a matrix presentation
%--------------------------------------------------------------------------  
%   For example,
%      if       input   = 1 2 2 3 1 3 2     k = 3
%
%      then     output  = 1 0 0 0 1 0 0
%                         0 1 1 0 0 0 1
%                         0 0 0 1 0 1 0 
%--------------------------------------------------------------------------
function output = label2matrix(input,k)
id = eye(k);
output = id(input,:);
end % function

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
end % function

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
end % function

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
        [maxoutput,consensus(1,i)] = max(sum(tmp,1));
    end
end

% if tie appears
tag = unique(consensus(consensus ~= 0));
label = 1:k;
for i = 1:length(tag)
    label = label(1,label ~= tag(i));
end
if ~isempty(find(consensus == 0, 1))
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
end % function