function [labelsCons, time] = HBGF(labelsEns, K)

% [labelsCons, time] = HBGF(labelsEns, K)
%--------------------------------------------------------------------------
% HBGF: Hybrid Bipartite Graph Formulation 
%
% Code available at: http://web.engr.oregonstate.edu/~xfern/HBGF_spec.m
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns		(matrix)	clustering ensemble; each column corresponds
%								to one clustering.
%	K				(scalar)	number of clusters in consensus clustering
%--------------------------------------------------------------------------
% OUTPUTS:
%   labelsCons		(vector)	ensemble clustering result - labels of
%								consensus clustering in a column vector
%	time			(scalar)	execution time in seconds
%
%------- LEGAL NOTICE -----------------------------------------------------
% This code is provided on "as is" basis for research use only.
% Copyright (C) 2004 Xiaoli Zhang Fern
% Part of Pepelka package.
%
%------- REFERENCE --------------------------------------------------------
% Fern, X. Z., & Brodley, C. E. (2004). Solving cluster ensemble problems
% by bipartite graph partitioning. 
% Proceedings of the 21st ICML (p. 36). New York, NY, USA: ACM.
%==========================================================================

% variables renaming
ceresults = labelsEns;
k = K;


% Use all the solutions in ensemble
sresults = ceresults;
[instnum,esize] = size(sresults);

tic;

% form the connectivity matrix instnum X totalk
w=formw(sresults, esize, instnum);

% compute the size of each cluster
csize=sum(w);

% normalize the vectors
w = w*diag(1./sqrt(csize));

% u is the eigenvector of normalized S(similarity matrix)
[u,s,v]= svds(w,k);
% normalize the rows
u = u./repmat(sqrt(sum(u.^2,2)), 1, k);
v = v./repmat(sqrt(sum(v.^2,2)), 1, k);

% find a good set of initial centers
idxperm=randperm(instnum);
centers(1,:) = u(idxperm(1),:);

c=zeros(instnum,1);
c(idxperm(1),1) = 2*k; %ensure this point is not to be selected again

for i=2:k

 c=c+abs(u*(centers(i-1,:))');
 [y,m] = min(c);
 centers(i,:) = u(m, :);
 c(m,1) = 2*k;
 
 end

% perform kmeans with the found initial centers
labelsCons = kmeans(vertcat(u,v), k, 'start', centers, 'maxiter', 200, 'emptyaction', 'singleton');
labelsCons(instnum+1:size(labelsCons,1), :) = [];

time=toc;
end

function w=formw(idxs, r, instnum)
% this function turns a cluster ensemble into a set of new binary features, one for each cluster
% if the object belongs to that cluster, the value is set to 1, and 0 otherwise
site=1;
for i=1:r
  idx = idxs(:,i);
  uniqueks=unique(idx(idx>0)); % 0 is considered to be the unclustered 
  t=zeros(instnum, length(uniqueks));
  t(repmat(idx, 1, length(uniqueks)) == repmat(uniqueks', instnum, 1))=1;
  w(:, site:site+length(uniqueks)-1) = t;
  site=site+length(uniqueks);
end
end

