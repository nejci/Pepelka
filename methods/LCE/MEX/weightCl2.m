function [wcl,pc]= weightCl2(E,no_allcl) 
%==========================================================================
% FUNCTION: wcl = weightCl(E)
% DESCRIPTION: This function computes weight for each pair of clusters using 
%              their shared members (Jaccard Coefficient)
%
% INPUTS:   E = N-by-M matrix of cluster ensemble
%
% OUTPUT: wcl = an weighted cluster matrix
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
% optimization for speed: Nejc Ilc, 2014
%==========================================================================

% N = size(E,1); %no of data points
% pc = zeros(N,no_allcl); % matrix indicates if data point belongs to the cluster (1=y, 0=n), row=data, col = cluster
% for i=1:N
%     pc(i,E(i,:))=1; % pc(i,j) = 1 if data i belongs to cluster j
% end
% 
% %find number of shared data points for each pair of clusters ==> intersect/union
% wcl = zeros(no_allcl,no_allcl);
% for i=1:no_allcl-1
%     for ii=i+1:no_allcl
%         pcSum = pc(:,i)+pc(:,ii);
%         tmp = sum(pcSum>0);
%         if tmp > 0
%             wcl(i,ii) = sum(pcSum==2) / tmp; %intersection/union
%         end
%     end
% end
% wcl = wcl + wcl';

N = size(E,1); %no of data points
pc = zeros(N,no_allcl); % matrix indicates if data point belongs to the cluster (1=y, 0=n), row=data, col = cluster
subInd = bsxfun(@plus,(E-1)*N,(1:N)');
pc(subInd(:)) = 1;
%find number of shared data points for each pair of clusters ==> intersect/union
wcl = zeros(no_allcl,no_allcl);
for i=1:no_allcl-1    
    ii=i+1:no_allcl;
    U = sum(bsxfun(@or,pc(:,i),pc(:,ii)),1); % union
    I = pc(:,i)' * pc(:,ii); % intersection
    tmp = I./U;
    tmp(~U) = 0;
    wcl(i,ii) = tmp;
end
wcl = wcl + wcl';