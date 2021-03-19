function [wcl,pc] = weightCl(E) 
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
%==========================================================================

N = size(E,1); %no of data points
no_allcl = max(max(E));
pc = zeros(N,no_allcl); % matrix indicates if data point belongs to the cluster (1=y, 0=n), row=data, col = cluster
for i=1:N
    pc(i,E(i,:))=1; % pc(i,j) = 1 if data i belongs to cluster j
end

%find number of shared data points for each pair of clusters ==> intersect/union
wcl = zeros(no_allcl,no_allcl);
for i=1:no_allcl-1
    for ii=i+1:no_allcl
        tmp = sum((pc(:,i)+pc(:,ii))>0);
        if tmp > 0
            wcl(i,ii) = sum((pc(:,i)+pc(:,ii))==2) / tmp; %intersection/union
        end
    end
end
wcl = wcl + wcl';