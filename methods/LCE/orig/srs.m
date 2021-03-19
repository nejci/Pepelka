function S = srs(E, dc, R)
%==========================================================================
% FUNCTION: S = srs(E, dc, R)
% DESCRIPTION: A funciton for computing SimRank Based Similarity matrix
%
% INPUT:  E = matrix of cluster ensemble
%        dc = decay factor, ranges [0,1]
%         R = the number of iterations for SimRank algorithm
%
% OUTPUT: S = SRS matrix
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

[n,M] = size(E); %no of data points and no of clusterings
[E, no_allcl] = relabelCl(E); % re-labelling clusters in the ensemble E

%==== initialize iteration r=0 =============================================
S = diag(ones(n,1)); %pairwise similarity between data points in the r-th iteration
C = diag(ones(no_allcl,1)); %pairwise similarity between clusters in the r-th iteration

%==== start iterations r1 to r(R-1) =====================================
for r = 1:R-1 %for each iteration
    
    %------ find similarity of each pair of data points --------------------------
    S1 = diag(ones(n,1)); %pairwise similarity between data points in the 'r+1'-th iteration
    for i = 1:n-1 %for each data point
        Ni = E(i,:); %clusters of data point i         
        for ii = i+1:n % for each other points
            sum_sim = 0;
            Nii = E(ii,:);%clusters of data point ii 
            for k = 1:M
                for kk =1:M
                    sum_sim = sum_sim + C(Ni(k),Nii(kk));
                end
            end
            S1(i,ii) = (dc/(M*M))*sum_sim;
        end
    end
    S1 = S1+S1';
    for i=1:n
        S1(i,i)=1;
    end    
    
    %-------find sim of each pair of clusters -----------------------------
    C1 = diag(ones(no_allcl,1)); %pairwise similarity between clusters in the 'r+1'-th iteration
    for i = 1:no_allcl-1 %for each cluster
        [Ni,col]=find(E==i); %data belong to cluster i 
        nki = length(Ni);
        for ii = i+1:no_allcl % for each other clusters
            sum_sim = 0;
            [Nii,col]=find(E==ii); %data belong to cluster ii            
            nkii = length(Nii);            
            for k = 1:nki
                for kk =1:nkii
                    sum_sim = sum_sim + S(Ni(k),Nii(kk));
                end
            end
            if (nki*nkii) > 0
                C1(i,ii) = (dc/(nki*nkii))*sum_sim;
            end
        end
    end
    C1 = C1+C1';
    for i=1:no_allcl
        C1(i,i)=1;
    end        
    %----------------------------------------------------------------------
    
    S = S1;
    C = C1;
end