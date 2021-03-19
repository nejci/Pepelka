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
% optimization for speed: Nejc Ilc, 2014
%==========================================================================

[n,M] = size(E); %no of data points and no of clusterings
[E, no_allcl, cumCl] = relabelCl(E); % re-labelling clusters in the ensemble E

%==== initialize iteration r=0 =============================================
S = eye(n); %pairwise similarity between data points in the r-th iteration
C = eye(no_allcl); %pairwise similarity between clusters in the r-th iteration


%==== start iterations r1 to r(R-1) =====================================
for r = 1:R-1 %for each iteration
    
    %------ find similarity of each pair of data points --------------------------
     %pairwise similarity between data points in the 'r+1'-th iteration
    
    S1 = eye(n);
    for i = 1:n-1 %for each data point
        Ni = E(i,:); %clusters of data point i         
        for ii = i+1:n % for each other points
            Nii = E(ii,:);%clusters of data point ii 
            sum_sim = sum(sum(C(Ni,Nii)));
            S1(i,ii) = (dc/(M*M))*sum_sim;
        end
    end
       
    S1 = S1+S1';
    S1(1:n+1:n^2) = 1;   
    
    %-------find sim of each pair of clusters -----------------------------   
    
    C1 = eye(no_allcl); %pairwise similarity between clusters in the 'r+1'-th iteration
    Eind = ones(1,no_allcl);
    for m=2:M
        Eind(cumCl(m-1)+1:cumCl(m)) = m; 
    end
    for i = 1:no_allcl-1 %for each cluster        
        Ni = E(:,Eind(i))==i; %data belong to cluster i 
        for ii = i+1:no_allcl % for each other clusters
            Nii = E(:,Eind(ii))==ii; %data belong to cluster ii            
            sum_sim = sum(sum(S(Ni,Nii)));
            prodElem = sum(Ni)*sum(Nii);
            if prodElem > 0
                C1(i,ii) = (dc/prodElem)*sum_sim;
            end
        end
    end
    
    C1 = C1+C1';
    C1(1:no_allcl+1:no_allcl^2) = 1;
        
    %----------------------------------------------------------------------
    
    S = S1;
    C = C1;
end