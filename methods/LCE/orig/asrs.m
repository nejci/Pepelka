function S = asrs(E, dc)
%==========================================================================
% FUNCTION: S = asrs(E, dc) 
% DESCRIPTION: A funciton for computing Approximated SimRank Based Similarity matrix
%
% INPUT:  E = matrix of cluster ensemble
%        dc = decay factor, ranges [0,1]
%
% OUTPUT: S = ASRS matrix
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

[n,M] = size(E); %no of data points and no of clusterings
[E, no_allcl] = relabelCl(E); % re-labelling clusters in the ensemble E
wcl = weightCl(E);

%---find pair-wise similarity of clusters in each clustering using connected triple algorithm-----
CS = zeros(no_allcl,no_allcl); % create matrix CS, pair-wise similarity matrix for each pair of clusters
for i=1:no_allcl-1 % for each cluster
    Ni = wcl(i,:);
    ni = length(Ni(Ni>0)); % no. of neightbors of i
    for j=i+1:no_allcl % for other clusters
        Nj = wcl(j,:);
        nj = length(Nj(Nj>0)); % no. of neightbors of j
        if (ni*nj) > 0
            CS(i,j) = sum(Ni.*Nj)/(ni*nj);
        end
    end
end
if max(max(CS)) > 0 
    CS = CS / max(max(CS));
end
CS=CS+CS';
for i=1:no_allcl
    CS(i,i)=1;
end

%---construct pairwise similarity of each pair of data points -------------
S = zeros(n,n); % create matrix S, pair-wise similairty matrix for each pair of data points
for i = 1:n-1 % for each row
    for ii = i+1:n % for other rows (below i)
        for j = 1:M
            for jj=1:M
                if CS(E(i,j),E(ii,jj)) == 1
                    S(i,ii) = S(i,ii) + 1;
                else
                    S(i,ii) = S(i,ii) + (CS(E(i,j),E(ii,jj))*dc);
                end
            end
        end
    end
end
S = S / (M*M);
S = S + S';
for i = 1:n
    S(i,i) = 1;
end
