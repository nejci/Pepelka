function S = cts0(E, dc) 
%==========================================================================
% FUNCTION: S = cts(E, dc) 
% DESCRIPTION: A funciton for computing Connected-Triple Based Similarity matrix
%
% INPUT:  E = matrix of cluster ensemble
%        dc = decay factor, ranges [0,1]
%
% OUTPUT: S = CTS matrix
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

[n,M] = size(E); %no. of data points and no. of clusterings
[E, no_allcl] = relabelCl(E); % re-labelling clusters in the ensemble E
wcl = weightCl0(E);

%---find pair-wise similarity of clusters in each clustering using connected triple algorithm-----
wCT = zeros(no_allcl,no_allcl); % create matrix wCT (weighted-connected trple of clusters), pair-wise similarity matrix for each pair of clusters
maxCl = max(E);
minCl = min(E);
for q = 1:M % for each clustering
    for i=minCl(q):maxCl(q)-1 %for each cluster
        Ni = wcl(i,:);
        for j=i+1:maxCl(q) %for other clusters
            Nj = wcl(j,:);
            wCT(i,j) = sum(min(Ni,Nj));
        end
    end
end
if max(max(wCT)) > 0 
    wCT = wCT / max(max(wCT));
end
wCT = wCT + wCT';
for i = 1:no_allcl
    wCT(i,i) = 1;
end

%---find pair-wise similarity of data points--------------------------------
S = zeros(n,n); % create matrix S, pair-wise similairty matrix for each pair of data points
for m = 1:M
    for i = 1:n-1 % for each row, start at row #2 
        for ii = i+1:n % for other rows (below i)
            if E(i,m) == E(ii,m)
                S(i,ii) = S(i,ii)+1;
            else
                S(i,ii) = S(i,ii) + (wCT(E(i,m),E(ii,m))*dc);
            end              
        end
    end
end
S = S/M;
S = S + S';
for i = 1:n
    S(i,i) = 1;
end
