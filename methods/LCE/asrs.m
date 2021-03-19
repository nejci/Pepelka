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
% optimization for speed: Nejc Ilc, 2014
%==========================================================================

[n,M] = size(E); %no of data points and no of clusterings
[E, no_allcl] = relabelCl(E); % re-labelling clusters in the ensemble E
wcl = weightCl(E,no_allcl);

%---find pair-wise similarity of clusters in each clustering using connected triple algorithm-----
nElem = sum(wcl>0,2);

CS = zeros(no_allcl,no_allcl); % create matrix CS, pair-wise similarity matrix for each pair of clusters
for i=1:no_allcl-1 % for each cluster
    Ni = wcl(i,:);
    for j=i+1:no_allcl % for other clusters
        Nj = wcl(j,:);
        prodNElem = nElem(i)*nElem(j);
        if prodNElem > 0
            CS(i,j) = sum(Ni.*Nj)/prodNElem;
        end
    end
end
maxCS = max(max(CS));
if maxCS > 0 
    CS = CS / maxCS;
end
CS=CS+CS';
CS(1:no_allcl+1:no_allcl^2)=1;

% CS = (wcl * wcl') ./ (nElem * nElem');
% maxCS2 = max(max(triu(CS,1)));
% if maxCS2 > 0 
%     CS = CS / maxCS2;
% end
% CS(1:no_allcl+1:no_allcl^2)=1;




%---construct pairwise similarity of each pair of data points -------------
S = zeros(n,n); % create matrix S, pair-wise similairty matrix for each pair of data points
for i = 1:n-1 % for each row
    for ii = i+1:n % for other rows (below i)
%         CSmat = CS(E(i,:),E(ii,:));
%         mask = CSmat == 1;
%         S(i,ii) = S(i,ii) + sum(sum(mask + (~mask.*(CSmat*dc))));
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
S(1:n+1:n^2) = 1;
