function S = cts2(E, dc) 
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
% optimization for speed: Nejc Ilc, 2014
%==========================================================================

[n,M] = size(E); %no. of data points and no. of clusterings
[E, nCls, C] = relabelCl(E); % re-labelling clusters in the ensemble E
% [wcl, pc] = weightCl2(E,nCls);
[wcl,pc]= weightCl(E,nCls);


%---find pair-wise similarity of clusters in each clustering using connected triple algorithm-----
wCT = zeros(nCls,nCls); % create matrix wCT (weighted-connected trple of clusters), pair-wise similarity matrix for each pair of clusters
q = 1;
for i=1:nCls-1 %for each cluster
    Ni = wcl(i,:);
    % determine index of ensemble member
    if i > C(q)
        q = q+1;
    end
    j=i+1:C(q); %for other clusters
    Nj = wcl(j,:);
    wCT(i,j) = sum(bsxfun(@min,Ni,Nj),2);    
end

% wCT = zeros(nCls,nCls);
% minCl = [1 C(1:end-1)+1];
% for q = 1:M % for each clustering
%     for i=minCl(q):C(q)-1 %for each cluster
%         Ni = wcl(i,:);
%         for j=i+1:C(q) %for other clusters
%             Nj = wcl(j,:);
%             wCT(i,j) = sum(min(Ni,Nj));
%         end
%     end
% end

maxWCT = max(wCT(:));
if maxWCT > 0 
    wCT = wCT / maxWCT;
end
wCT = wCT + wCT';
wCT(1:nCls+1:nCls^2) = 1;


%---find pair-wise similarity of data points--------------------------------
% tic();
% S2 = zeros(n,n); % create matrix S, pair-wise similairty matrix for each pair of data points
% for m = 1:M
%     for i = 1:n-1 % for each row, start at row #1 
%         ii = i+1:n; % for other rows (below i)
%         indEq = pc(ii,E(i,m));%E(i,m) == E(ii,m);
%         %indEq = bsxfun(@eq,E(i,:),E(ii,:));
%         
%         tmp = ones(1,length(ii));
%         tmp(~indEq) = wCT(E(i,m),E(ii(~indEq),m))*dc;
%         
%         S2(i,ii) = S2(i,ii) + tmp;
%     end
% end
% t2 = toc();
% 
% tic();
S = zeros(n,n);
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
% t1 = toc();
% 
% fprintf(1,'\n%d | %f, %f ---------------------------------\n', isequal(S,S2),t1,t2);

S = S/M;
S = S + S';
S(1:n+1:n^2) = 1;

