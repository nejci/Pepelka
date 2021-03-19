
labelsEns = [3 1 2; 2 2 1]';
Kcons = 2;
[labelsCons,numClust]= pplk_consEns(labelsEns,Kcons,'LCE2-SRS-SL',params);



% 
% pc = logical(pc);
% S2 = zeros(n,n); % create matrix S, pair-wise similairty matrix for each pair of data points
% 
%     for i = 1:n-1 % for each row, start at row #1 
%         ii = i+1:n; % for other rows (below i)
%         %indEq = pc(ii,E(i,m));%E(i,m) == E(ii,m);
%         indEq = pc(ii,E(i,:));%bsxfun(@eq,E(i,:),E(ii,:));
%         
%         tmp = ones(length(ii),M);
%         for m = 1:M
%             tmp(~indEq(:,m),m) = wCT(pc(i,:),E(ii(~indEq(:,m)),m))*dc;
%         end
%         S2(i,ii) = S2(i,ii) + tmp;
%     end