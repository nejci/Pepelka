function NCAI = getNCAI(baseCls)
% Huang Dong. 14 May, 2014.
% Get the normalized crowd agreement index.
% Optimized for speed by Nejc Ilc, 5 June 2014

cnt = size(baseCls,2);
ClsAfnyWrtNMI = zeros(cnt, cnt); % it is symmetric.
for i = 1 : cnt
    for j = i+1 : cnt
        ClsAfnyWrtNMI(i,j) = NMImax(baseCls(:,i), baseCls(:,j));
    end
end

ClsAfnyWrtNMI = ClsAfnyWrtNMI + ClsAfnyWrtNMI';

CAI = (sum(ClsAfnyWrtNMI)./(cnt-1))';
NCAI = CAI./max(CAI(:));