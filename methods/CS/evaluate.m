function [numErr]=evaluate(labels,target,K_true)

%sedaj je potrebno oceniti rezultat - koliko je bilo napak. Primerjamo
%pravo razdelitev (target) s programsko
%PROBLEM: poimenovanja skupin niso ista, torej skupina 1 v target ni nujno
%enaka skupini 1 v rešitvi s programom. Zato privzamemo, da je v programski
%rešitvi neka skupina doloèena z veèino pripadnosti.
%Primer: [2 2 1 2 3 2 2] doloèa skupino 2, napaki sta 2
%
%Avtor: Nejc Ilc

numErr=0;
avail=1:K_true;

for i=1:K_true
    
    j=1;
    L=[];
    
    tmp=labels(target==i);
    
    vec=(tmp==tmp(1));
    L(j,1)=tmp(1);
    L(j,2)=sum(vec);
    tmp(vec)=[];
    
    while length(tmp)>0
        j=j+1;
        vec=(tmp==tmp(1));
        L(j,1)=tmp(1);
        L(j,2)=sum(vec);
        tmp(vec)=[];        
    end
    [val,maxInd]=max(L(:,2));
    
    if(sum(L(maxInd,1)==avail)>0)
        avail(L(maxInd,1)==avail)=[];
    
        L(maxInd,:)=[];
    end
    numErr=numErr+sum(L(:,2));

end

    
end