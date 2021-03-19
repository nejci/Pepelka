function PCCA = PerClassAcc(C)
% Per-class accuracy (balanced classification accuracy)

a = diag(C)./sum(C,2);
a(isnan(a))=0;
PCCA = (1/size(C,1))*sum(a);
    
end