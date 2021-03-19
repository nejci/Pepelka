pool_init(6);

parfor i = 1:6
    
    t = getCurrentTask();
    ID = t.ID;
    fprintf(1,'%d\n',ID);
    
    labelsEns2 = [1 1 2 2; 1 2 1 2]';
    K = 2;
    labelsCons = pplk_consEns(labelsEns2,K,'STREHL-CSPA',[]);
    
end