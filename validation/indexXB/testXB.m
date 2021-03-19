% benchmark test

load test.mat

Kmax = size(centers,2)+1;
t1=0;
t2=0;
t3=0;
t4=0;

for iter = 1:1000
    
    tic;
    centers = cell(1,Kmax-1);
    U = cell(1,Kmax-1);
    for i=1:Kmax-1
        lbl = labels(:,i);
        E = eye(i+1);
        U{i} = logical(E(:,lbl));
        clsize = sum(U{i},2);
        centers{i} = bsxfun(@rdivide, U{i}*data, clsize);
    end 
    tC = toc;
    t2=t2+tC;
    t4=t4+tC;
    
    for k=1:Kmax-1
        
        tic;
        [XB1, XBmod1] = indexXB(data, labels(:,k));
        t1 = t1 + toc;
        
        tic;
        lblStruct.U = U{k};
        lblStruct.center = centers{k};
        [XB2, XBmod2] = indexXB(data, lblStruct);
        t2 = t2 + toc;
        
        tic;
        [XB3, XBmod3] = indexXB2(data, labels(:,k));
        t3 = t3 + toc;
        
        tic;
        lblStruct.U = U{k};
        lblStruct.center = centers{k};
        [XB4, XBmod4] = indexXB2(data, lblStruct);
        t4 = t4 + toc;
    end
end