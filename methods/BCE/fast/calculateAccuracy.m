function [micro_precision]=calculateAccuracy(true_labels,Ensemble_labels)
%
% Calculate  micro-precision given clustering results and true labels.
%
%   k = number of ensemble clusters
%   M = number of data points
%
% Input:
%   true_labels:        1*M, true class labels for the data points
%   Ensemble_labels:    1*M, labels obtained from BCE
    
% Output:
%   micro_precision:    micro-precision
%--------------------------------------------------------------------

k=length(unique(true_labels));
M=length(true_labels);
   
    for j=1:k
         for jj=1:k
            [xx,accurence(j,jj)]=size(find(((Ensemble_labels==jj)*j)==true_labels));
         end
    end 

    [rowm,coln]=size(accurence);
    amatrix=accurence;
    sumMax=0;
    while rowm>=1
        xx=max(max(amatrix));
        [x,y]=find(amatrix==xx);
        sumMax=sumMax+xx;                      
        iyy=1;
        temp=zeros(rowm,rowm-1);
        for iy=1:rowm
            if iy==y(1)
                continue;
            else                        
                temp(:,iyy)=amatrix(:,iy);
                iyy=iyy+1;
             end
         end  
         temp2=zeros(rowm-1,rowm-1);
         ixx=1;
         for ix=1:rowm
            if ix==x(1)
                continue;
            else                        
                temp2(ixx,:)=temp(ix,:);
                ixx=ixx+1;
             end
          end
          rowm=rowm-1;
          amatrix=zeros(rowm,rowm);
          amatrix=temp2;
                  
    end

   micro_precision=sumMax/M;
