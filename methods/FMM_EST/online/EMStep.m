function [mixEM,meansEM,varsEM]=EMStep(mixEM,meansEM,varsEM,dataNew,fAlpha)
%EMGaussianMixture biased OnLine
% [mix,means,vars]=EMStep(mix,means,vars,dataNew,fAlpha) 
% Recursive Unsupervised Learning of Finite Mixture Models
%
% Input:
% - current estimate:
%   mix     - 1 x nModes matrix
%   means   - nDimensions x nModes matrix
%   vars    - nDimensions x nDimensions x nModes array
%
% - new data sample:
%   dataNew - nDimensions x 1 vector
%
% - parameter:
%   fAlpha  - exponential envelope to limit the influence of the old data
% 
% Output:
% - updated estimates
%
% Example:
%   >> dataSet1;%creates data set: dataEM - nDimensions x nData matrix
%   >> nData=size(dataEM,2);
%   >> [mixEM,meansEM,varsEM]=EMRandomInit(dataEM,30);%start with 30 components
%   >> for iData=1:nData
%       [mixEM,meansEM,varsEM]=EMStep(mixEM,meansEM,varsEM,dataEM(:,iData),1/150);
%      end
%   >> PlotData(dataEM(:,nData-150:nData));
%   >> hold on
%   >> PlotGM(mixEM,meansEM,varsEM);
%
% Author:   Z.Zivkovic
% Date:     10-2-2003


nD=size(dataNew,1);
nModes=length(mixEM);
nModesUsed=length(find(mixEM));
%parameters
cPrune=-0.5*(nD+nD*(nD+1)/2)*fAlpha;
        
%%%%%%%%%%%
%algorithm - update step
%
%input:dataNew,mixEM,meansEM,varsEM
        
    
    %calculate P and Pi for the current data sample
    pdfPtot=0;
    pdfPi=zeros(1,nModes);
    for iModes=1:nModes
        if (mixEM(iModes)~=0)
        %gauss
        diff=(dataNew-meansEM(:,iModes));
        dist=0.5*diff'*inv(varsEM(:,:,iModes))*diff;
        cEM=(1/((2*pi)^(nD/2)*sqrt(det(varsEM(:,:,iModes)))));
        pdfPi(iModes)=mixEM(iModes)*cEM*exp(-dist);%!!!!mix already included !!!
        pdfPtot=pdfPtot+pdfPi(iModes);
        end
    end
        
    %update parameters
    for iModes=1:nModes
       if (mixEM(iModes)~=0)
        own=pdfPi(iModes)/pdfPtot;
        k=fAlpha*own/mixEM(iModes);
        %mixEM(iModes)=mixEM(iModes)+fAlfa*(own-mixEM(iModes))+fAlfa*cPrune;
        mixEM(iModes)=mixEM(iModes)+fAlpha*(own/(1+nModesUsed*cPrune)-mixEM(iModes))+fAlpha*cPrune/(1+nModesUsed*cPrune);
        %mixEM(iModes)=mixEM(iModes)+fAlfa*(own-mixEM(iModes))+fAlpha*cPrune;%try this
        
        if (mixEM(iModes)<0)%-fAlpha*cPrune/(1+nModesUsed*cPrune))%
            mixEM(iModes)=0;
            mixEM=mixEM/(sum(mixEM));%normalize
        else
            diff=dataNew-meansEM(:,iModes); 
            meansEM(:,iModes)=meansEM(:,iModes)+k*diff;
            if (k>20*fAlpha)%!!!! - to prevent the cov matrix from beeing singular
                k=20*fAlpha;
            end%if own was 1 this suppress update of the mixtures representing less than 1/20 of the data
            varsEM(:,:,iModes)=varsEM(:,:,iModes)+k*(diff*diff'-varsEM(:,:,iModes));
        end
       end
    end
    
    %output: new estimates: mixEM,meansEM,varsEM
    %
    %algorithm
    %%%%%%%%%%%


    