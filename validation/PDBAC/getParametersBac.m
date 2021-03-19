%This function get the parameters of the beta distrbution from all the
%classes present in the confusion matrix
%
% Usage:
%     getParametersBac(C)
%
% Arguments:
%     C - confusion matrix of classification outcomes
% Return:
%     rslt - Matrix with two column, containg alpha and beta of the beta
%     distribution of the classes, and as many rows as classes the
%     confussion matrix contains.
%
% Henry Carrillo, University of Zaragoza, Spain
% http://www.hcarrillo.com/
% Vectorized by Nejc Ilc
% -------------------------------------------------------------------------
function rslt = getParametersBac(C)
    
if( size(C,1) ~= size(C,2) )
    error('Confusion matrix is not square!!!')
end
numberClasses = size(C,1);
TP = diag(C);
cw = C.*~eye(numberClasses);
FN = sum(cw,1)';
rslt = [TP+1,FN+1];
end