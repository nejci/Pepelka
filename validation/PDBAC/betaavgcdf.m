% CDF of the average of 'n' independent random variables which are
% distributed according to Beta distributions.
%
% Literature:
%     K.H. Brodersen, C.S. Ong, K.E. Stephan, J.M. Buhmann (2010).
%     The balanced accuracy and its posterior distribution. In: Proceedings
%     of the 20th International Conference on Pattern Recognition.
%
% Edited by:
% Henry Carrillo, University of Zaragoza, Spain
% http://www.hcarrillo.com/
%
% Original coder:
% Kay H. Brodersen, ETH Zurich, Switzerland
% http://people.inf.ethz.ch/bkay/
% $Id: betaavgcdf.m 8245 2010-10-22 12:57:51Z bkay $
% -------------------------------------------------------------------------
function y = betaavgcdf(x, AB, res)
    
    y = betasumcdf(size(AB,1)*x, AB, res);
    
end
