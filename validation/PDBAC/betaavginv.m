% Inverse CDF of the sum of 'n' independent random variables which are
% distributed according to Beta distributions.
%
% If the optimization fails, function returns NaN.
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
% $Id: betaavginv.m 8245 2010-10-22 12:57:51Z bkay $
% -------------------------------------------------------------------------
function x = betaavginv(y, AB, res)
    
    try
        s=size(AB,1);
        x = fzero(@(z) betasumcdf(s*z,AB,res)-y, 0.5);
               
    catch err
        disp(['Error occurred in BETAAVGINV: ', err.message]);
        x = NaN;
    end
    
end
