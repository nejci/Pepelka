% PDF of the sum of 'n' independently distributed Beta distributions.
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
% $Id: betasumpdf.m 8245 2010-10-22 12:57:51Z bkay $
% -------------------------------------------------------------------------
function y = betasumpdf(x, AB, res) 
    % Compute the sum using the characteristic function
    c = betaChfSum(res, AB);
    
    % Prepare return value
    y = NaN(size(x));
    
    % Fill in return value
    % Set support
    supportPdf = size(AB,1);
    % - values outside support
    y(x<0 | x>supportPdf) = 0;
    % - values 
    % - all other values
    idx = int32(x/res+1);
    y(isnan(y)) = c(idx(isnan(y)));    
end
