% CDF of the sum of 'n' independent random variables which are distributed
% according to Beta distributions.
% 
% Note: When computing multiple values at once, 'betaconv' contains a more
% efficient implementation.
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
% $Id: betasumcdf.m 8246 2010-10-22 13:28:23Z bkay $
% -------------------------------------------------------------------------
function y = betasumcdf(x, AB, res)
    
    if ~(ismatrix(x) && (size(x,1)==1 || size(x,2)==1))
        error('only implemented for onedimensional input');
    end
    
    % Compute the PDF first (since we want the entire pdf rather than just
    % one value from it, using betaChfSum is computationally more efficient
    % than using betasumpdf)
    c = betaChfSum(res, AB);
    
    % Sum the PDF up to point x
    for i=1:length(x)
        idx = round(x(i)/res);
        if idx < 1
            y(i) = 0;
        elseif idx > length(c)
            y(i) = 1;
        else
            y(i) = trapz(c(1:idx)) * res;
        end
    end
    
end
