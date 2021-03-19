function [aggR, pval, rowNames] = aggregateRanks(R, N, method, complete, topCutoff)
%AGGREGATERANKS Aggregate ranked lists using traditional and robust methods
%   [aggR, pval, rowNames] = AGGREGATERANKS(R,N,method,complete,topCutoff) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
%   INPUTS:
%       R           -> matrix representation: numeric matrix with as many
%                   rows as the number of unique elements and with as many
%                   columns as the number of ranked lists. All entries have
%                   to be on interval [0,1]. Smaller numbers correspond to
%                   better ranking. If R is a column vector, nothing
%                   happens - it is returned in aggR.%                   
%                   -> list representation: cell array of
%                   cells with strings or vectors with numbers - in this
%                   case R is transformed into numeric rank matrix.
%       N           number of ranked elements, default is the number of
%                   unique elements in R
%       method      rank aggregation method. Could be one of the following:
%                   'min', 'median', 'mean', 'geom.mean', 'stuart', 'RRA'.
%                   Default is 'RRA' (Robust Rank Aggregation).
%       complete    1 - rankings are complete (there is perfect match
%                   between sets of rankings) 
%                   0 - default; rankings are incomplete.
%       topCutoff   vector of cutoff values that limit the number of
%                   elements in the input list.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       aggR        vector of aggregated ranks/scores (it equals pval for 
%                   methods 'stuart' and'RRA')
%       pval        p-values (relevant only for 'mean','stuart','RRA')
%       rowNames    if R contains lists, rowNames contains their unique
%                   names in the same order as the values of aggR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   EXAMPLES
%       % Lets have three ordered lists of names.
%       R = {   {'Joe', 'Bob', 'Lucy', 'Mary'}, ...
%               {'Bob', 'Joe', 'Lucy', 'Mary'}, ...
%               {'Lucy', 'Mary', 'Bob', 'Joe'}}
%
%       % We can also use numerical vectors instead of cell of strings.
%       % R = { [1,2,3,4], [2,1,3,4], [3,4,2,1] }
%       
%       % Obtain aggregated ranks with method 'RRA' (default).
%       [aggR, pval, rowNames] = aggregateRanks(R)
%
%       % Or, equivalently, use explicit parameters definition.
%       [aggR, pval, rowNames] = aggregateRanks(R, [], 'RRA')
%
%       % We can also compute a matrix with ranks first ...
%       rmat = rankMatrix(R)
%       % ... and then pass it to the aggregation method.
%       [aggR, pval, rowNames] = aggregateRanks(rmat)
%
%       % A case of incomplete lists.
%       R = {   {'Joe', 'Bob', 'Lucy', 'Mary'}, ...
%               {'Bob', 'Joe', 'Lucy',       }, ...
%               {'Lucy', 'Mary'              }}
%
%       % Lets compute mean ranks. Mind the fourth parameter, which
%       % indicates completeness of the lists. Note the return values; aggR
%       % contains average across the ranks, while pval contains the
%       % statistical significance (p-value) of mean ranks. 
%       [aggR, pval, rowNames] = aggregateRanks(R,[],'mean',0)
%
%       % We can also say that only top k elements are presented in data 
%       % by setting the topCutoff to [1,0.75,0.5].
%       [aggR, pval, rowNames] = aggregateRanks(R,[],'RRA',0,[1,0.75,0.5])
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   See also RANKMATRIX.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@fri.uni-lj.si> 
%   Based on GNU R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------

    rowNames = [];
    pval = NaN;
    
    if ~exist('R','var') || isempty(R)
        error('Input parameter R is missing!');
    end
    
    if ~exist('N','var')
        N=[];
    end
    
    if ~exist('complete','var') || isempty(complete)
        complete = 0;
    end

    % Input parameter R determination
    if iscell(R)
        [rmat, rowNames] = rankMatrix(R, N, complete);
    elseif ismatrix(R)
        if all(max(R,[],1)<=1) && all(min(R,[],1)>0)
            rmat = R;            
        else
            error('Columns of matrix R can only contain numbers from interval (0,1].');
        end
    else
        error('R should be cell (of lists) or matrix (of ranks).');
    end
    
    % if there is only one ranking, nothing can be done, return it
    if isvector(rmat)
        aggR = rmat;
        return;
    end
    
    
    if ~exist('method','var') || isempty(method)
        method = 'RRA';
    end

    if ~exist('topCutoff','var') || isempty(topCutoff)
        topCutoff = NaN;
    end   

    switch lower(method)

        case 'min'
            aggR = min(rmat,[],2);
            
        case 'median'
            aggR = nanmedian(rmat,2);
            
        case 'geom.mean'
            aggR = exp(nanmean(log(rmat),2));

        case 'mean'
            aggR = nanmean(rmat, 2);
            n = sum(~isnan(rmat),2);
            pval = normcdf(aggR, 0.5, sqrt(1/12./n));

        case 'stuart'
            aggR = stuart(rmat);
            pval = aggR;

        case 'rra'
            aggR = rhoScores(rmat, topCutoff);
            pval = aggR;

        otherwise
            error('Method should be one of:  "min", "geom.mean", "mean", "median", "stuart" or "RRA"');
    end
end

%--------------------------------------------------------------------------
% aggregateRanks helper functions

function [rmat, rowNames] = rankMatrix(glist, N, complete)
%RANKMATRIX Transform ranked lists into rank matrix
%   [rmat, rowNames] = RANKMATRIX(glist, N, complete) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
%   INPUTS:
%       glist       cell array of cells with strings or numerical vectors.
%       N           number of ranked elements, default is the number of
%                   unique elements in R
%       complete    1 - rankings are complete (there is perfect match
%                   between sets of rankings) 
%                   0 - default; rankings are incomplete.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       rmat        numeric matrix with as many rows as the number of
%                   unique elements and with as many columns as the number 
%                   of ranked lists. All entries have to be on interval 
%                   [0,1]. Smaller numbers correspond to better ranking.
%       rowNames    if R contains lists, rowNames contains their unique
%                   names in the same order as the values of aggR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   EXAMPLES
%       % Lets have three ordered lists of names.
%       R = {   {'Joe', 'Bob', 'Lucy', 'Mary'}, ...
%               {'Bob', 'Joe', 'Lucy', 'Mary'}, ...
%               {'Lucy', 'Mary', 'Bob', 'Joe'}}
%
%       % We can also use numerical vectors instead of cell of strings.
%       % R = { [1,2,3,4], [2,1,3,4], [3,4,2,1] }
%
%       % Compute a matrix with ranks.
%       rmat = rankMatrix(R)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   See also AGGREGATERANKS.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@fri.uni-lj.si> 
%   Based on GNU R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------
    
    if ~exist('complete','var') || isempty(complete)
        complete = 0;
    end
    
    nCell = length(glist);
    NperCell = cellfun(@length, glist);
    
    % convert each cell in glist to row vector
    glist = cellfun(@(c) c(:)', glist, 'UniformOutput',false);
    glist_all = [glist{:}];
  
    % find unique elements and map strings to numbers (ib_map)
    [u_map, ~, ib_map]= unique(glist_all);
    realLen = length(u_map);
    
    if ~exist('N','var') || isempty(N)
        N = realLen;
    end
    
    % complete each list with NaNs to match its length with realLen
    M = nan(realLen, nCell);
    start = 1;
    for c=1:nCell
        stop = start + NperCell(c)-1;
        M(1:NperCell(c),c) = ib_map(start:stop);
        start = stop + 1;
    end

    [U,~,ib] = unique(M);
    fullLen = length(U);
    
    if ~complete
        rmat = ones(fullLen, nCell)*N;
        N = ones(1,nCell)*N;
    else
        rmat = nan(fullLen, nCell);
        N = NperCell;
    end

    Umat = reshape(ib,realLen,nCell);
    v = (1:realLen)';
    
    for i=1:nCell
        u=Umat(:,i);
        rmat(u,i) = v;
    end
    
    rmat = rmat(1:realLen,:);
    rmat = bsxfun(@rdivide,rmat,N);

    rowNames = u_map(1:realLen)';
end

function rho = rhoScores(r, topCutoff)
%RHOSCORES Compute Rho scores for rank vector
%   rho = RHOSCORES(r, topCutoff) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
%   INPUTS:
%       r           vector of normalized rank values on interval [0,1]
%       topCutoff   a vector of cutoff values used to limit the number of 
%                   elements in the input lists
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       rho         a vector of rho values, corrected against bias from
%                   multiple hypothesis testing (Bonferroni correction).
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   See also CORRECTBETAPVALUES, THRESHOLDBETASCORE, AGGREGATERANKS.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@fri.uni-lj.si> 
%   Based on GNU R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------
    if ~exist('topCutoff','var') || isempty(topCutoff)
        topCutoff = NaN;
    end
    
    if isvector(r)
        rows = 1;
        r = r(:)'; % force row vector form
    else
        rows = size(r,1);
    end
    
    rho = nan(rows,1);
    
    for rInd = 1:rows
        r1 = r(rInd,:);
        
        if(isnan(topCutoff(1)))
            x = betaScores(r1);
            % Correct using Bonferroni method.
            rho(rInd) = correctBetaPvalues( min(x), sum(~isnan(x)));
        else            
            r1 = r1(~isnan(r1));
            r1(r1 == 1) = nan;
            % Consider thresholds in topCutoff vector.
            x = thresholdBetaScore(r1,[],[],topCutoff);
            % Correct using Bonferroni method.
            rho(rInd) = correctBetaPvalues( min(x), length(r1));
        end
    end
end

function pval = correctBetaPvalues(p,k)
%CORRECTBETAPVALUES Compute p-values based on Beta distribution
%   pval = CORRECTBETAPVALUES(p,k) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@fri.uni-lj.si> 
%   Based on GNU R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------
    pval = betacdf(p,1,k);
end

function [Beta,names] = thresholdBetaScore(r, k, n, sigma)
%THRESHOLDBETASCORE Compute p-values based on Beta distribution
%   [Beta,names] = THRESHOLDBETASCORE(r, k, n, sigma) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@fri.uni-lj.si> 
%   Based on GNU R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------

% Input variables validation
    rLen = length(r);

    if ~exist('k','var') || isempty(k)
        k = 1:rLen;
    end
    if ~exist('n','var') || isempty(n)
        n = rLen;
    end
    if ~exist('sigma','var') || isempty(sigma)
        sigma = ones(1,n);
    end

    if(length(sigma) ~= n)
        error('The length of sigma does not match n!');
    end
    if(length(r) ~= n)
        error('The length of p-values does not match n!');
    end
    if(min(sigma)< 0 || max(sigma) > 1)
        error('Elements of sigma are not in the range [0,1]!');
    end
    if(any(~isnan(r) & r > sigma))
        error('Elements of r must be smaller than elements of sigma!');
    end
%--------------------------------------------------------------------------

    x = sort(r(~isnan(r)));
    sigma = sort(sigma, 'descend');
    Beta = nan(1, length(k));

    for i = 1:length(k)

        if(k(i) > n)
            Beta(i) = 0;
            continue;
        end
        if(k(i) > length(x))
            Beta(i) = 1;
            continue;
        end
        if(sigma(n) >= x(k(i)))
            Beta(i) = betacdf( x(k(i)), k(i), n + 1 - k(i));
            continue;
        end
        
        % Non-trivial cases
        % Find the last element such that sigma(n0) <= x(k(i))
        n0 = find(sigma < x(k(i)));
        n0 = n0(1) - 1;

        % Compute beta score vector beta(n,k) for n = n0 and k = 1..k(i)
        if(n0 == 0)
            B = [1, zeros(1, k(i))];
        elseif(k(i) > n0)
            B = [1, betacdf(x(k(i)), 1:n0, n0:-1:1), zeros(1, k(i) - n0)];
        else
            B = [1, betacdf(x(k(i)), 1:k(i), n0+1-(1:k(i)) )];
        end

        % In the following update steps sigma < x(k(i))
        z = sigma( (n0+1) : n );
        for j = 1:(n - n0)
            B( 2:(k(i)+1)) = (1-z(j)) * B(2:(k(i)+1)) + z(j) * B(1:k(i));
        end

        Beta(i) = B(k(i)+1);
    end

    names = k;

end

function p = betaScores(r)
%BETASCORES Compute Beta scores for rank vector
%   p = BETASCORES(r) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
%   INPUTS:
%       r   vector of normalized rank values on interval [0,1]
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       p   a vector of p-values that corresponds to the sorted input 
%           vector. The NaN-s are moved to the end.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   See also CORRECTBETAPVALUES, THRESHOLDBETASCORE, AGGREGATERANKS.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@fri.uni-lj.si> 
%   Based on GNU R package RobustRankAggreg written by Raivo Kolde. 
%   Reference:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012).
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------
    n = sum(~isnan(r));
    p = nan(1,length(r));
    % Sort the values.
    r = sort(r);
    % Get the order statistics and calculates p-values for each of the
    % order statistics. These are based on their expected distribution
    % under the null hypothesis of uniform distribution.
    p(1:n) = betacdf(r(1:n),1:n,n:-1:1);
end

function aggR = stuart(rmat)
%STUART Compute aggregated rank with Stuart-Aerts method
%   aggR = STUART(rmat) 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
%   INPUTS:
%       rmat    numeric matrix with as many rows as the number of
%               unique elements and with as many columns as the number 
%               of ranked lists. All entries have to be on interval 
%               [0,1]. Smaller numbers correspond to better ranking. 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   OUTPUTS:
%       aggR    vector of aggregated ranks (p-values)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   See also    AGGREGATERANKS.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Copyright (2013) Nejc Ilc <nejc.ilc@fri.uni-lj.si> 
%   Based on GNU R package RobustRankAggreg written by Raivo Kolde. 
%   References:
%     Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012). 
%     Robust rank aggregation for gene list integration and meta-analysis.
%     Bioinformatics, 28(4), 573-580
%
%     Stuart, J. M., Segal, E., Koller, D., & Kim, S. K. (2003). 
%     A gene-coexpression network for global discovery of conserved genetic
%     modules. Science, 302(5643), 249-55
%
%     Aerts, S., Lambrechts, D., Maity, S., Van Loo, P., Coessens, B., De
%     Smet, F., Tranchevent, L.-C., et al. (2006). Gene prioritization
%     through genomic data fusion. Nature biotechnology, 24(5), 537-44
%   
%   Revision: 1.0 Date: 2013/05/16
%--------------------------------------------------------------------------
	rmat = sort(rmat, 2);
    N=size(rmat,1);
	aggR = zeros(N,1);
    for ai = 1:N
       aggR(ai) = qStuart(rmat(ai,:)); 
    end
end

function q=qStuart(r)
% Stuart-Aerts method helper function
	N = sum(~isnan(r));
	v = ones(1, N+1);
    for k = 1:N
        v(k+1) = sumStuart( v(1:k), r(N-k+1));
    end
	q = factorial(N) * v(N+1);
end

function s = sumStuart(v, r)
% Stuart-Aerts method helper function
	k = length(v);
	l_k = 1:k;
	ones = (-1).^(l_k + 1);
	f = factorial(l_k);
	p = r.^l_k;
	s = ones * ( v(end:-1:1) .* p ./ f)';
end

%--------------------------------------------------------------------------
% Auxiliary functions NANMEAN, NANMEDIAN from Statistics Toolbox

function m = nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:15:50 $

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end
end

function y = nanmedian(x,dim)
%NANMEDIAN Median value, ignoring NaNs.
%   M = NANMEDIAN(X) returns the sample median of X, treating NaNs as
%   missing values.  For vector input, M is the median value of the non-NaN
%   elements in X.  For matrix input, M is a row vector containing the
%   median value of non-NaN elements in each column.  For N-D arrays,
%   NANMEDIAN operates along the first non-singleton dimension.
%
%   NANMEDIAN(X,DIM) takes the median along the dimension DIM of X.
%
%   See also MEDIAN, NANMEAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:15:51 $

if nargin == 1
    y = prctile(x, 50);
else
    y = prctile(x, 50,dim);
end
end

