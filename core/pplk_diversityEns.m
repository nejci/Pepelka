function [dvrsMean, dvrsSpread, dvrsComb1, dvrsComb2, accu] = pplk_diversityEns(labelsEns,labelsRef,method,measures)
% [dvrsMean, dvrsSpread, dvrsComb1, dvrsComb2, accu] = pplk_diversityEns(labelsEns,labelsRef,method,measures)
% Computes diversity measures from partitions in the ensemble.
%
% INPUTS		
%   labelsEns    
%       Cluster ensemble - labels are stored in the columns of an N-by-E
%       matrix, where N is number of data points and E is number of
%       clusterings.
%
%   labelsRef    
%       Reference partitioning, relevant only for methods 'NPW'.
%
%   method       
%       PW
%           Pairwise.
%       NPW 
%           Non-pairwise, required labelsRef.
%
%   measures 
%       External validation indeces names (one or more).
%
%
% OUTPUTS
%   dvrsMean     
%       Mean of PW or NPW diversity measure across ensemble members.
%
%   dvrsSpread   
%       Sum of standard deviation of NPW diversity.
%
%   dvrsComb1    
%       Combined measure: 1/2*(1-dvrsMean + dvrsSpread).
%
%   dvrsComb2    
%       Combined measure: dvrsSpread/dvrsMean.
%
%   accu         
%       Accuracy of ensemble members.
%
%
% ACKNOWLEDGEMENTS AND REFERENCES
%   X. Fern and C. Brodley, "Random projection for high dimensional
%     data clustering: A cluster ensemble approach," in International
%     Conference on Machine Learning, 2003. 
%   S. T. Hadjitodorov, L. I. Kuncheva, and L. P. Todorova, "Moderate
%     diversity for better cluster ensembles," Information Fusion, vol. 7,
%     no. 3, pp. 264-275, Sep. 2006.
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka

callDir=chdir(pplk_homeDir());

[N,E] = size(labelsEns);

if ~exist('measures','var') || isempty(measures)
    measures = 'NMI';
end
if ~iscell(measures)
   measures = {measures}; 
end
M = length(measures);

if ~exist('method','var') || isempty(method)
    method = 'PW';
end
if strcmpi(method,'NPW')
    if ~exist('labelsRef','var') || isempty(labelsRef)
        error('labelsRef must be passed to compute non-pairwise diversity.');
    else
        [r c] = size(labelsRef);
        if r == 1 && c > 1
            labelsRef = labelsRef';
            r = c;
            c = 1;
        end
        
        assert(c == 1,'Reference partitioning must have only 1 column.');
        assert(r == N,'Reference partitioning must have the same number of points as ensemble members.');
    end
end




if nargout > 1 && strcmpi(method,'PW')
    error('Cannot assign value to output arguments other than dvrsMean for method PW.');
end

switch method
    case 'PW'
        % Compute accuracy measure between each pair of labels
        Smat = zeros(E,E,M);
        for i=1:E-1
            for j=i+1:E
                [~,list] = pplk_validExt(labelsEns(:,i),labelsEns(:,j),measures);
                Smat(i,j,:) = 1-list;
            end
        end
        D = squeeze(sum(sum(Smat,1),2))';
        dvrsMean = D / ((E^2-E)*0.5);
        
    case 'NPW'
        % Compute accuracy measure between each ensemble member and
        % reference partitioning
        Smat = zeros(E,M);
        for i=1:E
            [~,list] = pplk_validExt(labelsRef,labelsEns(:,i),measures);
            Smat(i,:) = 1-list;            
        end
        dvrsMean = mean(Smat,1);
        dvrsSpread = sqrt(sum(bsxfun(@minus,Smat,dvrsMean).^2,1)/(E-1)); %std(Smat,0,1);
        dvrsComb1 = 0.5 * (1-dvrsMean + dvrsSpread);
        dvrsComb2 = max(0,dvrsSpread./dvrsMean);
        accu = 1-Smat;
        
end



chdir(callDir);