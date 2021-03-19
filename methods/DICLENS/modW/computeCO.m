function CO = computeCO(P,numPartitions,inMode,outMode,weights)
% CO = computeCO(P,numPartitions,inMode,outMode,weights) 
% computes co-occurence accumulation matrix from binary representations or 
% integer labels of clusterings in the ensemble
%--------------------------------------------------------------------------
% INPUTS
% P:            - binary representation of ensemble (nClusters X N)
%               - ensemble labels (N X numPartitions)
%
% numPartitions: number of partitions (size of the ensemble). Can be empty
%                for ensemble labels.
%
% inMode:       - 'binary': input P is binary matrix (default)
%               - 'labels': input P is matrix with labels of ensemble (each
%                           column in a partition)
%
% outMode:      - 'full'  : output is N X N matrix (default)
%               - 'vec'   : output is N*(N-1)/2 long vector of pairs
%               - 'sparse': output is sparse matrix
% weights:      - vector of weights for each partition if inMode is labels; 
%               - vector of weights for each cluster if inMode is binary;
%               - if empty, we will use a vector of ones.
%
% OUTPUTS
% CO:           co-occurence evidence accumulation values in a form
%               determined by outMode parameter
%--------------------------------------------------------------------------
% EXAMPLE
%
% % generate random labels; 1000 data points, 10 partitions with max 3
% % clusters in each 
% P = randi(3,1000,10);
% 
% % Compute co-occurence matrix in two different ways:
% % A. Use labels
% CO_A = computeCO(P);
% 
% % B. First compute BA matrix and then use it
% BA = computeBA(P);
% CO_B = computeCO(BA,10,'binary');
% 
% % Test for equality
% isequal(CO_A, CO_B)
%
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 11-July-2013 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if ~exist('inMode','var') || isempty(inMode)
    inMode = 'labels';
end

if ~exist('outMode','var') || isempty(outMode)
    outMode = 'full';
end

switch inMode
    case 'binary'
        compareFun = @and; % for binary comparison
        if ~exist('numPartitions','var') || isempty(numPartitions)
            error('Parameter numPartitions is not defined!');
        end
        M = numPartitions;
        wSize = size(P,1);
        
    case 'labels'
        compareFun = @eq; % for integer comparison
        P = P';
        if ~exist('numPartitions','var') || isempty(numPartitions)
            M = size(P,1);
        else
            M = numPartitions;
        end
        wSize = M;
        
    otherwise
        error('Wrong inMode parameter!');
end

N = size(P,2);

if ~exist('weights','var') || isempty(weights)
    weights = ones(wSize,1);
end



switch outMode
    case 'full'
        CO = zeros(N,N);
        
        for i=1:N-1
            m = bsxfun(compareFun,P,P(:,i));
            m = bsxfun(@times,m,weights);
            CO(i,:)= sum(m,1)./M;
        end
        CO = triu(CO,1); % upper part of the matrix without diagonal
        CO = CO + CO';
        
    case 'vec'
        nPairs = N*(N-1)/2;
        CO = zeros(1,nPairs);
        ind=[0 cumsum(N-1:-1:1)];
        
        for i=1:N-1
            m = bsxfun(compareFun,P,P(:,i));
            m = bsxfun(@times,m,weights);
            s = sum(m,1)./M;
            CO( ind(i)+1 : ind(i+1)) = s(i+1:end);
        end        
        
    case 'sparse'
        % try if sparse matrix would fit into memory
        AMvec = zeros(1,N*(N-1)/2);
        V = nchoose2(1:N);
        I = V(:,1)';
        J = V(:,2)';
        ind=[0 cumsum(N-1:-1:1)];
        
        for i=1:N-1
            m = bsxfun(compareFun,P,P(:,i));
            m = bsxfun(@times,m,weights);
            s = sum(m,1)./M;
            AMvec( ind(i)+1 : ind(i+1)) = s(i+1:end);
        end
        CO = sparse(I,J,AMvec,N,N);
        CO = CO + CO';
        
    otherwise
        error('Wrong outMode value!');
end