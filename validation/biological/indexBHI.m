function BHI = indexBHI(genenames, labels, annotations, aspect, ignoreEvidence)

% Function computes Biological Homogenity Index (Datta and Datta 2006)
%-------------------------------------------------------------------------
% INPUTS
%   genenames       (cell)     strings with gene identifiers (must match 
%                              with identifiers in annotations)
%
%   labels          (vector)   cluster labels
%
%   annotations     (struct)   structure array with fields: 
%                               .geneID         (gene identifier, e.g., Affymetrix probe ID)
%                               .functionID     (functional class identifier, e.g., GO ID)
%                               .aspect         (optional - when used with GO)
%                               .evidence       (optional - when used with GO)
%
%   aspect          (cell)      which aspects to include:
%                               'BP' | 'CC' | 'MF' (or any subset of them); 
%                               leave empty to compute all of them.
%
%   ignoreEvidence  (cell)      which evidence to ignore evidence code,
%                               i.e.: 'EXP', 'IDA', 'IEA',... or any subset;
%                               leave empty to consider all evidences.
%--------------------------------------------------------------------------
% OUTPUTS:
%   value			(scalar)	value of BHI index
%
%--------------------------------------------------------------------------
% REQUIRES:
%    FastSet routines (fast_unique, fast_ismember_sorted) provided by 
%    Lev Muchnik at: http://www.levmuchnik.net/Content/ProgrammingTips/MatLab/FastSet/FastSet.html
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2012  Nejc Ilc
% Part of Pepelka package based on R package clValid.
%
%------- REFERENCE --------------------------------------------------------
% Datta, S., & Datta, S. (2006). 
% Methods for evaluating clustering algorithms for gene expression data 
% using a reference set of functional classes. 
% BMC bioinformatics, 7, 397. doi:10.1186/1471-2105-7-397
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 6-Nov-2012 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <myName.mySurname@fri.uni-lj.si>
%==========================================================================
% tic
K = max(labels); % cluster ids

if exist('ignoreEvidence','var') && ~isempty(ignoreEvidence)
    evidenceMask = ~ismember({annotations.evidence}, ignoreEvidence);
else
    evidenceMask = 1;
end

if exist('aspect','var') && ~isempty(aspect)
    aspectMask = ismember({annotations.aspect}, aspect);
else
    aspectMask = 1;
end

geneMask = ismember({annotations.geneID}, genenames);

mask = (geneMask & evidenceMask & aspectMask);
geneID = {annotations(mask).geneID};
functionID = uint32(str2double({annotations(mask).functionID}));

BHI = zeros(1,K);

% prepare data once to speedup for-loops
goList = cell(length(genenames),1);
for g=1:length(genenames)
    goList{g} = fast_unique(functionID(strcmp(genenames{g},geneID)));
end


% for each cluster c
Bj ={};

for k = 1:K
    cluster = find(labels == k);
    
    BHI_cluster = 0;
    nk = 0;
    %for each gene i in c
    for i = 1:length(cluster)-1
        % count how many other genes from c share any of gene i's
        % annotation
        Bi = goList{cluster(i)};
        
        % number of annotated genes in cluster
        nk = nk + ~isempty(Bi);
        
        for j = (i+1) : length(cluster)
            Bj = goList{cluster(j)};
            I = any(fast_ismember_sorted(Bi,Bj));
            BHI_cluster = BHI_cluster + I;
        end
        
        
    end
    nk = nk + ~isempty(Bj); % consider also the last cluster
    BHI(k) = 2* BHI_cluster / (nk*(nk-1)); % divide by the number of pairs
    
end
BHI = mean(BHI(:));
% toc
