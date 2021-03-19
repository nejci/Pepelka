function [validExt,list,moreInfo] = pplk_validExt(labelsT,labelsCons,methods,options)
% [validExt,list,moreInfo] = pplk_validExt(labelsT,labelsCons,methods,options)
%
%   External validity indices for comparing two partitions
%   (ground truth vs. arbitrary partition).
%
%
%   labelsT - ground truth clustering
%   labelsCons - clusterer's output = compared clustering
%
% methodList is a cell containing one or more method's identifier.
%   if empty, all of them are included into result validExt.
%
% Implemented methods:
%   - Rand Index (method='RI'),
%   - Adjusted Rand Index (method='ARI'),
%   - Jaccard Index (method='JI'),
%   - Fowlkes and Mallows (method='FM')
%   - Clustering Accuracy (method='CA')
%   - Balanced Clustering Accuracy (method='BCA')
%   - Posterior Distribution Balanced Clustering Accuracy (method='PDBCA')
%   - Variation of information (method='VOI')
%   - ADCO (method='ADCO', needed: options.data, options.bins)
%   - Normalized Mutual Information (NMIsqrt) (method='NMI')
%   - MAX Normalized Mutual Information (method='NMIMAX')
%   - Adjusted Mutual Information (AMImax) (method = 'AMI')
%	- V-measure (method='VM')
%	- B-Cubed measures (method='B3C' for B-Cubed F-measure averaged over Clusters;
%                       method='B3E' for B-Cubed F-measure averaged over Elements)
%
% Acknowledgements:
% [1]   Cluster Validity Analysis Platform (CVAP) (Version 3.4) 2006-2007 by Kaijun Wang.
%       (http://www.mathworks.com/matlabcentral/fileexchange/14620)
%
% [2]   SpectraLIB - Package for symmetric spectral clustering written by
%       Deepak Verma (http://www.stat.washington.edu/spectral/#code)
%
%==========================================================================

callDir=chdir(pplk_homeDir());

moreInfo = [];

%===========Arguments checking================================
if nargin < 2; error('Too few arguments!'); end

% unify labels
[uT,~,iT] = unique(labelsT);
labelsT = iT;
numClustersTrue = length(uT);

[uP,~,iP] = unique(labelsCons);
labelsCons = iP;
numClustersPredicted = length(uP);


%===========User defined functions============================
f=fopen(['..',filesep,'validation',filesep,'validExt.info'],'r');
usrfile=textscan(f,'%s %s %s %*[^\n]','delimiter','\t','commentStyle','%');
fclose(f);

usrMeth=usrfile{1};
usrAbbr=usrfile{2};
usrDesc=usrfile{3};

% remove duplicated entries
[~,ia,~]=unique(usrAbbr);
if length(ia) ~= length(usrMeth)
    fprintf(1,'WARNING: Found duplicates in user list, considerind only lines: %s.\n',sprintf('%d, ',sort(ia)));
end

usrMeth = usrMeth(ia);
usrAbbr = usrAbbr(ia);
usrDesc = usrDesc(ia);

if isempty(methods)
    disp('No method given - computing all of them, incl. user-defined!');
    methods=[{'RI','ARI','JI','FM','CA','BCA','PDBCA','VOI','ADCO','NMI',...
        'NMIMAX','AMI','VM','B3C','B3E'}, usrAbbr'];  
end

nMethods=length(methods);


%============Common values=============================
%form contingency table (or confusion matrix)
%Required by: RI, ARI, JI, FM, AMI, CA, BCA, PDBCA, VM
if any(ismember(methods,{'RI','ARI','JI','FM','AMI','CA','BCA','PDBCA','VM'}))
    C = getcm(labelsT,labelsCons);
    n = length(labelsT);
    nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
    njs=sum(sum(C,1).^2);		%sum of squares of sums of columns    
    ns = n*(n-1)/2;	         %total number of pairs  (n 2)
    sumC=sum(sum(C.^2));	 %sum over rows & columnns of nij^2
    sumij = nis+njs;    
    R = ns+sumC-sumij*0.5;    %no. agreements; disagreements -t2+sumij*0.5;
    
    % some methods work with confusion matrix that is adjusted (cluster aligned
    % to classed by majority)
    if any(ismember(methods,{'CA','BCA','PDBCA'}))
        oldPath=chdir(['..',filesep,'validation',filesep,'PDBAC']);
        [Caligned,~,hungarianCost] = getcmClust(C);
        moreInfo.ConfMat = Caligned;
        chdir(oldPath);
    end
end

% Precompute I(x;y), H(x), H(y) for information-based methods
if any(ismember(methods,{'AMI','NMI','NMIMAX','VOI'}))
    [NMI,MI,Hx,Hy] = nmi(labelsT,labelsCons);
    Hx = full(Hx);
    Hy = full(Hy);
end

validExt = struct();
list = zeros(1,nMethods);

for ind=1:nMethods
    currMethod=upper(methods{ind});
    
    switch(currMethod)
        
        case 'RI'
            %============= Rand index =====================================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % W. M. Rand, �Objective Criteria for the Evaluation of
            % Clustering Methods,� Journal of the American Statistical
            % Association, vol. 66, no. 336, pp. 846�850, 1971.
            validExt.RI = R/ns;
            list(ind) = validExt.RI;
            
        case 'ARI'
            %============= Adjusted Rand index ============================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % L. Hubert and P. Arabie, �Comparing partitions,� Journal
            % of Classification, vol. 2 , pp. 193�218, 1985.
            nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));
            if ns==nc
                ARI=0;			   %to avoid division by zero; if k=1, define Rand = 0
            else
                ARI=(R-nc)/(ns-nc);
            end
            validExt.ARI=ARI;
            list(ind) = validExt.ARI;
            
        case 'JI'
            %============= Jaccard index ==================================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % Jain, A., Dubes, R.: Algorithms for Clustering Data.
            % Prentice-Hall, Englewood Cliffs (1988)
            validExt.JI = (sumC-n)/(sumij-sumC-n);
            list(ind) = validExt.JI;
            
        case 'FM'
            %============= Fowlkes and Mallows ============================
            % 0 - labelsT is identical to labelsCons
            % 1 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % E. B. Fowlkes and C. L. Mallows. A method for comparing two
            % hierarchical clusterings. Journal of the American Statistical
            % Association, 78(383):553�569, 1983.
            ni = sum(C,2);
            ni = ni.*(ni-1)/2;
            nis = sum(ni);
            nj = sum(C,1);
            nj = nj.*(nj-1)/2;
            njs = sum(nj);
            validExt.FM = 1 - 0.5*(sumC-n)/sqrt(nis*njs);
            list(ind) = validExt.FM;
            
        case 'CA'
            %============= Clustering Accuracy ============================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % M. Meil�, �Comparing clusterings -- an information based
            % distance,� Journal of Multivariate Analysis, vol. 98, pp.
            % 873�895, 2007.
            validExt.CA = -hungarianCost/n;
            list(ind) = validExt.CA;
            
        case 'BCA'
            %============= Balanced Clustering Accuracy ===================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % H. Carrillo, K. H. Brodersen, and J. A. Castellanos,
            % �Probabilistic Performance Evaluation for Multiclass
            % Classification Using the Posterior Balanced Accuracy,� in
            % ROBOT2013: First Iberian Robotics Conference SE - 25, vol.
            % 252, M. A. Armada, A. Sanfeliu, and M. Ferre, Eds. Springer
            % International Publishing, 2014, pp. 347�361.
            a = diag(Caligned)./sum(Caligned,2);
            a(isnan(a))=0;
            % NAPAKA! 
            %validExt.BCA = (1/size(Caligned,1))*sum(a);
            % namesto size(Caligned,1) bi moralo biti stevilo
            % razredov v labelsT
            validExt.BCA = (1/numClustersTrue)*sum(a);
            list(ind) = validExt.BCA;
            
            
        case 'PDBCA'
            %=== Posterior Distribution Balanced Classification Accuracy ==
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % H. Carrillo, K. H. Brodersen, and J. A. Castellanos,
            % �Probabilistic Performance Evaluation for Multiclass
            % Classification Using the Posterior Balanced Accuracy,� in
            % ROBOT2013: First Iberian Robotics Conference SE - 25, vol.
            % 252, M. A. Armada, A. Sanfeliu, and M. Ferre, Eds. Springer
            % International Publishing, 2014, pp. 347�361.
            alpha = 0.05;
            res = 0.001;
            computeCI = 0;
            show = 0;
            
            if exist('options','var') && isstruct(options)
                if(isfield(options,'alpha'))
                    alpha = options.alpha;
                end
                if(isfield(options,'res'))
                    res = options.res;
                end
                if(isfield(options,'computeCI'))
                    computeCI = options.computeCI;
                end
                if(isfield(options,'show'))
                    show = options.show;
                end
            end
            
            oldPath=chdir(['..',filesep,'validation',filesep,'PDBAC']);
            
            if computeCI
                [bmean,CI] = PDBAC(Caligned,[],alpha,show,res);
                moreInfo.PDBAC_CI = CI;
            else
                bmean = PDBAC(Caligned,[],alpha,show,res);
            end
            
            validExt.PDBCA = bmean;
            list(ind) = validExt.PDBCA;
            
            chdir(oldPath);
            
        case 'VOI'
            %============ Variation Of Information ========================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            % original VOI index is distance metric: 1-validExt.VOI
            %--------------------------------------------------------------
            % M. Meil�, �Comparing clusterings - an information based
            % distance,� Journal of Multivariate Analysis, vol. 98, pp.
            % 873�895, 2007.
            VOI = 1 - (Hx + Hy - 2*MI)/log(length(labelsT));
            validExt.VOI = VOI;
            list(ind) = validExt.VOI;
            
        case 'ADCO'
            %===== Attribute Distribution Clustering Orthogonality ========
            % 0 - labelsT is identical to labelsCons
            % 1 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % E. Bae, J. Bailey, and G. Dong, �Clustering Similarity
            %Comparison Using Density Profiles,� in in AI 2006: Advances in
            %Artificial Intelligence, vol. 4304 , A. Sattar and B.-H. Kang,
            %Eds. Springer Berlin / Heidelberg , 2006, pp. 342�351.
            
            if (exist('options','var'))
                if(isfield(options,'data'))
                    if(isfield(options,'bins'))
                        validExt.ADCO=ADCO(options.data, labelsT, labelsCons,options.bins);
                    else
                        %default value for bins is 10
                        validExt.ADCO=ADCO(options.data, labelsT, labelsCons,10);
                    end
                    list(ind) = validExt.ADCO;
                else
                    disp('Warning: ADCO -> Missing options.data field! Ignoring ADCO!');
                end
            else
                disp('Warning: ADCO -> Missing options variable! Ignoring ADCO!');
            end
            
        case 'NMI'
            %============ Normalized Mutual Information ===================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % A. Strehl and J. Ghosh, �Cluster ensembles - a knowledge
            % reuse framework for combining multiple partitions,� The
            % Journal of Machine Learning Research, vol. 3, pp. 583�617,
            % 2003.
            validExt.NMI=NMI;
            list(ind) = validExt.NMI;
            
        case 'NMIMAX'
            %============ Normalized Mutual Information Max ===============
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % NMIMAX proposed in:
            % A. Kraskov, H. St�gbauer, R. G. Andrzejak, and P.
            % Grassberger, �Hierarchical clustering using mutual
            % information,� Europhysics Letters, vol. 70, no. 2, pp.
            % 278�284, 2005.
            NMIMAX = MI/max(Hx,Hy);
            if abs(NMIMAX)<eps
                NMIMAX = 0;
            end
            validExt.NMIMAX=NMIMAX;
            list(ind) = validExt.NMIMAX;
            
        case 'AMI'
            %============ Adjusted Mutual Information =====================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % N. X. Vinh, J. Epps, and J. Bailey, �Information Theoretic
            % Measures for Clusterings Comparison: Variants, Properties,
            % Normalization and Correction for Chance,� Journal of Machine
            % Learning Research, vol. 11, pp. 2837�2854, Dec. 2010.
            AMI = ami(C,[],MI,Hx,Hy);
            validExt.AMI = AMI;
            list(ind) = validExt.AMI;
            
        case 'VM'
            %============= V-Measure ======================================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % A. Rosenberg and J. Hirschberg, �V-Measure: A Conditional
            % Entropy-Based External Cluster Evaluation Measure,� in
            % Proceedings of the 2007 Joint Conference on Empirical Methods
            % in Natural Language Processing and Computational Natural
            % Language Learning, 2007, pp. 410�420.
            beta = 1.0;
            if (exist('options','var'))
                if(isfield(options,'beta'))
                    beta = options.beta;
                end
            end
            
            oldPath=chdir(['..',filesep,'validation',filesep,'vmeasure']);
            validExt.VM = vmeasure(C,beta);
            list(ind) = validExt.VM;
            chdir(oldPath);
            
        case 'B3C'
            %============= B-Cubed Cluster ================================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % A. Bagga and B. Baldwin, �Algorithms for Scoring Coreference
            % Chains,� in In The First International Conference on Language
            % Resources and Evaluation Workshop on Linguistics Coreference,
            % 1998, pp. 563�566.
            oldPath=chdir(['..',filesep,'validation',filesep,'BCUBED']);
            validExt.B3C = bCubedCluster(labelsT,labelsCons);
            list(ind) = validExt.B3C;
            chdir(oldPath);
            
        case 'B3E'
            %============= B-Cubed Element ================================
            % 1 - labelsT is identical to labelsCons
            % 0 - labelsT is totally different to labelsCons
            %--------------------------------------------------------------
            % A. Bagga and B. Baldwin, �Algorithms for Scoring Coreference
            % Chains,� in In The First International Conference on Language
            % Resources and Evaluation Workshop on Linguistics Coreference,
            % 1998, pp. 563�566.
            oldPath=chdir(['..',filesep,'validation',filesep,'BCUBED']);
            validExt.B3E = bCubedElement(labelsT,labelsCons);
            list(ind) = validExt.B3E;
            chdir(oldPath);
            
        otherwise
            % Check for user-defined imported functions, relating file
            % validExt.info .
            
            % search is made in the list of abbreviations, stored in
            % usrAbbr string array.
            idx = find(strcmp(currMethod,upper(usrAbbr)));
            if length(idx) > 1
                disp(['Multiple matches for ', methods{ind}, '! Considering first occurence.']);
                idx = idx(1);
            end
            
            if ~isempty(idx)
                oldPath=chdir(['..',filesep,'validation',filesep]);
                disp(['USER-DEFINED :: Executing function: ',usrAbbr{idx},'=',usrMeth{idx},'(...) -- ',usrDesc{idx}])
                validExt.(usrAbbr{idx})=feval(str2func(usrMeth{idx}),labelsT,labelsCons);
                list(ind) = validExt.(usrAbbr{idx});
                chdir(oldPath);
            else
                disp(['Non-existing method name: ', methods{ind}, '! Ignoring.']);
                validExt  = NaN;
                list = NaN;
            end
    end
end

chdir(callDir);

end

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

function confus = getcm(target,pred)

% GETCM : gets confusion matrices, precision, recall, and F scores
% [confus,numcorrect,precision,recall,F] = getcm (actual,pred,[classes])
%
% actual is a N-element vector representing the actual classes
% pred is a N-element vector representing the predicted classes
%
% dinoj@cs.uchicago.edu , Apr 2005, modified July 2005
% [confus,numcorrect,precision,recall,F] = getcm ([1 2 3 2 2 3 1 2 3],[1 2 3 3 2 1 1 2 3],1:3)

% if size(target,1) ~= size(pred,1)
%     pred=pred';
% end
% classes = 1:max(max(target),max(pred));
% 
% nC = length(classes);
% confus = zeros(nC);
% for i=1:nC
%     a = classes(i);
%     d = target==a;     % d ones where points are with class a
%     for j=1:nC
%         confus(i,j) = length(find(pred(d)==classes(j)));
%     end
% end

if size(target,2) == 1
    target = target';
end
if size(pred,2) == 1
    pred = pred';
end

classes = 1:max(max(target),max(pred));
nC = classes(end);
confus = zeros(nC);
for i=1:nC
    d = target==i;     % d ones where points are with class a
    confus(i,:) = sum(bsxfun(@eq,pred(d)',classes),1);
end

end


function simADCO=ADCO(data, c1, c2,q)
%Implemented according to the paper
%Eric Bae, James Bailey, Guozhu Dong: Clustering Similarity Comparison using density Profiles
%
%Coding by: Nejc Ilc

if ~exist('q','var')
    q=10;
end

%1. atribute vsakega vzorca (feature) razdelimo v q celic. Za vsak atribut
%najprej izra�unamo min in max ter nato ta interval razdelimo linearno na q
%podintervalov.

%N je �tevilo vzorcev, A je �tevilo atributov
[~,A]=size(data);

K=length(unique(c1));

if(K~=length(unique(c2)))
    error('Number of clusters in c1 and c2 must be the same!');
end


densC1=zeros(K,q*A);
densC2=zeros(K,q*A);

for k=1:K
    for a=1:A
        min_v=min(data(:,a))-1e-6;
        max_v=max(data(:,a))+1e-6;
        grid=linspace(min_v,max_v,q+1);
        
        for bin=1:q
            %pre�tejemo to�ke, ki padejo v vsak interval pri dani razvrstitvi (c1 in
            %c2)
            p1=data(c1==k,a);
            p2=data(c2==k,a);
            densC1(k,(a-1)*2+bin)=sum( p1 >= grid(bin) & p1 < grid(bin+1) );
            densC2(k,(a-1)*2+bin)=sum( p2 >= grid(bin) & p2 < grid(bin+1) );
        end
        
    end
    
end

%ra�unanje podobnosti kot skalarni produkt vektorjev dens - ra�unati moramo
%z vsemi permutacijami in kot rezultat vzamemo tisto z najve�jo vrednostjo

% number of permutations of K clusters
nP=factorial(K);

% all permutations of K clusters
allPerms=perms(1:K);
maxDotProd=0;

for p=1:nP
    
    dotProd=sum(sum(densC1 .* densC2(allPerms(p,:),:)));
    
    if dotProd > maxDotProd
        maxDotProd=dotProd;
    end
    
end
% maximal possible similarity of cluster
maxSim=max(sum(sum(densC1.*densC1)),sum(sum(densC2.*densC2)));
simADCO=1- maxDotProd/maxSim;
end

function [NMI,MI,Hx,Hy] = nmi(x, y)
% Compute nomalized mutual information I(x,y)/sqrt(H(x)*H(y)).
% Written by Michael Chen (sth4nth@gmail.com).
% Modified by Nejc Ilc (log2 -> log).

assert(numel(x) == numel(y));
n = numel(x);
x = reshape(x,1,n);
y = reshape(y,1,n);

l = min(min(x),min(y));
x = x-l+1;
y = y-l+1;
k = max(max(x),max(y));

idx = 1:n;
Mx = sparse(idx,x,1,n,k,n);
My = sparse(idx,y,1,n,k,n);
Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
Hxy = -dot(Pxy,log(Pxy+eps));

Px = mean(Mx,1);
Py = mean(My,1);

% entropy of Py and Px
Hx = -dot(Px,log(Px+eps));
Hy = -dot(Py,log(Py+eps));

% mutual information
MI = Hx + Hy - Hxy;

% normalized mutual information
NMI = real(sqrt((MI/Hx)*(MI/Hy))) ;
end

function [AMI] = ami(true_mem,mem,MI,Ha,Hb)
%Program for calculating the Adjusted Mutual Information (AMI_max) between
%two clusterings, tested on Matlab 7.0 (R14)
%(C) Nguyen Xuan Vinh 2008-2010
%Contact: n.x.vinh@unsw.edu.au
%         vthesniper@yahoo.com
% Changed by Nejc Ilc, 2014
%--------------------------------------------------------------------------
%**Input: a contingency table T
%   OR
%        cluster label of the two clusterings in two vectors
%        eg: true_mem=[1 2 4 1 3 5]
%                 mem=[2 1 3 1 4 5]
%        Cluster labels are coded using positive integer.
%**Output: AMI: adjusted mutual information  (AMI_max)
%
%**Note: In a prevous published version, if you observed strange AMI
%results, eg. AMI>>1, then it's likely that in these cases the expected MI
%was incorrectly calculated (the EMI is the sum of many tiny elements, each
%falling out the precision range of the computer). However, you'll likely
%see that in those cases, the upper bound for the EMI will be very tiny,
%and hence the AMI -> NMI (see [3]). It is recommended setting AMI=NMI in
%these cases, which is implemented in this version.
%--------------------------------------------------------------------------
%References:
% [1] 'A Novel Approach for Automatic Number of Clusters Detection based on
%       Consensus Clustering',
%       N.X. Vinh, and Epps, J., in Procs. IEEE Int. Conf. on
%       Bioinformatics and Bioengineering (Taipei, Taiwan), 2009.
% [2] 'Information Theoretic Measures for Clusterings Comparison: Is a
%	    Correction for Chance Necessary?', N.X. Vinh, Epps, J. and Bailey, J.,
%	    in Procs. the 26th International Conference on Machine Learning (ICML'09)
% [3] 'Information Theoretic Measures for Clusterings Comparison: Variants,
%       Properties, Normalization and Correction for Chance', N.X. Vinh,
%       Epps, J. and Bailey, J., Journal of Machine Learning Research,
%       11(Oct), pages 2837-2854, 2010

if ~exist('mem','var') || isempty(mem)
    T=true_mem; %contingency table pre-supplied
    % remove zero cols and rows
    mask = T == 0;
    colInd = all(mask,1);
    rowInd = all(mask,2);
    T = T(~rowInd,~colInd);
else
    %build the contingency table from membership arrays
    R=max(true_mem);
    C=max(mem);
    %n=length(mem);N=n;
    
    %identify & removing the missing labels
    list_t=ismember(1:R,true_mem);
    list_m=ismember(1:C,mem);
    T=Contingency(true_mem,mem);
    T=T(list_t,list_m);
end
n=sum(sum(T));
N=n;

%update the true dimensions
[R C]=size(T);
if C>1
    a=sum(T,2)';
else
    a=T';
end
if R>1
    b=sum(T);
else
    b=T;
end

% %calculating the Entropies
% Ha=-(a/n)*log(a/n)';
% Hb=-(b/n)*log(b/n)';

% %calculate the MI (unadjusted)
% MI=0;
% for i=1:R
%     for j=1:C
%         if T(i,j)>0
%             MI=MI+T(i,j)*log(T(i,j)*n/(a(i)*b(j)));
%         end
%     end
% end
% MI=MI/n;

%-------------correcting for agreement by chance---------------------------
AB=a'*b;
bound=zeros(R,C);

E3=(AB/n^2).*log(AB/n^2);

EPLNP=zeros(R,C);
LogNij=log((1:min(max(a),max(b)))/N);
for i=1:R
    for j=1:C
        nij=max(1,a(i)+b(j)-N);
        X=sort([nij N-a(i)-b(j)+nij]);
        if N-b(j)>X(2)
            nom=[(a(i)-nij+1:a(i)) (b(j)-nij+1:b(j)) (X(2)+1:N-b(j))];
            dem=[(N-a(i)+1:N) (1:X(1))];
        else
            nom=[(a(i)-nij+1:a(i)) (b(j)-nij+1:b(j))];
            dem=[(N-a(i)+1:N) (N-b(j)+1:X(2)) (1:X(1))];
        end
        p0=prod(nom./dem)/N;
        
        sumPnij=p0;
        
        EPLNP(i,j)=nij*LogNij(nij)*p0;
        p1=p0*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);
        
        for nij=max(1,a(i)+b(j)-N)+1:1:min(a(i), b(j))
            sumPnij=sumPnij+p1;
            EPLNP(i,j)=EPLNP(i,j)+nij*LogNij(nij)*p1;
            p1=p1*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);
            
        end
        CC=N*(a(i)-1)*(b(j)-1)/a(i)/b(j)/(N-1)+N/a(i)/b(j);
        bound(i,j)=a(i)*b(j)/N^2*log(CC);
    end
end

EMI_bound=sum(sum(bound));
EMI=sum(sum(EPLNP-E3));

AMI=(MI-EMI)/(max(Ha,Hb)-EMI);
NMI=real(MI/sqrt(Ha*Hb));

if isnan(NMI)
    NMI = 0;
end

%If expected mutual information negligible, use NMI.
if abs(EMI)>EMI_bound
    %fprintf('The EMI is small: EMI < %f, setting AMI=NMI',EMI_bound);
    AMI=NMI;
end
end

function Cont=Contingency(Mem1,Mem2)

if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
    error('Contingency: Requires two vector arguments');
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
    Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
end

