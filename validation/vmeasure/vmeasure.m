function VM = vmeasure(arg1, arg2, beta)
% VM = vmeasure(arg1, arg2, beta)
%--------------------------------------------------------------------------
% Calculates external measure of matching between two paritions using
% V-measure (Rosenberg2007).
%--------------------------------------------------------------------------
% arg1: labels of true classes or contingency matrix
% arg2: labels of clustering or parameter beta (if arg1 is contingency matrix)
% beta: parameter beta, if arg1 and arg2 are labels vector
%--------------------------------------------------------------------------


% function can accept labels or contingency matrix
if nargin == 1
	cT = arg1;
	beta = 1.0;
	
elseif nargin == 2
	if isvector(arg1)
		cT = contingency(arg1,arg2);
		beta = 1.0;
	else
		cT = arg1;
		beta = arg2;
	end
	
else 
	cT = contingency(arg1,arg2);
end

homogeneity = getHomogeneity(cT);
completeness = getCompleteness(cT);

if (homogeneity + completeness == 0.0) 
	VM = 0;
else
	VM = (1 + beta) * homogeneity * completeness / ((beta * homogeneity) + completeness);
end
end %function vmeasure



% =========================================================================
function Cont=contingency(Mem1,Mem2)

if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments');
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
end
	

function compl = getCompleteness(cT)
% compl = 1-H(K|C)/H(K)
% cT - confusion matrix

class_v = sum(cT,1); % colSum
clust_v = sum(cT,2); % rowSum

n = sum(class_v);

h_k_c = 0;
h_k = calcEntropy(clust_v);
joint = calcJointEntropy(cT);

[numRows,numCols] = size(cT);

for i=1:numRows
  for j=1:numCols
	if (cT(i,j) ~= 0)
	  h_k_c = h_k_c - (cT(i,j) / n) * (log(cT(i,j) / class_v(j)) / log(2));
	end
  end
end

% disp(['h(k):', num2str(h_k)]);
% disp(['h(k|c):', num2str(h_k_c)]);
% disp(['joint:', num2str(joint)]);

if (h_k < 0) 
	h_k = 0;
end
if (h_k_c < 0) 
	h_k_c = 0;
end

if ((h_k_c / h_k > 1.0) || (1 - (h_k_c / h_k) < 0.000000001)) 
	compl = 0;
else 
	if h_k == 0
		compl = 1;
	else
		compl = 1 - (h_k_c / h_k);
	end
end

end % function getCompleteness



function h = calcEntropy(vec)

h = 0;
if all(vec == 0)
	ar = vec;
else
	ar = vec./sum(vec);
end

for i=1:length(ar)
  if (ar(i) ~= 0) 
	  h = h - ar(i) * (log(ar(i)) / log(2));
  end
end
end %function calcEntropy

function joint = calcJointEntropy(cT)

joint = 0.0;
n = sum(sum(cT));
[numRows,numCols] = size(cT);

for i=1:numRows
  for j=1:numCols
	if (cT(i,j) ~= 0)
	  joint = joint - (cT(i, j) / n) * (log(cT(i, j) / n) / log(2));
	end
  end
end
end % function calcJointEntropy



function hom = getHomogeneity(cT)
% homogeneity = 1 - H(C|K)/H(C)
h_c_k = 0;

class_v = sum(cT,1); % colSum
clust_v = sum(cT,2); % rowSum
n = sum(class_v);

h_c = calcEntropy(class_v);
joint = calcJointEntropy(cT);

[numRows,numCols] = size(cT);

for i=1:numRows
  for j=1:numCols
	if (cT(i, j) ~= 0)
	  h_c_k = h_c_k - (cT(i, j) / n) * (log(cT(i, j) / clust_v(i)) / log(2));
	end
  end
end

% disp(['h(c):', num2str(h_c)]);
% disp(['h(c|k):', num2str(h_c_k)]);
% disp(['joint:', num2str(joint)]);

if (h_c < 0) 
	h_c = 0;
end
if (h_c_k < 0) 
	h_c_k = 0;
end

% handle numerical error that can occur when h_c_k and h_c are nearly identical.
if ((h_c_k / h_c > 1) || (1 - (h_c_k / h_c) < 0.000000001)) 
	hom = 0.0;
else
	if (h_c == 0) 
		hom = 1;
	else
		hom = 1 - (h_c_k / h_c);
	end
end
end % function getHomogeneity

