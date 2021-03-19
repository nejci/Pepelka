function [fig] = pplk_validPlot(validResults, numClust, validMethods)

% validResults  (struct)    structure with index name as field name and
%                           index values as its content (numeric vector)
%
%               (matrix)    [numValidMethods X length(numClust)] matrix
%                           with values of indeces in rows
%
% numClust      (vector)    vector of integers; for each index value 
%                           (columns of validResults) it gives the number
%                           of clusters
%
% validMethods  (cell)      if validResults is a matrix, user must specify
%                           used validation methods in this argument 
%                           as cell of strings.
%

numClust_len=length(numClust);

if isnumeric(validResults)
    if ~exist('validMethods','var') || isempty(validMethods)
        [numValidMethods, numClust_len2]=size(validResults);
        
        validMethods = cellstr([repmat('index ',numValidMethods,1),  num2str([1:numValidMethods]')]);
    end
else
    validMethods = fieldnames(validResults);
    numValidMethods = length(validMethods);
    numClust_len2 = length(validResults.(validMethods{1}));
    validResults = reshape(struct2array(validResults), numClust_len2, numValidMethods)';
end

assert(numClust_len == numClust_len2);

% max 4 subplots in a row of panel
nCols = min(numValidMethods,4);

nRows=floor(numValidMethods/nCols);
if rem(numValidMethods,nCols)~=0
	nRows=nRows+1;
end

fig=figure();
% set(0,'DefaultAxesLineStyleOrder','-|--|:|-.');
% set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1]);

   
for methodInd = 1:numValidMethods

    subplot(nRows,nCols,methodInd);
    plot(numClust, validResults(methodInd,:));
    
%     ht=title(['Internal validation']);
%     set(ht,'Interpreter','none');
%     hl=legend(clustMethods,'Location','NorthEast');	
%     set(hl,'Interpreter','none');
   

    xlim([min(numClust),max(numClust)]);
    xlabel('k');
    ylabel(validMethods{methodInd});
    
end