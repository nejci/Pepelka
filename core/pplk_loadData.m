function [data,labelsT,info] = pplk_loadData(datasetName)
% [data,labelsT,info] = pplk_loadData(datasetName)
% Function loads data matrix from the file dataName.mat and stores its
% content into variable named 'data'. If file dataNameT.mat exists, its
% content are loaded into 'labelsT' variable, otherwise the latter becomes
% an empty vector [].
%
% INPUTS       
%   datasetName
%       - Filename in datasets folder.
%       - Name of dataset (file is selected by matching in reference table).
%       - ID number of dataset (backward compatibility).
%
%
% OUTPUTS      
%   data
%       A nPoints-by-nDimensions data matrix.
%
%   labelsT      
%       A nPoints-by-1 true labels vector stored in 
%       '..\datasets\dataNameT.mat' 
%
%
% EXAMPLES
%   [data,target,info] = pplk_loadData('real\UCI\iris');
%   [data,target,info] = pplk_loadData('iris');
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka

callDir=chdir(pplk_homeDir());


% FILENAME construction
%-------------------------------------------------------------------------
% Legacy - addressing dataset with numeric ID (GSOM pack)
% If datasetName is a number, it means ID of dataset. Added to ensure
% compatibility with older beta version. Filename of such dataset consists
% of 'dataG' and ID, eg. 'dataG2.mat'
if length(datasetName) <= 4 || strcmp(datasetName(end-3:end),'.mat')
    suffix = '';
else
    suffix = '.mat';
end

if isnumeric(datasetName)
    fileName=['..',filesep,'datasets',filesep,'dataG',num2str(datasetName),suffix];
    fileNameT=['..',filesep,'datasets',filesep,'dataG',num2str(datasetName),'T',suffix];
    
else
    % lookup for matches in dataset names table
    data_path = data_name2path(datasetName);
    if isempty(data_path)
        fileName=['..',filesep,'datasets',filesep,datasetName,suffix];
        fileNameT=['..',filesep,'datasets',filesep,datasetName,'-T',suffix];
    else
        fileName = ['..',filesep,'datasets',filesep,data_path];
        fileNameT=['..',filesep,'datasets',filesep,data_path(1:end-4),'-T',suffix];
    end
end

% LOAD data and rename variable to 'data', size [nPoints x nDimensions]
%-------------------------------------------------------------------------
tmpFields = [];
info = [];


if exist(fileName,'file')
    tmp=load(fileName);
    %Name of the variable that contains the data has to be 'data'. Next lines
    %ensure this.
    if isstruct(tmp)
        tmpFields=fieldnames(tmp);
        if any(strcmp(tmpFields,'data'))
            data = tmp.data;
        else
            data=tmp.(tmpFields{1});
        end
        if any(strcmp(tmpFields,'info'))
            info = tmp.info;
        end
    else
        data=tmp;
    end
else
    error(['File: ', fileName, ' does not exist!']);
end

% LOAD target vector (if available) and rename variable to 'labelsT', size
% [nPoints x 1]
%-------------------------------------------------------------------------
if exist(fileNameT,'file')
    tmpT=load(fileNameT);
    %Name of the variable that contains the target vector has to be 'labelsT'. Next lines
    %ensure this.
    if isstruct(tmpT)
        tmpFields=fieldnames(tmpT);
        labelsT=tmpT.(tmpFields{1});
    else
        labelsT=tmpT;
    end
else
    if any(strcmp(tmpFields,'target'))
        labelsT = tmp.target;
    else
        warning('pplk:noFile','There is no target vector.');
        labelsT=[];
    end
end
% force labelsT to be column vector
[r,c] = size(labelsT);
if ~isvector(labelsT)
   error('labelsT must be vector!'); 
end
if r==1
    labelsT = labelsT';
end

chdir(callDir);

end

% Helper functions for matching datasets names with their paths.
function data_path = data_name2path(data_name)

data_nameTable = data_getNameTable(['..',filesep,'datasets',filesep,'datasets_names.txt']);
mask = strcmpi(data_name,data_nameTable(:,2));
data_path = [];
if sum(mask)==1
    data_path = data_nameTable{mask,1};
end
end

function data_name = data_path2name(data_path)

data_nameTable = data_getNameTable(['..',filesep,'datasets',filesep,'datasets_names.txt']);
mask = strcmpi(data_path,data_nameTable(:,1));
data_name = [];
if sum(mask)==1
    data_name = data_nameTable{mask,2};
end
end

function data_nameTable = data_getNameTable(fileName)
fid = fopen(fileName,'r');
data_nameTable = textscan(fid,'%s %s','Delimiter','\t','EndOfLine','\r\n');
data_nameTable = [data_nameTable{1},data_nameTable{2}];
fclose(fid);
end