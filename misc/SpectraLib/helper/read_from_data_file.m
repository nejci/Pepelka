function [similarity,cluster_assignments,points]=read_from_data_file(filePrefix,directory)
% function [similarity,cluster_assignments,points]=read_from_data_file(filePrefix,directory)
% Reads from the data file that has been written using
% 'write_from_data_file' 
  if nargin==1
	directory='data'; 
  end 
  filePrefix=strcat(directory,'/',filePrefix); 
  
% function [similarity,cluster_assignments,points]=read_from_data_file(filePrefix)

  clear('similarity____local', 'cluster_assignment____local','points____local');

  dataFile=strcat(filePrefix,'.mat'); 
  load(dataFile  , 'similarity____local', 'cluster_assignment____local','points____local');
  similarity=similarity____local; 
  cluster_assignments=cluster_assignment____local;
  points=points____local;

  %sanity check 
  clear('similarity____local', 'cluster_assignment____local','points____local');
  
  
