function write_data_file2(filePrefix, similarity,cluster_assignments,points,param_string)
% function write_data_file(filePrefix, similarity,cluster_assignments,points,param_string)
% Writes to datafile the specified inputs.   
  dataFile=strcat(filePrefix,'.mat'); 
  paramFile=strcat(filePrefix,'.param');
  
  % write the parameter file if the param_string is not empty. 
  if (~isempty(param_string))
	fid=fopen(paramFile,'w');
	fprintf(fid, param_string);
	fclose(fid);
  end
  
  similarity____local=similarity; 
  cluster_assignment____local=cluster_assignments;
  points____local =points;
  save(dataFile, 'similarity____local', 'cluster_assignment____local','points____local');
  clear('similarity____local', 'cluster_assignment____local','points____local');
  
  
	
  
