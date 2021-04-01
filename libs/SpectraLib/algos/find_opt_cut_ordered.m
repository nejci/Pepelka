function[min_cut_point, min_cut_val]= find_opt_cut_ordered(S,optimization_option)
% function[min_cut_point, min_cut_val]= find_opt_cut_ordered(S)
%
% Helper function for cluster_kvv
% Finds the optimal cut in S assuming that S is "ordered" i.e. the cut
% can only be of the form 1:k, k+1:n 
  
  n=length(S); 
  min_cut_val=inf; 
  min_cut_point=-1; 
  
  d=sum(S,2);

  for cuti=1:n-1
	num = sum(sum(S(1:cuti,(cuti+1):n)));
	i1=1:cuti;
	i2=(cuti+1):n;
	switch (optimization_option) 
	  case {'ncut'}
	   mult_factor=(1/sum( d( i1 ))+1/sum( d( i2)));
	  case {'conductance'}
	   mult_factor=(1/min( sum(d(i1)), sum(d(i2)) ));
	 otherwise
	  warning('******** WRONG OPTIMIZATION_OPTION IN CLUSTER_SHI_R ********'); 
	  mult_factor=1; 
	end
	
	
%	den = min( sum(d(1:cuti)) , sum(d((cuti+1):n)) );
	cut_val=num*mult_factor; 
	if cut_val<min_cut_val
	  min_cut_val=cut_val; 
	  min_cut_point=cuti; 
	end
  end
  
  
