function uu=njw_mapping(S,k)
%function uu=njw_mapping(S,k,use_gen)
%Returns the vectors using the NJW algorithm

  %  Compute Laplacian L=D^-1/2 S D^-1/2
	D=diag(sum(S)); 
	Dsqrt=sqrt(D); 
	L=Dsqrt\S/Dsqrt;

	% the top k EV 
	[uu, dummy]=myeigs(S,k); 
	uu=uu'; 
	uu=normalize_2nd(uu); 
