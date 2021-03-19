+-------------------------------------------------------------------------+
|             CLUSTER VALIDATION WITH MODIFIED DUNN'S INDEX               |
+-------------------------------------------------------------------------+

Copyright (C) 2019, Nejc Ilc, version 1.0
This software is free and distributed under BSD license.


------- BRIEF -------------------------------------------------------------
This package computes the modified Dunn's cluster validity index 
using shortest paths in the neighbourhood graph (Ilc, 2012). It also includes 
Generalized Dunn's index (Pal & Biswas, 1997) and a bunch of algorithms for
graph building.

Best cluster partition maximizes the index value.


------- HOW TO USE -------------------------------------------------------
1)	Run script file 'demo' to see illustrative example.

2)	Run function 'indexDNs' (or 'indexDNg') with specified parameters.
    This function builds a graph on input data. If the graph has been 
    already built, you should run 'indexDNs_graph' (or 'indexDNg_graph')
    instead to bypass (slow) computing of the graph.

3)  Type 'help indexDNs' or 'help indexDNg' to get a list of the parameters 
    and their description.
 

------- FILES and FOLDERS -------------------------------------------------
  indexDNs.m		: main function for Modified Dunn's index
  indexDNg.m		: main function for Generalized Dunn's index
  indexDNs_graph.m	: compute Modified Dunn's index on an existing graph
  indexDNg_graph.m	: compute Generalized Dunn's index on an existing graph
  
  demo.m		    : script file with illustrative example 
                      (requires Statistic toolbox for single-linkage alg.)
  
  graph_create.m    : wrapper function that calls other functions for graph
                      construction
  graph_(method).m  : constructs graph on data points according to one of 
                      the methods:
                        -> Directed KNN graph (1)
                        -> Euclidean Minimum Spanning Tree graph
                        -> Epsilon-neighborhood undirected graph (1)
                        -> Gabriel graph (2)
                        -> Mutual KNN graph (1)
                        -> Relative Neighborhood Graph (2)
                        -> Symmetric KNN graph (1)
  putbnd, putbnds.m : helper functions written by Richard E. Strauss (2)

Credits for graph-constructing methods (labeled with (1) and (2)) go to:

(1) Matthias Hein and Ulrike von Luxburg (2007). Code was used at the 
practical sessions of the Machine Learning Summer School 2007 at the Max 
Planck Institute for Biological Cybernetics in Tuebingen, Germany 
http://www.mlss.cc/tuebingen07/

(2) Strauss, R. E. (2012). Matlab statistical functions [computer software]. 
Retrieved August 28, 2012, from 
http://www.faculty.biol.ttu.edu/Strauss/Matlab/matlab.htm.

Authors are appreciated for their valuable contribution to this package.  

------- CITING OUR WORK ---------------------------------------------------
If you find modified Dunn's index useful, please cite the paper:

    Ilc, N. (2012). Modified Dunn�s cluster validity index based on graph
    theory. Przeglad Elektrotechniczny (Electrical Review), (2), 126-131.

------- CONTACT ----------------------------------------------------------
Bug reports, comments and questions are appreciated.
Please write to: 
	Nejc Ilc <nejc.ilc@fri.uni-lj.si>

------- LEGAL NOTICE ------------------------------------------------------
Copyright (c) 2019, Nejc Ilc
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

LICENSING THE LIBRARY FOR PROPRIETARY PROGRAMS
It is not possible to include this software library in a commercial
proprietary program without written permission from the owners of the
copyright. If you wish to obtain such permission, you can reach us by
paper mail:
    Nejc Ilc
    Faculty of computer and information science
    Vecna pot 113
    1000 Ljubljana
    Slovenia

or by email:
    nejc.ilc@fri.uni-lj.si
---------------------------------------------------------------------------
