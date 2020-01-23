% PolygonClip - This is a polygon clipper that takes advantage of gpc
%
%   The Polygon Clipper (based on the gpc-library) is used to perform
%   algebraic operations on two polygons.
%
%   Given two arbitrary polygons (which may self intersect, may contain holes, 
%   may be constructed of several contours) the Polygon Clipper is used to 
%   calculate the resulting polygon for the operations diff, union, AND, XOR.
%
%   All polygons are specified as structures:
%
%       S.x:    X Coordinates of the vertices
%       S.y:    y Coordinates of the vertices
%       P.hole: Boolean to determine whether it is a hold definition
%
%   These structures can be multi-element to indicate multiple polygons
%
%   NOTE:   If this function doesn't run, try recompiling it via the source in
%           the mex/polygonclipper folder (compile.m)
%
% USAGE:
%   P3 = PolygonClip(P1,P2,type)
%
% INPUTS:
%   P1,P2:  Structure, Structures describing polygons as shown above
%   type:   Integer, Must be an integer which indicates which of the following
%           to perform:
%               
%               0:  P1 - P2
%               1:  P1 and P2
%               2:  XOR(P1,P2)
%               3:  union(P1,P2)
%
% OUTPUTS:
%   P3:     Structure, Describes the polygon as detailed above
%
% Source: 
%   http://www.mathworks.com/matlabcentral/fileexchange/8818-polygon-clipper
%
% %FOOTER%

% % General
% -------
%    This folder contains the files needed to implement the Polygon Clipper "gpc" in Matlab.

%    Credit for the gpc-library goes to ... 

% 	Alan Murta
% 	Advanced Interfaces Group
% 	Department of Computer Science
% 	University of Manchester
% 	Manchester M13 9PL, UK
% 	http://www.cs.man.ac.uk/~toby/alan/software//
% 	
% 	The gateway to Matlab has been designed to work for the library version gpc 2.32 (included in this folder).
% 	If compiled with past resp. future versions of gpc, the gateway might also work, however I can not assure this and
% 	I will not give any support for versions other than 2.32 ...
% 	
% Description
% -----------
% 	The Polygon Clipper is used to calculate ...
% 		... difference
% 		... union
% 		... XOR
% 		... AND
% 	between two arbitrary polygons. Each polygon may be constructed of several contours and may contain holes (s.examples).
% 	
% Content
% -------
% 	Files in this folder are ...
% 	gpc.*:               Original Polygon Clipper files taken from the above home page (used Version: 2.32)
% 	gpc_mexfile.*:       File containing the mexfile, i.e. the gateway between Matlab and the gpc library.
% 	PolygonClip.dll:     Comiled library (mex-file for windows), which can be called from Matlab (s. Example.m).
% 	                     If this file does not work for you (e.g. under Linux), you will have to recompile.
% 	Example.m:           Simple example file, which demonstrates how to use the library in Matlab.
% 	ReadMe.txt:          This file
% Examples
% --------
% 	Try out the file Example.m for a demonstration in Matlab.
% 	Alan Murta's home page offers general examples (e.g. for Windows: ftp://ftp.cs.man.ac.uk/pub/toby/gpc/gpc_test.zip).
% How to compile
% -------------
% 	If you need to recompile the code, e.g. for Linux, try the following command, which worked for me:
% 	
% 		mex gpc.c gpc_mexfile.c -O -output PolygonClip              % optimized
% 		mex gpc.c gpc_mexfile.c -argcheck -output PolygonClip       % with argument checking
% 		mex gpc.c gpc_mexfile.c -g -output PolygonClip              % for debugging
% 		
% 	You need a C-compiler installed on your computer for this to work.
% 	
% Bugs
% ----
%     ??? 1) It seems that the routine is sometimes passing back empty arrays (0 vertex). This should be handled. ???
% PS ...
% ------
% 	I'm not exactly an experienced C-programmer. If you have any suggestions for improving the gateway, i.e. the files 
% 	gpc_mexfile.* (and only these files!!!), I will incorporate necessary changes.
% 	
% Author (for gateway)
% --------------------
% 	Dipl. Geophys. Sebastian Hölz
% 	Geophysik
% 	TU Berlin
% 	Germany
% 	hoelz_*_geophysik.tu-berlin.de   
% 	(replace _*_ with @)
