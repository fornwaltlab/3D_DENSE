% qtClustering - mex-file Implementation of the QT clustering algorithm
%
%   An improvement upon qtclusteuclid by Misha Koshelev
%
%   Performs QT clustering as desribed in:
%
%   Heyer, L. J., Kruglyak, S., Yooseph, S. (1999). Exploring expression
%   data: Identification and analysis of coexpressed genes. Genome Research
%   9, 1106–1115.
%
%   http://genome.cshlp.org/content/9/11/1106.full
%   http://genome.cshlp.org/content/9/11/1106/F5.large.jpg
%
%   For this particular implementation, euclidean distance is used if type == 0
%   and RMSE if type == 1.
%
% USAGE:
%   Ndx = qtClustering(G,dia,type)
%
% INPUTS:
%   G:      [M x N] Array, Data to cluster where each row is a separate data
%           point i.e. [X1 Y1 Z1; X2 Y2 Z2]
%   dia:    Scalar, Diameter cutoff for QT clustering
%   type:   Scalar, Indicates which parameter we want to use. They are:
%           0 => Euclidean distance [Default]
%           1 => Root Mean Squared Error
%
% OUTPUTS:
%   Ndx:    [M x 1] Vector, The index of the cluster to which each data point
%           has been assigned via the clustering.
%
% See also:
%   qtclusteuclid, qtClustLoc
%
% Last Modified: 03-05-2012
% Modified By: Suever (suever@gmail.com)
