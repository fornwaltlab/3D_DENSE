classdef RBFInterpolator < hgsetget
% RBFInterpolator - Interpolator that uses radial basis functions
%
%   The RBFInterpolator class allows the user to easily compute radial
%   basis functions for the supplied data. This class is extrenemely
%   flexible and can be used for N-D data and can simultaneously
%   interpolate several data associated with the specified points
%
%   You can create an interpolator from a series of points and values:
%
%   DEMO:
%       % Create some x y points and their corresponding z values
%       [x, y] = meshgrid(-10:10, -10:10);
%       z = x.^2 + y.^2;
%
%       % Construct the interpolator
%       RBF = RBFInterpolator([x(:), y(:)], z(:));
%
%       % Interpolate the value at unknown points
%       [xq, yq] = meshgrid(linspace(-10, 10), linspace(-10, 10));
%
%       value = RBF(xq, yq);
%
%       % Now plot the results to show that the interpolation is decent
%       plot3(xq(:), yq(:), value(:), '.')
%
% USAGE:
%   R = RBFInterpolator(points, values, options)
%
% INPUTS:
%   points:     [M x N] Array, An N-dimensional set of points at which the
%               values are defined.
%
%   values:     [M x 1] Array, Values corresponding to each row in the
%               POINTS input matrix.
%
%   options:    Parameter/Value Pairs, Other properties of the
%               RBFInterpolator object that should be specified as
%               parameter and value pairs (optional)
%
% OUTPUTS:
%   R:      Object, Handle to the RBFInterpolator object which can be used
%           to query the interpolator at a specific point.
%
% Last Modified: 10-14-2013
% Modified By: Jonathan Suever (suever@gmail.com)

    properties (SetAccess = 'protected')
        Constant            % Constant to be used for some basis functions
        Points      = []    % Locations (coordinates) of corresponding values
        Data        = []    % Values associated with each location
        Normalized  = false % Is this a normalized RBF system

        % Type of radial basis function to use (function handle)
        Type        = RBFInterpolator.LINEAR
    end

    properties (Access = 'protected')
        Weights    = []    % Computed weights for basis functions
    end

    properties (Dependent)
        NPoints             % Number of data points
    end

    properties (Constant)
        GAUSSIAN        = @(r,const)RBFInterpolator.gaussian(r,const);
        THINPLATE       = @(r,const)RBFInterpolator.thinplate(r);
        CUBIC           = @(r,const)RBFInterpolator.cubic(r);
        MULTIQUADRICS   = @(r,const)RBFInterpolator.multiquadrics(r, const);
        LINEAR          = @(r,const)RBFInterpolator.linear(r);
    end

    methods
        function self = RBFInterpolator(pts, vals, varargin)
            % RBFInterpolator - RBFInterpolator constructor
            %
            %   This function creates the radial basis function which can
            %   then be queried at any point to yield a value.
            %
            %   The resulting RBFInterpolator can be indexed at any point
            %   to yield the interpolated value at that point.
            %
            %   Example:
            %
            %       R = RBFInterpolator([1 1; 2 2], [0; 1]);
            %       R(1.5, 1.5)
            %
            %       ans =
            %           0.5
            %
            % USAGE:
            %   R = RBFInterpolator(pts, vals, props/val)
            %
            % INPUTS:
            %   pts:    [M x N] Array, Coordinates associated with each
            %           data point where they are N-D coordinates.
            %   vals:   [M x S] Array, Values associated with the PTS.
            %           Multiple types of data can be associated with each
            %           point using the columns of VALS
            %
            %   Props:
            %
            %       Constant:   Constant value to use in the calculation of
            %                   certain basis functions. For example, for
            %                   the GAUSSIAN basis function, Constant is
            %                   the standard deviation of the Gaussian. The
            %                   constant is supplied as the second input to
            %                   the TYPE function handle
            %
            %       Type:       Determines the function to use as the basis
            %                   function. This function must be a function
            %                   of the radius (distance) from the origin
            %                   and optionally a constant. The following
            %                   are already provided:
            %
            %                       RBFInterpolator.GAUSSIAN
            %                       RBFInterpolator.THINPLATE
            %                       RBFInterpolator.CUBIC
            %                       RBFInterpolator.MULTIQUADRICS
            %                       RBFInterpolator.LINEAR
            %
            %                   Alternately, the user may specify a custom
            %                   basis function using an anonymous function
            %                   of the form:
            %
            %                       @(radius, constant)function()
            %
            %     Normalized:   Determines whether to use normalized RBFs
            %                   (true) or not (false, the default).
            %
            % OUTPUTS:
            %   R:      RBFInterpolator, object handle to be used to
            %           perform interpolation.

            self.Points = pts;
            self.Data = vals;

            % Make sure that the data format is correct
            assert(size(self.Points, 1) == size(self.Data, 1),...
                'Data must be the same size as your points');

            % Process all input arguments and assign the necessary props
            if numel(varargin)
                set(self, varargin{:});
            end

            % Compute the average distance between nodes
            if isempty(self.Constant)
                npts = size(pts, 1);
                tmp = bsxfun(@minus, permute(pts, [3 2 1]), pts);
                distances = reshape(sqrt(sum(tmp.^2, 2)), npts, npts);
                self.Constant = mean(distances(:));
            end

            % Go ahead and solve according to our data
            self.solve();
        end

        function solve(self)
            % solve - Determine the coefficients of the basis functions

            % Compute the radius of all points from all other points
            R = dist(self.Points');

            % Compute the value of the basis function using these radii
            A = self.Type(R, self.Constant);

            % Normalize if required
            if self.Normalized
                A = bsxfun(@rdivide, A, sum(A,2));
            end

            dims = size(self.Points, 2);

            % Polynomial part
            P = [ones(self.NPoints, 1) self.Points];
            A = [A  P; P' zeros(dims + 1, dims + 1)];
            B = [self.Data; zeros(dims + 1, size(self.Data, 2))];

            self.Weights = A \ B;
        end

        function res = query(self, varargin)
            % query - Query the interpolant at a particular point
            %
            % USAGE:
            %   r = query(R, pts)
            %   r = query(R, x, y)
            %
            % INPUTS:
            %   R:      RBFInterpolator, Object handle
            %   pts:    [M x N] Array, N-D coordinates at which we want the
            %           value of the interpolant
            %   x:      [M x 1] Array, X Coordinates of the points to query
            %   y:      [M x 1] Array, Y Coordinates of the points to query
            %
            % OUTPUTS:
            %   r:      [M x S] Array, Data values at the interpolated points

            if numel(varargin) == 1
                queries = varargin{1};
            else
                % Make sure that we flatten all arrays
                varargin = cellfun(@(x)x(:), varargin, 'UniformOutput', 0);
                queries = cat(2, varargin{:});
            end

            % Make sure that the dimensionality is correct
            assert(size(queries, 2) == size(self.Points, 2), ...
                ['The dimensionality of your query points does not ', ...
                 'match the dimensionality of your data.']);

            nNodes = self.NPoints;

            % Compute the distance from every point to every other point
            tmp = bsxfun(@minus, permute(self.Points, [3 2 1]), queries);
            distances = reshape(sqrt(sum(tmp.^2, 2)), size(tmp, 1), nNodes);

            % Compute the value of the function for all of these distances
            val = self.Type(distances, self.Constant);

            % Normalize the RBFs if necessary
            if self.Normalized
                val = bsxfun(@rdivide, val, sum(val, 2));
            end

            bsx = val * self.Weights(1:nNodes, :);
            S = bsxfun(@plus, bsx, self.Weights(nNodes + 1,:));
            res = (queries * self.Weights(nNodes + 2 : end,:)) + S;
        end

        function B = subsref(self, S)
            % subsref - Overloaded subscript referencing
            %
            %   This is overloaded because we want to be able to use
            %   indices to interpolate the RBF

            switch S(1).type
                case '()'
                    B = query(self, S(1).subs{:});
                otherwise
                    B = builtin('subsref', self, S);
            end
        end
    end

    %-- Get / Set Methods --%
    methods
        function res = get.NPoints(self)
            res = size(self.Points, 1);
        end
    end

    methods (Hidden, Static)
        %-- Built-in Basis Functions --%
        function u = linear(r, ~)
            u = r;
        end

        function u = cubic(r, ~)
            u = r .^ 3;
        end

        function u = multiquadrics(r, const)
            u = sqrt(1 + r.^r / const^2);
        end

        function u = thinplate(r, ~)
            u = r.^2 .* log(r + 1);
        end

        function u = gaussian(r, const)
            u = exp((-0.5 / (const^2)) * r.^2);
        end
    end
end
