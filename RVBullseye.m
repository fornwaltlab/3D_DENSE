classdef RVBullseye < plugins.dense3D_plugin.Bullseye

    properties (SetObservable)
        MajorAxis   = 2
    end

    properties (Dependent)
        MinorAxis
        Septum
    end

    methods
        function res = get.Septum(self)
            if ismember(self.Segments, [12, 13])
                res = 'off';
            else
                res = 'on';
            end
        end

        function res = get.MinorAxis(self)
            res = self.Radius;
        end

        function set.MinorAxis(self, value)
            self.Radius = value;
        end

        function resetRVCache(self)
            self.rvcache = struct();
        end
    end

    methods (Access = 'protected')

        function resetCache(self, propname)
            old = self.cache;

            resetCache@plugins.dense3D_plugin.Bullseye(self, propname)

            fields = {'ahalines', 'surfcoords', 'centers'};

            if strcmpi(propname, 'CData')
                tocheck = isfield(old, fields);
                fields = fields(tocheck);

                for k = 1:numel(fields)
                    self.cache.(fields{k}) = old.(fields{k});
                end
            end
        end

        function bool = hasApex(self)
            bool = ismember(self.Segments, [13, 19]);
        end

        function [X, Y] = getSpokes(self, nPoints)

            if ~exist('nPoints', 'var')
                nPoints = 100;
            end

            switch self.Segments
                case {12, 13}
                    nSegments = 5;
                case {18, 19}
                    nSegments = 6;
                otherwise
                    error(sprintf('%s:InvalidSegments', mfilename), ...
                        'Invaild number of segments');
            end

            [tx, ty] = self.getSampledPoints(self.MinorAxis, ...
                self.MajorAxis, nPoints, nSegments);

            tx(end+1,:) = NaN;
            ty(end+1,:) = NaN;

            X = tx(:);
            Y = ty(:);

            [tx, ty] = self.getSampledPoints(self.MinorAxis, ...
                self.MajorAxis, nSegments - 1, 100, 0);

            tx = tx.';
            ty = ty.';

            tx(end+1,:) = NaN;
            ty(end+1,:) = NaN;

            X = cat(1, X, tx(:));
            Y = cat(1, Y, ty(:));

            if strcmpi(self.Apex, 'on')
                % Get the outers
                [xx,yy] = self.getSampledPoints(self.MinorAxis, self.MajorAxis, nPoints, 2, 1);

                xx(end+1,:) = NaN;
                yy(end+1,:) = NaN;

                X = cat(1, X, xx(:));
                Y = cat(1, Y, yy(:));
            end
        end

        function [X, Y, Z] = ahalines(self)

            if ~isfield(self.cache, 'ahalines')
                [X, Y] = self.getSpokes(10);
                self.cache.ahalines = [X(:), Y(:)];
            else
                tmp = self.cache.ahalines;
                X = tmp(:,1); Y = tmp(:,2);
            end

            Z = zeros(size(X)) + self.ZData + 0.01;
        end

        function centers = segmentCenters(self)

            if ~isfield(self.cache, 'centers')
                [tx, ty] = self.getSampledPoints(self.MinorAxis, ...
                    self.MajorAxis, 7, 9);

                tx(1:2:end,:) = [];
                ty(1:2:end,:) = [];
                tx(:,1:2:end) = [];
                ty(:,1:2:end) = [];

                tx = tx.';
                ty = ty.';

                centers = [tx(:), ty(:)];

                % Fill in some NaN values
                pad = nan(max(self.validSegments) - size(centers, 1), 2);

                centers = cat(1, centers, pad);

                self.cache.centers = centers;
            else
                centers = self.cache.centers;
            end
        end

        function res = validSegments(~)
            res = [13, 12];
        end
    end

    methods
        function self = RVBullseye(varargin)
            % Change the default number of segments

            if numel(varargin) && ~ischar(varargin{1})
                varargin = cat(2, {'CData'}, varargin);
            end

            % TODO: Consider passing only parent parameters here and then
            % set the others later
            self@plugins.dense3D_plugin.Bullseye('Segments', 12, 'AngularOffset', 0, 'CData', 0);

            % Now specify all other parameters
            if numel(varargin)
                set(self, varargin{:})
            end
        end

        function [X, Y] = computeCoordinates(self, sz, varargin)
            % Determine the coordinates for the surface
            if ~isfield(self.cache, 'surfcoords') || ~isequal(self.cache.surfcoords{3}, sz)
                [X, Y] = self.getSampledPoints(self.MinorAxis, self.MajorAxis, ...
                                sz(1) + 1, sz(2) + 1, strcmpi(self.Apex, 'on'));

                X = rot90(X, 2);
                Y = rot90(Y, 2);

                self.cache.surfcoords = {X, Y, sz};
            else
                [X, Y, ~] = self.cache.surfcoords{:};
            end
        end
    end

    methods (Static, Hidden)
        function template = getTemplate(nSegments, sz)

            template = [
                reshape(repmat(1:4, 12, 1), 1, []);
                reshape(repmat(1:4, 12, 1), 1, []) + 4;
                reshape(repmat(1:4, 12, 1), 1, []) + 8;
            ];

            if ismember(nSegments, [13, 19])
                template(end+1,:) = nSegments;
            end

            if exist('sz', 'var')
                template = imresize(template, sz, 'nearest');
            end
        end


        function [X, Y] = getSampledPoints(minor, major, nLayers, nPoints, allin)

            if ~exist('allin', 'var')
                allin = false;
            end

            % Find the RV insertion point
            theta = atan(-sqrt(3) * (major / minor)) + pi;
            thetas = mod([theta, -theta], 2*pi);

            thetaRange = linspace(thetas(1), thetas(2), 100);

            ellipse = getEllipse(major, minor, 0, thetaRange);

            % XXX Assumes that our origin is at 0,0
            bullseyeRadius = norm(ellipse(1,:));
            bullseyeOrigin = [0 0];

            % The origin for everything is [0 0]
            if allin
                shiftRange = [0 5/6];
            else
                shiftRange = [0 2/3];
            end

            shifts = linspace(shiftRange(1),shiftRange(2), nLayers);

            X = nan(numel(shifts), nPoints);
            Y = nan(numel(shifts), nPoints);

            for k = 1:numel(shifts)
                e = getEllipse(major - shifts(k), minor - shifts(k), 0, thetaRange);

                % Now crop it with the circle
                [xx,yy] = ellipseCircleIntersect(e(:,1), e(:,2), ...
                                                 bullseyeOrigin, ...
                                                 bullseyeRadius);


                % Find all of the angles here
                [angles, ~] = cart2pol(xx, yy);
                [angles, inds] = unique(angles, 'stable');
                xx = xx(inds);
                yy = yy(inds);
                angles = mod(angles, 2*pi);

                [~, argmin] = min(angles);
                [~, argmax] = max(angles);

                % Find all the points that are radii
                tokeep = ismember([xx(:), yy(:)], e, 'rows');

                % Ensure that we keep the entire range and don't miss any
                tokeep([argmin argmax]) = true;
                angles = angles(tokeep);

                xx = xx(tokeep);
                yy = yy(tokeep);

                if isempty(angles)
                    continue;
                end

                tt = linspace(min(angles), max(angles), nPoints);
                vals = interp1(angles, [xx(:), yy(:)], tt);

                X(k,:) = vals(:,1);
                Y(k,:) = vals(:,2);
            end

            if allin
                X(end+1,:) = bullseyeOrigin(1) - bullseyeRadius;
                Y(end+1,:) = bullseyeOrigin(2);
            end
        end
    end
end


function ellipse = getEllipse(major, minor, theta, t)
    ellipse(:,2) = minor * sin(t);
    ellipse(:,1) = major * cos(t);

    % Rotate by theta if needed
    ellipse = ellipse * [cos(theta), -sin(theta); sin(theta), cos(theta)];
end

function [xx,yy] = ellipseCircleIntersect(x, y, origin, r)
    ellipse.x = x;
    ellipse.y = y;

    t = linspace(0, 2*pi, 100);

    circle.x = r * cos(t) + origin(1);
    circle.y = r * sin(t) + origin(2);

    p3 = PolygonClip(ellipse, circle, 0);

    if isempty(p3)
        xx = [];
        yy = [];
        return
    end

    ind = (mean(p3(1).x) > 1) + 1;

    xx = p3(ind).x;
    yy = p3(ind).y;
end
