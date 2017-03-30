classdef Bullseye < plugins.dense3D_plugin.HGParrot

    properties (SetObservable)
        AngularOffset   = 2*pi/3    % Angular offset of bullseye
        ButtonDownFcn   = []        % Overloaded ButtonDownFcn
        CData                       % Data to display
        Colormap                    % Custom colormap to use for this
        CLim                        % Color limits
        Data
        Display = 'raw'             % 'raw' | 'average' | 'median'
        LabelFormat = '%0.2f'       % Either a string or func handle
        Radius = 1;                 % Radius to use
        ZData = 0;                  % Height of the plotted object
        Segments = 17               % Number of AHA Segments (16 | 17)
    end

    properties (Hidden)
        haha        % Handle to the AHA lines
        hlistener   % Listener to SelectSegment event
        hsurf       % Handle to the surface object
        htext       % Labels of each of the segments
        htextgroup
        hlisteners
        xdata
        ydata
        zdata
        cache
    end

    properties (Dependent)
        Apex
        SegmentAverages
    end

    events
        SelectSegment
    end

    %--- Get/Set Methods ---%
    methods

        function res = get.Colormap(self)
            if isempty(self.Colormap)
                % Use the colormap of the axes
                res = colormap(ancestor(self.Handle, 'figure'));
            else
                res = self.Colormap;
            end
        end

        function set.CData(self, value)
            if ~isnumeric(value)
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'CData must be numeric')
            end

            self.CData = value;
            self.resetCache();
            self.refresh();
        end

        function set.SegmentAverages(self, value)

            if ~ismember(numel(value), [16 17])
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'There must be either 16 or 17 segment averages.')
            end

            % Compute the CData from this
            template = self.getTemplate(numel(value));
            self.CData = value(template);
            self.Segments = numel(value);
        end

        function res = get.SegmentAverages(self)
            if ~isfield(self.cache, 'segmentAverages')
                self.cache.segmentAverages = self.segment(self.CData, ...
                    self.Segments, @mean);
            end

            res = self.cache.segmentAverages;
        end

        function set.ZData(self, value)
            if ~isnumeric(value) || ~isscalar(value)
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'ZData must be a scalar number')
            end

            self.ZData = self.value;
        end

        function res = get.Apex(self)
            opts = {'off', 'on'};
            tf = ismember([16 17], self.Segments);
            res = opts{tf};
        end

        function set.Segments(self, value)
            if ~isnumeric(value) || ~isscalar(value) || ~ismember(value, [16 17])
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'Segments must be either 16 or 17');
            end

            self.Segments = value;
            self.resetCache();
            self.refresh();
        end

        function set.Apex(self, value)

            if ~ischar(value) || ~any(strcmpi(value, {'on', 'off'}))
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'Apex value must be ''on'' or ''off''')
            end

            if strcmpi(value, 'on')
                self.Segments = 17;
            else
                self.Segments = 16;
            end
        end
    end

    methods
        function self = Bullseye(varargin)

            if numel(varargin) && ~ischar(varargin{1})
                varargin = cat(2, {'CData'}, varargin);
            end

            ip = inputParser();
            ip.KeepUnmatched = true;
            ip.addParamValue('Parent', [], @ishghandle);
            ip.parse(varargin{:})

            if ~isempty(ip.Results.Parent)
                parent = ip.Results.Parent;
            else
                parent = gca;
            end

            self@plugins.dense3D_plugin.HGParrot(hggroup('Parent', parent))
            self.Type = 'Bullseye';

            set(self.Handle, 'ButtonDownFcn', @(s,e)buttonDown(self));

            % Add some linked properties
            initUI(self);

            if numel(varargin)
                set(self, varargin{:});
            end

            self.hlisteners = addlistener(self, ...
                {'LabelFormat', 'AngularOffset', 'Display', 'Radius'}, ...
                'PostSet', @(s,e)refresh(self));

            refresh(self)
        end

        function initUI(self)
            % Surface plot containing the actual data
            self.hsurf = surf(zeros(2), ...
                'Parent',           self.Handle, ...
                'EdgeColor',        'none', ...
                'HitTest',          'on', ...
                'ButtonDownFcn',    @(s,e)buttonDown(self), ...
                'HandleVisibility', 'off', ...
                'Tag',              'bullseye.surf');

            t = linspace(0, 2*pi, 100);
            x = self.Radius * cos(t);
            y = self.Radius * sin(t);

            % AHA Lines
            self.haha = line(x, y, ...
                'Parent',           self.Handle, ...
                'Color',            'k', ...
                'HitTest',          'off', ...
                'ButtonDownFcn',    @(s,e)buttonDown(self), ...
                'HandleVisibility', 'off', ...
                'Tag',              'bullseye.lines');

            % Text objects in the middle of the segments
            self.htext = gobjects(17, 1);

            self.htextgroup = hggroup( ...
                'Parent',           self.Handle, ...
                'ButtonDownFcn',    @(s,e)buttonDown(self), ...
                'Visible',          'off');

            for k = 1:numel(self.htext)
                self.htext(k) = text(NaN, NaN, '', ...
                    'Parent',               self.htextgroup, ...
                    'ButtonDownFcn',        @(s,e)buttonDown(self), ...
                    'HorizontalAlignment',  'center', ...
                    'VerticalAlignment',    'middle', ...
                    'HitTest',              'off', ...
                    'HandleVisibility',     'off', ...
                    'FontUnits',            'normalized', ...
                    'FontSize',             0.05, ...
                    'Tag',                  'bullseye.labels');
            end

            % Link up the necessary attributes here
            addLinkedProperty(self, 'AHA', 'Visible', self.haha);
            addLinkedProperty(self, 'LineColor', 'Color', self.haha);
            addLinkedProperty(self, 'LineWidth', 'LineWidth', self.haha);

            % Text properties
            addLinkedProperty(self, 'FontColor', 'Color', self.htext);
            addLinkedProperty(self, 'FontName', 'FontName', self.htext);
            addLinkedProperty(self, 'FontSize', 'FontSize', self.htext);
            addLinkedProperty(self, 'FontWeight', 'FontWeight', self.htext);
            addLinkedProperty(self, 'FontUnits', 'FontUnits', self.htext);
            addLinkedProperty(self, 'Labels', 'Visible', self.htextgroup);
        end

        function delete(self)
            if ishghandle(self.Handle)
                delete(self.Handle);
            end
        end

        function refresh(self)

            if strcmpi(self.Apex, 'on')
                rlims = [0 self.Radius];
            else
                rlims = [0.25 1] * self.Radius;
            end

            % Redraw the AHA lines
            t = [linspace(0, 2*pi, 100) NaN].';
            rings = (0.25:0.25:1) * self.Radius;

            X = bsxfun(@times, rings, cos(t));
            Y = bsxfun(@times, rings, sin(t));

            % Now worry about the 6 * spokes
            t = linspace(0, 2*pi, 6 + 1).' + self.AngularOffset;
            tmpx = bsxfun(@times, [rings([2 4]) NaN], cos(t)).';
            tmpy = bsxfun(@times, [rings([2 4]) NaN], sin(t)).';

            X = cat(1, X(:), tmpx(:));
            Y = cat(1, Y(:), tmpy(:));

            t = linspace(0, 2*pi, 4 + 1).' + (self.AngularOffset + pi/12);
            tmpx = bsxfun(@times, [rings(1:2) NaN], cos(t)).';
            tmpy = bsxfun(@times, [rings(1:2) NaN], sin(t)).';

            X = cat(1, X(:), tmpx(:));
            Y = cat(1, Y(:), tmpy(:));

            Z = zeros(size(X)) + self.ZData + 0.01;

            set(self.haha, 'XData', X(:), 'YData', Y(:), 'ZData', Z(:));

            sz = max(size(self.CData), [120, 12]);

            % Get the coordinates of the surface
            [X,Y] = self.computeCoordinates(sz, self.AngularOffset, rlims);
            Z = self.ZData * zeros(size(X));

            switch lower(self.Display)
                case 'raw'
                    if ~isempty(self.CData)
                        cdata = imresize(self.CData, sz, 'nearest');
                    else
                        cdata = Z;
                    end
                case 'average'
                    img = self.getTemplate(self.Segments, sz);
                    averages = self.SegmentAverages;
                    cdata = averages(img);
                otherwise
                    error(sprintf('%s:NotImplemented', mfilename), ...
                        '''%s'' display is not implemented yet.', self.Display)
            end

            cdata = rot90(cdata, 2);
            set(self.hsurf, 'XData', X, 'YData', Y, 'ZData', Z, 'CData', cdata);

            % Update the label positions if necessary
            centers = self.segmentCenters;
            centers(:,end+1) = self.ZData + 0.01;

            % If the LabelFormat is a function handle
            if isa(self.LabelFormat, 'function_handle')
                vals = self.segment(self.CData, self.Segments);
                labels = cellfun(self.LabelFormat, vals, 'Uniform', 0);
            else
                func = @(x)sprintf(self.LabelFormat, x);
                labels = arrayfun(func, self.SegmentAverages, 'Uniform', 0);
            end

            if numel(labels) == 16
                labels{end+1} = '';
            end

            set(self.htext, ...
                {'Position'}, num2cell(centers, 2), ...
                {'String'}, labels);

        end
    end

    methods (Access = 'protected')
        function centers = segmentCenters(self)
            % segmentCenters - Compute the center of each AHA segment
            %
            % USAGE:
            %   centers = segmentCenters(self)
            %
            % OUTPUTS:
            %   centers:    [17 x 2] Array, X and Y coordinates for the
            %               center of each segment. If this is a 16-segment
            %               bullseye, then the last X/Y pair is NaN

            thetas = linspace(0, 2*pi, 7).';
            thetas = thetas(1:end-1) - pi/6 + self.AngularOffset;

            x = bsxfun(@times, [7/8 5/8], cos(thetas));
            y = bsxfun(@times, [7/8 5/8], sin(thetas));

            thetas = linspace(0, 2*pi, 5).';
            thetas = thetas(1:end-1) + self.AngularOffset - pi/6;

            x = cat(1, x(:), (3/8) * cos(thetas), 0);
            y = cat(1, y(:), (3/8) * sin(thetas), 0);

            centers = [x, y];

            if self.Segments == 16
                centers(end,:) = NaN;
            end
        end

        function buttonDown(self)
            % 1) Compute which segment the click happened in
            % 2) Create a custom event data object containing this info
            % 3) Trigger the ButtonDownFcn

            % Make sure that hittest is on
            if strcmpi(get(self.Handle, 'HitTest'), 'off')
                return;
            end

            point = get(self.Parent, 'CurrentPoint');

            % Determine which segment was clicked
            data.X = point(1,1);
            data.Y = point(1,2);
            data.Segment = self.whichSegment(data.X, data.Y, ...
                self.Radius, self.AngularOffset);

            import plugins.dense3D_plugin.*

            if ~isempty(self.ButtonDownFcn)
                hgfeval(self.ButtonDownFcn, self, CustomEventData(data));
            end
        end

        function resetCache(self)
            self.cache = struct();
        end
    end

    methods (Static, Hidden)
        function segment = whichSegment(x, y, radius, offset)
            % whichSegment - Determine which segment a point lies in
            %
            % USAGE:
            %   segment = whichSegment(x, y, radius, offset)
            %
            % INPUTS:
            %   x:      Double, X Coordinates of the points to query
            %   y:      Double, Y Coordinates of the points to query
            %   radius: Scalar, Radius of the bullseye
            %   offset: Scalar, Angular offset (in radians) of the bullseye

            % Initialize the output
            segment = nan(size(x));

            [theta, R] = cart2pol(x, y);

            % Determine which circumferential ring we belong to
            ringnum = ceil(R ./ radius ./ 0.25);
            [outer, inds] = ismember(ringnum, [3 4]);

            % Now figure out the circumferential position
            outerthetas = ceil(mod(theta(outer) - offset + pi/3, 2*pi) / (pi/3));
            segment(outer) = (2 - inds(outer)) * 6 + outerthetas;

            inner = ringnum == 2;
            innerthetas = ceil(mod(theta(inner) - offset + 5*pi/12, 2*pi) / (pi/2));
            segment(inner) = innerthetas + 12;

            segment(ringnum == 1) = 17;
        end

        function segments = segment(data, nSegments, operation)
            % segment - Segments the provided data into AHA segment data
            %
            % USAGE:
            %   segments = segment(data, nSegments)
            %   values = segment(data, nSegments, operation)
            %
            % INPUTS:
            %   data:       [M x N] Matrix, Data to be carved up into AHA
            %               segments
            %
            %   nSegments:  Scalar, Indicates the number of AHA segments.
            %               Valid values are either 16 or 17.
            %
            %   operation:  Function Handle (optional), Operation to
            %               perform on the data for each segment
            %
            % OUTPUTS:
            %   segments:   [1 x nSegments] Cell Array, A cell array where
            %               each element contains samples from the input
            %               data that fell into a particular segment.
            %
            %   values:     [1 x nSegments] Array, An array (either a cell
            %               array or numeric array) that is the result of
            %               applying the specified operation to the data
            %               from each segment

            import plugins.dense3D_plugin.*
            template = Bullseye.getTemplate(nSegments);

            % Now scale the image and the mask to the same size
            dims = lcm(size(template), size(data));

            data = imresize(data, dims, 'nearest');
            mask = imresize(template, dims, 'nearest');

            if ~exist('operation', 'var')
                operation = @(x){x};
            end

            segments = accumarray(mask(:), data(:), [], operation);
        end

        function template = getTemplate(nSegments, sz)
            % getTemplate - Method for generating an AHA template image
            %
            % USAGE:
            %   template = Bullseye.getTemplate(nSegments)
            %
            % INPUTS:
            %   nSegments:  Integer, Indicates the number of segments to
            %               use. This must be either 16 or 17.
            %
            %   sz:         [1 x 2] Array, Indicates the desired dimensions
            %               of the output. If not provided, the smallest
            %               possible matrix will be used to accurately
            %               represent the segments.
            %
            % OUTPUTS:
            %   template:   [M x N] Array, A template of the desired size.

            template = [
                reshape(repmat([2:6 1], 12, 1), 1, []);
                reshape(repmat([2:6 1], 12, 1), 1, []) + 6;
                circshift(reshape(repmat(1:4, 18, 1), 1, []), [0 -15]) + 12;
            ];

            if nSegments == 17
                template(end+1,:) = 17;
            end

            if exist('sz', 'var')
                template = imresize(template, sz, 'nearest');
            end
        end

        function [X, Y] = computeCoordinates(sz, offset, rlim)
            % computeCoordinates - Static method to compute coordinates
            %
            %   The surface requires X,Y,Z coordinates of all of the
            %   vertices so this function does all of that calculation,
            %   provided the angular offset, the size of the CData, as well
            %   as the radius limits.

            % Get the dimensions of the CData
            nl = sz(2);
            ns = sz(1);

            % Scale the radii to be in the radius range provided
            r = linspace(rlim(1), rlim(2), ns + 1).';

            % Compute the angle and take into account the offset
            theta = linspace(2*pi + offset, offset, nl + 1);
            theta(end) = [];

            X = r * cos(theta);
            Y = r * sin(theta);

            % Append the first coordinate to the end because of the extra
            % row in the CData
            X = X(:, [1:end 1]);
            Y = Y(:, [1:end 1]);
        end

        function h = demo()
            % demo - Demo showing the basics of using the Bullseye class
            %
            % USAGE:
            %   h = Bullseye.demo()
            %
            % OUTPUTS:
            %   h:  Handle, Handle to the Bullseye object

            import plugins.dense3D_plugin.*
            h = Bullseye(rand(10), ...
                'Labels', 'on', ...
                'ButtonDownFcn', @(s,e)disp('click'));
            axis tight
            axis equal
        end
    end
end
