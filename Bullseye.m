classdef Bullseye < plugins.dense3D_plugin.HGParrot & matlab.mixin.Heterogeneous

    properties (SetObservable)
        AngularOffset   = 2*pi/3    % Angular offset of bullseye
        ButtonDownFcn   = []        % Overloaded ButtonDownFcn
        CData                       % Data to display
        Data
        Display = 'raw'             % 'raw' | 'average' | 'median'
        LabelFormat = '%0.2f'       % Either a string or func handle
        Radius = 1;                 % Radius to use
        ZData = 0;                  % Height of the plotted object
        Segments                    % Number of AHA Segments (16 | 17)
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

        alphadata_
        colormap_
        clim_
    end

    properties (Dependent)
        AlphaData                   % Transparency
        Apex
        Colormap                    % Custom colormap to use for this
        CLim                        % Color limits
        SegmentAverages
    end

    events
        SelectSegment
    end

    %--- Get/Set Methods ---%
    methods

        function res = get.Colormap(self)
            if isempty(self.colormap_)
                % Use the colormap of the axes
                res = colormap(ancestor(self.Handle, 'figure'));
            else
                res = self.colormap_;
            end
        end

        function res = get.AlphaData(self)
            res = self.alphadata_;
        end

        function set.AlphaData(self, value)
            self.alphadata_ = value;
            refresh(self);
        end

        function set.Colormap(self, value)
            self.colormap_ = value;
            self.CLim = self.CLim;
        end

        function res = get.CLim(self)
            if isempty(self.clim_)
                res = get(ancestor(self.Handle, 'axes'), 'CLim');
            else
                res = self.clim_;
            end
        end

        function set.CLim(self, value)
            if isempty(self.colormap_)
                set(ancestor(self.Handle, 'axes'), 'CLim', value);
            else
                self.clim_ = value;
                refresh(self);
            end
        end

        function set.CData(self, value)

            if isequal(value, self.CData)
                return
            end

            if ~isnumeric(value)
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'CData must be numeric')
            end

            self.CData = value;
            self.resetCache('CData');
            self.refresh();
        end

        function set.SegmentAverages(self, value)

            if ~ismember(numel(value), self.validSegments())
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'Invalid number of segments specified')
            end

            % Compute the CData from this
            template = self.getTemplate(numel(value));
            self.CData = value(template);
            self.Segments = numel(value);
        end

        function res = get.SegmentAverages(self)
            if ~isfield(self.cache, 'segmentAverages')
                self.cache.segmentAverages = self.segment(@mean);
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
            if self.hasApex()
                res = 'on';
            else
                res = 'off';
            end
        end

        function set.Segments(self, value)
            if isequal(value, self.Segments)
                return
            end

            self.Segments = value;
            self.resetCache('Segments');
            self.refresh();
        end

        function set.Apex(self, value)

            if ~ischar(value) || ~any(strcmpi(value, {'on', 'off'}))
                error(sprintf('%s:InvalidValue', mfilename), ...
                    'Apex value must be ''on'' or ''off''')
            end

            self.Segments = floor(self.Segments / 2) * 2 + strcmpi(value, 'on');
        end
    end

    methods (Sealed)
        function disp(self)
            disp@plugins.dense3D_plugin.HGParrot(self)
        end
    end

    methods (Access = 'protected')
        function bool = hasApex(self)
            bool = self.Segments == 17;
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

            % Set the default segment
            segments = self.validSegments();
            self.Segments = segments(1);

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

            % AHA Lines
            self.haha = line(NaN, NaN, NaN, ...
                'Parent',           self.Handle, ...
                'Color',            'k', ...
                'HitTest',          'off', ...
                'ButtonDownFcn',    @(s,e)buttonDown(self), ...
                'HandleVisibility', 'off', ...
                'Tag',              'bullseye.lines');

            % Text objects in the middle of the segments
            self.htext = gobjects(max(self.validSegments), 1);

            self.htextgroup = hggroup( ...
                'Parent',           self.Handle, ...
                'ButtonDownFcn',    @(s,e)buttonDown(self), ...
                'Visible',          'off');

            for k = 1:numel(self.htext)
                self.htext(k) = text(NaN, NaN, '', ...
                    'Parent',               self.htextgroup, ...
                    'ButtonDownFcn',        @(s,e)buttonDown(self), ...
                    'Visible',              'off', ...
                    'HorizontalAlignment',  'center', ...
                    'VerticalAlignment',    'middle', ...
                    'HitTest',              'off', ...
                    'HandleVisibility',     'on', ...
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
            addLinkedProperty(self, 'Labels', 'Visible', [self.htextgroup; self.htext(:)]);
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
            [X, Y] = self.ahalines();
            Z = zeros(size(X)) + self.ZData;

            % Need to apply a z offset so that OpenGL renderer doesn't
            % screw things up on HG1
            if ~ishg2
                Z = Z + 0.6;
            end

            set(self.haha, 'XData', X(:), 'YData', Y(:), 'ZData', Z(:));

            sz = max(size(self.CData), [40, 120]);

            % Get the coordinates of the surface
            [X,Y] = self.computeCoordinates(sz, self.AngularOffset, rlims);
            Z = self.ZData * zeros(size(X));

            switch lower(self.Display)
                case 'raw'
                    if ~isempty(self.CData)
                        cdata = self.CData;

                        if ~isempty(self.AlphaData)
                            adata = self.AlphaData;
                        else
                            adata = 1;
                        end
                    else
                        cdata = Z;
                        adata = 1;
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

            % Scale the colordata if necessary
            if ~isempty(self.colormap_)
                % Check if we need to apply scaling
                if ~isempty(self.CLim)
                    cdata = (cdata - self.CLim(1)) ./ self.CLim(2);

                    % Chop the upper/lower values
                    cdata(cdata < 0) = 0;
                    cdata(cdata > 1) = 1;
                end

                ind = gray2ind(cdata, size(self.colormap_, 1));
                cdata = ind2rgb(ind, self.colormap_);
            end

            inputs = {};

            if exist('adata', 'var')
                adata = rot90(adata, 2);
                if isscalar(adata)
                    inputs = [inputs, {'FaceAlpha', adata}];
                else
                    inputs = [inputs, {'AlphaData',         adata, ...
                                       'FaceAlpha',         'texturemap', ...
                                       'AlphaDataMapping',  'none'}];
                end
            else
                inputs = [inputs, {'FaceAlpha', 1}];
            end

            set(self.hsurf, ...
                'XData', X, ...
                'YData', Y, ...
                'ZData', Z, ...
                'FaceColor', 'texturemap', ...
                'CData', cdata, ...
                inputs{:});

            % Update the label positions if necessary
            centers = self.segmentCenters;
            centers(:,end+1) = self.ZData + 0.01;

            % If the LabelFormat is a function handle
            if isa(self.LabelFormat, 'function_handle')
                vals = self.segment();
                labels = cellfun(self.LabelFormat, vals, 'Uniform', 0);
            else
                func = @(x)sprintf(self.LabelFormat, x);
                labels = arrayfun(func, self.SegmentAverages, 'Uniform', 0);
            end

            if strcmpi(self.Apex, 'off')
                labels{end+1} = '';
            end

            set(self.htext, ...
                {'Position'}, num2cell(centers, 2), ...
                {'String'}, labels);

        end
    end

    methods (Access = 'protected')

        function res = validSegments(~)
            res = [17, 16];
        end

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

        function [X, Y] = ahalines(self)

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
            data.Segment = self.whichSegment(data.X, data.Y);

            import plugins.dense3D_plugin.*

            if ~isempty(self.ButtonDownFcn)
                hgfeval(self.ButtonDownFcn, self, CustomEventData(data));
            end
        end

        function resetCache(self, propname)
            self.cache = struct();
        end
    end

    methods (Hidden)
        function segments = segment(self, operation)
            % segment - Segments the provided data into AHA segment data
            %
            % USAGE:
            %   segments = segment()
            %   values = segment(operation)
            %
            % INPUTS:
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

            template = self.getTemplate(self.Segments);
            data = self.CData;

            if isempty(data)
                data = NaN;
            end

            % Now scale the image and the mask to the same size
            dims = lcm(size(template), size(data));

            data = imresize(data, dims, 'nearest');
            mask = imresize(template, dims, 'nearest');

            if ~exist('operation', 'var')
                operation = @(x){x};
            end

            segments = accumarray(mask(:), data(:), [], operation);
        end

        function [X, Y] = computeCoordinates(~, sz, offset, rlim)
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

        function segment = whichSegment(self, x, y)
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


            X = get(self.hsurf, 'XData');
            Y = get(self.hsurf, 'YData');

            template = rot90(self.getTemplate(self.Segments, size(X)), 2);
            [~, ind] = min((X(:) - x).^2 + (Y(:) - y).^2);

            segment = template(ind);
        end
    end

    methods (Static, Hidden)

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
