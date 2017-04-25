classdef DENSE3Dviewer < DataViewer

    % Data to show
    %
    %   Mesh (3D Display)
    %   Polar Strains (Global curve, regional curves, bullseye)
    %   Principal Strains (Global curve, regional curves, bullseye)
    %   Shear Strains (Global curve, regional curves, bullseye)
    %   Torsion (Global curve, regional curves, bullseye)
    %   Peak Polar Strains (Bullseye)
    %   Peak Principal Strains (Bullseye)
    %   Peak Shear Strains (Bullseye)
    %   CURE/RURE/LURE (global curve)
    %   Delay Times (Bullseye)

    properties
        Handles
        Data
        ActiveSegments = {true(1, 17); true(1, 12)};
        LineAppearance = struct('LineWidth', 2);
        AxesAppearance = struct( ...
            'LineWidth',    2, ...
            'FontWeight',   'bold', ...
            'FontSize',     14);
        BullseyeAxesAppearance = struct( ...
            'LineWidth',    2, ...
            'FontWeight',   'bold', ...
            'FontUnits',    'normalized', ...
            'FontSize',     0.03)
    end

    properties (SetAccess = 'private')
        isAllowExportMat = false;
        isAllowExportExcel = false;
    end

    properties (Hidden)
        api
        options

        template = struct( ...
            'initFcn',      [], ...
            'resizeFcn',    [], ...
            'playbackFcn',  [], ...
            'deleteFcn',    []);

        dispclr = [78, 101, 148] / 255;

        cache = struct()
        hlistener
    end

    properties (Hidden)
        % Provides a mapping from each type of strain to it's name,
        % scaling, display format, and peak determination function
        mappings = struct( ...
            'RR',       {{'Radial Strain (%)',                  @(x)x * 100,        '%.0f', @max}}, ...
            'CC',       {{'Circum. Strain (%)',                 @(x)x * 100,        '%.0f', @min}}, ...
            'LL',       {{'Long. Strain (%)',                   @(x)x * 100,        '%.0f', @min}}, ...
            'p1',       {{'1st Principal Strain (%)',           @(x)x * 100,        '%.0f', @max}}, ...
            'p2',       {{'2nd Principal Strain (%)',           @(x)x * 100,        '%.0f', @min}}, ...
            'p3',       {{'3rd Principal Strain (%)',           @(x)x * 100,        '%.0f', @min}}, ...
            'RL',       {{'RL Shear Strain (%)'                 @(x)x * 100,        '%.0f', @min}}, ...
            'CL',       {{'CL. Shear Strain (%)',               @(x)x * 100,        '%.0f', @min}}, ...
            'CR',       {{'CR Shear Strain (%)',                @(x)x * 100,        '%.0f', @min}}, ...
            'CURE',     {{'CURE',                               @(x)x,              '%.3f', @min}}, ...
            'RURE',     {{'RURE',                               @(x)x,              '%.3f', @min}}, ...
            'LURE',     {{'LURE',                               @(x)x,              '%.3f', @min}}, ...
            'TORSION',  {{'Torsion (degrees)',                  @(x)x.CLShearAngle, '%.1f', @max}}, ...
            'DELAY',    {{'Regional Delay (%CC)',               @(x)x.DelayTimes * 100,           '%.1f', @min}});
    end

    properties (Dependent = true)
        Frames
    end

    properties (SetAccess = 'private', GetAccess = 'public')

    end

    methods
        function res = get.Frames(self)
            AI = [self.Data.Data.AnalysisInfo];

            frames = arrayfun(@(x)x.FramesForAnalysis(:).', AI, 'Uniform', 0);
            frames = cat(1, frames{:});
            res = [max(frames(:,1)), min(frames(:,2))];
        end

        function strains = regionalStrains(self)
            function result = getStrains()
                result = self.Data.regionalStrains();

                % If this is BiVentricular remove the segments we don't
                % want to keep
                if numel(result) > 1
                    allfields = fieldnames(result);
                    allfields(ismember(allfields, {'CURE', 'RURE', 'LURE'})) = [];

                    tokeep = [1:4, 7:10, 13:16];

                    for m = 1:numel(allfields)
                        result(2).(allfields{m}) = result(2).(allfields{m})(tokeep,:);
                    end
                end
            end

            strains = self.fetchCache('regionalStrain', @getStrains);
        end

        function clearCache(self)
            self.cache = struct();
        end

        function res = fetchCache(self, key, func)
            if ~isfield(self.cache, key)
                self.cache.(key) = func();
            end

            res = self.cache.(key);
        end

        function newDataCallback(self)
            self.clearCache();

            % Check to see if we want to enable the buttons
            if ~isempty(self.Data.Data)
                set([self.options.Handle], 'Enable', 'on')
            else
                set([self.options.Handle], 'Enable', 'off')
            end

            set([self.options.Handle], 'Value', 0)

            % Now hide everything
            reset(self)

            self.isAllowExportMat = ~isempty(self.Data);
        end

        function varargin = exportExcel(varargin)
            error(sprintf('%s:NotImplemented', mfilename), ...
                'exportExcel not implemented yet');
        end

        function [file, data] = exportMat(obj, startpath)

            if ~exist('startpath', 'var')
                startpath = pwd;
            end

            if ~exist(startpath, 'dir')
                [origstartpath, ~, ext] = fileparts(startpath);
                file = startpath;
                startpath = origstartpath;

                if isempty(ext)
                    file = [file, '.mat'];
                end
            end

            if ~exist('file', 'var')
                [uifile, uipath] = uiputfile('*.mat', [], startpath);
                if isequal(uifile, 0)
                    file = [];
                    return
                end

                file = fullfile(uipath, uifile);
                [~, ~, e] = fileparts(file);

                if ~isequal(e, '.mat')
                    file = [file, '.mat'];
                end
            end

            if isempty(obj.Data.Strains)
                obj.computeStrains();
            end

            hwait = waitbartimer();
            hwait.String = 'Generating output file...';
            hwait.WindowStyle = 'modal';
            hwait.AllowClose = false;
            hwait.Visible = 'on';
            start(hwait);

            drawnow

            hclean = onCleanup(@()delete(hwait(isvalid(hwait))));

            data = obj.Data.save();
            save(file, '-struct', 'data')
        end


        function self = DENSE3Dviewer(data, varargin)

            fakedata = DENSEdata;
            self = self@DataViewer(struct([]), fakedata, varargin{:});

            % Assign the data object that we actually care about here
            self.Data = data;

            self.hlistener = addlistener(data, 'NewData', @(s,e)self.newDataCallback());

            % Now set stuff up like we care
            self.Handles.hbuttongroup = uibuttongroup( ...
                'Parent', self.hcontrol, ...
                'BorderType', 'none', ...
                'Units', 'normalized', ...
                'position', [0 0 1 1], ...
                'SelectionChangeFcn', @(s,e)self.changeDisplay(e));

            fcn = @(str, varargin)uicontrol( ...
                'Parent', self.Handles.hbuttongroup, ...
                'Style', 'radiobutton', ...
                'Units', 'normalized', ...
                'Enable', 'off', ...
                'Value', 0, ...
                'String', str, varargin{:});

            buttons = struct([]);

            buttons(end+1).String = '3D Mesh';
            buttons(end).Fcn = @(s,e)showMesh(self);

            buttons(end+1).String = 'Global Polar Strains';
            buttons(end).Fcn = @(s,e)showCurves(self, {'RR', 'CC', 'LL'});

            buttons(end+1).String = 'Global Principal Strains';
            buttons(end).Fcn = @(s,e)showCurves(self, {'p1', 'p2', 'p3'});

            buttons(end+1).String = 'Global Shear Strains';
            buttons(end).Fcn = @(s,e)showCurves(self, {'CL', 'CR', 'RL'});

            buttons(end+1).String = 'Global Torsion';
            buttons(end).Fcn = @(s,e)showCurves(self, {'TORSION', 'CL'});

            buttons(end+1).String = 'Polar Strain Bullseye';
            buttons(end).Fcn = @(s,e)showBullseye(self, {'RR', 'CC', 'LL'});

            buttons(end+1).String = 'Principal Strain Bullseye';
            buttons(end).Fcn = @(s,e)showBullseye(self, {'p1', 'p2', 'p3'});

            buttons(end+1).String = 'Shear Strain Bullseye';
            buttons(end).Fcn = @(s,e)showBullseye(self, {'CL', 'CR', 'RL'});

            buttons(end+1).String = 'Torsion Bullseye';
            buttons(end).Fcn = @(s,e)showBullseye(self, {'TORSION'});

            buttons(end+1).String = 'Regional Polar Strain';
            buttons(end).Fcn = @(s,e)showCurves(self, {'RR', 'CC', 'LL'}, 1);

            buttons(end+1).String = 'Regional Principal Strain';
            buttons(end).Fcn = @(s,e)showCurves(self, {'p1', 'p2', 'p3'}, 1);

            buttons(end+1).String = 'Regional Shear Strain';
            buttons(end).Fcn = @(s,e)showCurves(self, {'CL', 'CR', 'RL'}, 1);

            buttons(end+1).String = 'Regional Torsion';
            buttons(end).Fcn = @(s,e)showCurves(self, {'TORSION', 'CL'}, 1);

            buttons(end+1).String = 'Peak Polar Strain';
            buttons(end).Fcn = @(s,e)showBullseye(self, {'RR', 'CC', 'LL'} , false);

            buttons(end+1).String = 'Peak Principal Strain';
            buttons(end).Fcn = @(s,e)showBullseye(self, {'p1', 'p2', 'p3'} , false);

            buttons(end+1).String = 'Peak Torsion';
            buttons(end).Fcn = @(s,e)showBullseye(self, {'TORSION'}, false);

            buttons(end+1).String = 'CURE / RURE / LURE';
            buttons(end).Fcn = @(s,e)showCurves(self, {'CURE', 'RURE', 'LURE'} , false);

            buttons(end+1).String = 'Regional Delay Times';
            buttons(end).Fcn = @(s,e)showBullseye(self, 'DELAY', false);

            bottoms = linspace(0, 1, numel(buttons) + 1);
            bottoms(end) = [];

            buttons = flip(buttons);

            for k = 1:numel(buttons)
                pos = [0.05 bottoms(k) 0.9 diff(bottoms(1:2))];
                buttons(k).Handle = fcn(buttons(k).String, 'Position', pos);
            end

            self.options = buttons;

            set([self.options.Handle], 'Value', 0);

            self.api = self.template;

            self.redrawenable = true;
            redraw(self);

            set(self.hdisplay, 'BackgroundColor', 'k')
        end

        function api = showBullseye(self, type, dynamic)

            if ~exist('dynamic', 'var')
                dynamic = true;
            end

            if isempty(self.Data.Strains)
                self.computeStrains();
            end

            if ischar(type)
                type = {type};
            end


            hax = NaN;
            hpan = NaN;
            hpanning = pan(gcbf);
            hzoom = zoom(gcbf);
            hrot = rotate3d(gcbf);

            import plugins.dense3D_plugin.*
            hbull = Bullseye.empty();

            api.initFcn = @initFcn;

            if dynamic
                api.playbackFcn = @playbackFcn;
            else
                api.playbackFcn = [];
            end

            strains = self.regionalStrains();

            api.deleteFcn = @deleteFcn;
            api.resizeFcn = @resizeFcn;

            function resizeFcn()
                pos = getpixelposition(self.hdisplay);

                pos(2) = pos(2) + 50;
                pos(4) = pos(4) - 60;
                pos(1) = 0;
                pos(3) = diff(pos([1 3]));

                setpixelposition(hpan, pos);
            end


            function deleteFcn()
                if ishandle(hpan)
                    delete(hpan)
                end
            end

            function playbackFcn()
                fr = self.hplaybar.Value;

                for k = 1:numel(type)
                    values = strains.(type{k});
                    mx = max(abs(values(:)));
                    values = values(:, fr);

                    set(hax(k), 'CLim', [-mx mx])
                    set(hbull(k), 'SegmentAverages', values);
                end
            end

            function initFcn()
                import plugins.dense3D_plugin.*
                self.exportaxes = true;

                set(self.hfigure_display, 'renderer', 'opengl')

                hpan = uipanel( ...
                    'Parent', self.hdisplay, ...
                    'BorderType', 'none', ...
                    'BackgroundColor', get(self.hdisplay, 'BackgroundColor'));

                nAxes = numel(type);

                for k = 1:nAxes
                    hax(k) = axes( ...
                        'Parent',       hpan, ...
                        self.BullseyeAxesAppearance);

                    set(hax(k), ...
                        'LooseInset',       [0.05 0.05 0.05 0.05], ...
                        'OuterPosition',    [(k-1)/nAxes, 0 (1/nAxes) 1])

                    mapping = self.mappings.(type{k});

                    if ~isfield(strains, type{k})
                        strains.(type{k}) = mapping{2}(strains);
                    else
                        strains.(type{k}) = mapping{2}(strains.(type{k}));
                    end

                    title(hax(k), mapping{1}, 'Color', self.dispclr)
                    values = mapping{4}(strains.(type{k}), [], 2);

                    hbull(k) = Bullseye( ...
                        'Parent', hax(k), ...
                        'SegmentAverages', values, ...
                        'Labels', 'on', ...
                        'LabelFormat', mapping{3});

                    %set(hbull, 'UIContextMenu', context);

                    mx = max(abs(values));
                    set(hax(k), 'CLim', [-mx mx], 'xtick', [], 'ytick', []);

                    if ishg2
                        cb = colorbar(handle(hax(k)));
                        set(cb, 'Color', self.dispclr);
                    else
                        cb = colorbar('peer', hax(k));
                        set(cb, 'XColor', self.dispclr, ...
                            'YColor', self.dispclr)
                    end

                    set(cb, ...
                            'FontWeight', 'bold', ...
                            'FontSize', 12, ...
                            'LineWidth', 2, ...
                            'Location', 'southoutside')
                end

                hpanning = pan(gcbf);
                hzoom = zoom(gcbf);
                hrot = rotate3d(gcbf);

                hpanning.setAllowAxesPan(hax, false);
                hzoom.setAllowAxesZoom(hax, false);
                hrot.setAllowAxesRotate(hax, false);
                hrot.Enable = 'off';

                colormap(bwr(65))

                if dynamic
                    frames = self.Frames;
                    self.hplaybar.Min = frames(1);
                    self.hplaybar.Max = frames(2);
                    self.hplaybar.Enable = 'on';
                    self.hlisten_playbar.Enabled = true;
                    playbackFcn();
                else
                    self.hplaybar.Enable = 'off';
                    self.hlisten_playbar.Enabled = false;
                end

                axis(hax, 'equal')
                axis(hax, 'tight')
                axis(hax, 'off')

                % Allow exports
                self.isAllowExportImage = true;
                if (self.hplaybar.Max - self.hplaybar.Min) > 1
                    self.isAllowExportVideo = true;
                end
            end
        end

        function computeDisplacementSplines(self)

            hwait = waitbartimer();
            hwait.String = 'Fitting Radial Basis Functions...';
            hwait.WindowStyle = 'modal';
            hwait.AllowClose = false;
            hwait.Visible = 'on';
            start(hwait);

            drawnow

            hclean = onCleanup(@()delete(hwait(isvalid(hwait))));

            function statusFcn(x, total)
                fmt = 'Fitting Radial Basis Functions... (%d/%d)';
                hwait.String = sprintf(fmt, x, total);
                drawnow
            end

            self.Data.interpolate(@statusFcn)
        end


        function computeStrains(self)

            hwait = waitbartimer();
            hwait.String = 'Computing Strains...';
            hwait.WindowStyle = 'modal';
            hwait.AllowClose = false;
            hwait.Visible = 'on';
            start(hwait);

            drawnow

            hclean = onCleanup(@()delete(hwait(isvalid(hwait))));

            if isempty(self.Data.EndocardialMesh)
                hwait.String = 'Rendering 3D Mesh...';
                drawnow
                self.Data.generateMeshes();
            end

            function statusFcn(x, total)
                fmt = 'Fitting Radial Basis Functions... (%d/%d)';
                hwait.String = sprintf(fmt, x, total);
                drawnow
            end

            if isempty(self.Data.Interpolants)
                self.Data.interpolate(@statusFcn);
            end

            hwait.String = 'Parameterizing the ventricle...';

            self.Data.parameterize();

            drawnow

            hwait.String = 'Computing Strains...';

            self.Data.computeStrains();

            % Go ahead and fetch the regional strains
            drawnow
        end

        function api = showCurves(self, type, regional)

            if ~exist('regional', 'var')
                regional = false;
            end

            if isempty(self.Data.Strains)
                self.computeStrains();
            end

            if ischar(type)
                type = {type};
            end

            if regional
                strains = self.regionalStrains();
            else
                strains = self.Data.Strains;
            end

            isBiV = numel(strains) > 1;

            hax = NaN;
            hpan = NaN;
            hplot = {};
            hbullax = NaN;
            hbull = [];
            hrvbull = [];

            api.initFcn = @initFcn;
            api.playbackFcn = @playbackFcn;
            api.deleteFcn = @deleteFcn;
            api.resizeFcn = @resizeFcn;

            function playbackFcn()

                if isempty(hbull) || ~isvalid(hbull)
                    return
                end

                for k = 1:size(hplot, 2)
                    nPlots = numel(hplot{1,k});

                    colors = num2cell(hsv(nPlots), 2);

                    visible = repmat({'on'}, nPlots, 1);
                    visible(~self.ActiveSegments{k}) = {'off'};

                    if k == 1
                        bullhandle = hbull;
                    else
                        bullhandle = hrvbull;
                    end

                    cdata = bullhandle.CData;
                    tf = ismember(cdata, find(~self.ActiveSegments{k}));

                    % Set the disabled segments to be partly transparent
                    alphadata = ones(size(tf));
                    alphadata(tf) = 0.3;

                    set(bullhandle, 'AlphaData', alphadata);

                    for m = 1:size(hplot, 1)
                        set(hplot{m,k}, {'Color'}, colors, {'Visible'}, visible)
                    end
                end
            end

            function initFcn()
                self.exportaxes = true;
                set(self.hfigure_display, 'Renderer', 'opengl')

                hpan = uipanel( ...
                    'Parent', self.hdisplay, ...
                    'BorderType', 'none', ...
                    'BackgroundColor', get(self.hdisplay, 'BackgroundColor'));

                nAxes = numel(type);

                hflow = uiflowcontainer('v0', ...
                    'Parent', hpan, ...
                    'FlowDirection', 'lefttoright', ...
                    'BackgroundColor', get(self.hdisplay, 'BackgroundColor'));

                hzoom = zoom(gcbf);
                hpanning = pan(gcbf);
                hrot = rotate3d(gcbf);

                if regional

                    if isBiV
                        xlims = [-2, 1.2];
                        ylims = [-1.2, 1.2];
                    else
                        xlims = [-1.2, 1.2];
                        ylims = xlims;
                    end

                    hbullax = axes( ...
                        'Parent', hflow, ...
                        'XTick', [], ...
                        'YTick', [], ...
                        'XLim', xlims, ...
                        'YLim', ylims, ...
                        'Color', 'none', ...
                        'XColor', get(hflow, 'BackgroundColor'), ...
                        'YColor', get(hflow, 'BackgroundColor'), ...
                        'Visible', 'on');

                    hzoom.setAllowAxesZoom(hbullax, false);
                    hpanning.setAllowAxesPan(hbullax, false);
                    hrot.setAllowAxesRotate(hbullax, false);

                    set(hbullax, 'WidthLimits', [250, 250])

                    title(hbullax, 'Segment Model', ...
                        'Color', self.dispclr, ...
                        'FontWeight', 'bold');

                    import plugins.dense3D_plugin.*

                    hbull = Bullseye( ...
                        'SegmentAverages', 1:17, ...
                        'Labels', 'on', ...
                        'LabelFormat', '%d', ...
                        'Parent', hbullax, ...
                        'Colormap', hsv, ...
                        'CLim', [1 17], ...
                        'ButtonDownFcn', @(s,e)clickBullseye(e, 1));

                    if isBiV
                        hrvbull = RVBullseye( ...
                            'SegmentAverages', 1:12, ...
                            'Labels', 'on', ...
                            'LabelFormat', '%d', ...
                            'Parent', hbullax, ...
                            'Colormap', hsv, ...
                            'CLim', [1 12], ...
                            'ButtonDownFcn', @(s,e)clickBullseye(e, 2));
                    end

                    axis(hbullax, 'equal')
                    axis(hbullax, 'tight')

                    % Add some instructions
                    msg = {
                        'Click on a segment to toggle or'
                        'right-click to toggle all segments'
                    };
                    text(mean(xlims), -1.3, msg, 'Parent', hbullax, 'HorizontalAlignment', 'center', ...
                        'Color', self.dispclr);

                    colormap(gcbf, hsv);
                end

                subpan = uipanel( ...
                    'Parent', hflow, ...
                    'BorderType', 'none', ...
                    'BackgroundColor', get(self.hdisplay, 'BackgroundColor'));

                hplot = cell(nAxes, numel(strains));
                hax = gobjects(nAxes, 1);

                colors = [self.dispclr; self.dispclr];
                linestyles = {'-', '--'};

                for k = 1:nAxes

                    mapping = self.mappings.(type{k});
                    hax(k) = subplot(nAxes, 1, k, 'Parent', subpan);
                    set(hax(k), self.AxesAppearance);

                    ylabel(hax(k), mapping{1});
                    xlabel(hax(k), 'Frame')

                    for n = 1:numel(strains)
                        strain = strains(n);

                        if ~isfield(strain, type{k})
                            strain.(type{k}) = mapping{2}(strain);
                        else
                            strain.(type{k}) = mapping{2}(strain.(type{k}));
                        end

                        values = strain.(type{k});

                        if ~regional
                            values = mean(values(~any(isnan(values), 2), :), 1);
                        end

                        hplot{k,n} = line(1:size(values, 2), values, 'Parent', hax(k), self.LineAppearance);

                        nPlots = numel(hplot{k,n});

                        if nPlots > 1
                            set(hplot{k,n}, {'Color'}, num2cell(hsv(nPlots), 2), 'LineStyle', linestyles{n});
                        else
                            set(hplot{k,n}, 'Color', colors(n,:), 'LineStyle', linestyles{n})
                        end
                    end

                    if numel(strains) > 1 && ~regional
                        hleg = legend({'LV', 'RV'});
                        set(hleg, 'EdgeColor', 'none', ...
                            'Color', 'none', 'TextColor', self.dispclr)
                    end
                end

                set(hax, ...
                    'Color', 'none', ...
                    'box',  'on', ...
                    'xcolor', self.dispclr(1,:), ...
                    'xlim', [1, 1 + diff(self.Frames)], ...
                    'ycolor', self.dispclr(1,:));

                hzoom.setAllowAxesZoom(hax, false);
                hpanning.setAllowAxesPan(hax, false);
                hrot.setAllowAxesRotate(hax, false);
                hrot.Enable = 'off';

                self.isAllowExportImage = true;
                playbackFcn();
            end

            function clickBullseye(evnt, index)
                switch lower(get(gcf, 'SelectionType'))
                    case 'alt'
                        sz = size(self.ActiveSegments{index});
                        if ~all(self.ActiveSegments{index})
                            self.ActiveSegments{index} = true(sz);
                        else
                            self.ActiveSegments{index} = false(sz);
                        end
                    otherwise
                        self.ActiveSegments{index}(evnt.Segment) = ~self.ActiveSegments{index}(evnt.Segment);
                end


                playbackFcn();
            end

            function resizeFcn()
                % Get the position of the parent object
                pos = getpixelposition(self.hdisplay);

                pos(2) = pos(2) + 50;
                pos(4) = pos(4) - 60;
                pos(1) = 0;

                pos(3) = diff(pos([1 3]));

                setpixelposition(hpan, pos)
            end

            function deleteFcn()

                if ishandle(hpan)
                    delete(hpan)
                end

                if ishandle(hbullax)
                    delete(hbullax)
                end
            end
        end

        function api = showMesh(self)

            hax = NaN;
            htitle = NaN;
            hendo = NaN;
            hepi = NaN;

            meshes = cell(0,1);

            api.initFcn     = @initFcn;
            api.playbackFcn = @playbackFcn;
            api.deleteFcn   = @deleteFcn;
            api.resizeFcn   = @resizeFcn;

            function initFcn()
                self.exportaxes = true;

                set(self.hfigure_display, 'renderer', 'opengl')

                hax = axes('Parent', self.hdisplay);
                htitle = title(hax, '3D Mesh', 'Color', self.dispclr);

                set(hax, ...
                    'Color', 'none', ...
                    'box',  'on', ...
                    'Visible', 'off', ...
                    'xtick', [], ...
                    'ytick', [], ...
                    'Clipping', 'off', ...
                    'FontWeight', 'bold', ...
                    'xcolor', self.dispclr(1,:), ...
                    'zcolor', self.dispclr(1,:), ...
                    'ycolor', self.dispclr(1,:));

                grid(hax, 'on')

                if isempty(self.Data.EpicardialMesh)
                    self.Data.generateMeshes();

                    hwait = waitbartimer();
                    hwait.String = 'Rendering 3D Mesh...';
                    hwait.WindowStyle = 'modal';
                    hwait.AllowClose = false;
                    start(hwait);

                    clean = onCleanup(@(x)delete(hwait(isvalid(hwait))));
                    pause(3)
                end

                meshes = {[self.Data.EpicardialMesh, self.Data.EndocardialMesh]};

                for k = 1:numel(meshes)
                    hmeshes(k) = patch(meshes(k), 'FaceColor', 'w', 'FaceAlpha', 0.5, 'Parent', hax);
                    hold(hax, 'on')
                end

                axis(hax, 'equal');

                set(htitle, 'FontSize', 12, ...
                    'Color', self.dispclr(1,:), ...
                    'FontWeight', 'bold');

                view(hax, 3);

                set(hax, 'Position', [0 0 1 1])

                % Enable rotation
                hrot = rotate3d(self.hfigure_display);
                hrot.setAllowAxesRotate(hax, true);
                hrot.Enable = 'on';

                % Disable Contrast
                self.hcontrast.setAllowAxes(hax, false);

                % Set base zoom level
                zoom(hax, 'reset');

                frames = self.Frames;
                self.hplaybar.Min = frames(1);
                self.hplaybar.Max = frames(2);
                self.hplaybar.Enable = 'on';

                playbackFcn();

                self.hlisten_playbar.Enabled = true;

                self.isAllowExportImage = true;
                if (self.hplaybar.Max - self.hplaybar.Min) > 1
                    self.isAllowExportVideo = true;
                end
            end

            function playbackFcn()
                frame = self.hplaybar.Value;

                if numel(meshes) < frame || isempty(meshes{frame})

                    if isempty(self.Data.Interpolants)
                        % Compute displacement splines as needed
                        self.computeDisplacementSplines()
                    end

                    interp = self.Data.Interpolants(frame);

                    % Interpolate the epicardial surface
                    epi = self.Data.EpicardialMesh;
                    epi.vertices = epi.vertices + interp.query(epi.vertices);

                    endo = self.Data.EndocardialMesh;
                    endo.vertices = endo.vertices + interp.query(endo.vertices);

                    meshes{frame} = [epi, endo];
                end

                msh = meshes{frame};

                set(hepi, msh(1))
                set(hendo, msh(2));

                drawnow
            end

            function deleteFcn()
                delete(hax(ishandle(hax)));
            end
        end

        function redraw(self)
            redraw@DataViewer(self);
        end

        function changeDisplay(self, evnt)

            if isequal(evnt.NewValue, evnt.OldValue)
                return;
            end

            % Reset
            reset(self)

            [tf, ind] = ismember(evnt.NewValue, [self.options.Handle]);

            if ~tf
                return;
            end

            % Initialize the object
            try
                API = self.options(ind).Fcn(self);
                API.initFcn();
            catch ERR
                if exist('API', 'var') && isstruct(API) && ~isempty(API.deleteFcn)
                    API.deleteFcn();
                end

                self.api = self.template;

                rethrow(ERR);
            end

            self.api = API;

            if isempty(API.playbackFcn)
                self.hplaybar.Visible = 'off';
            else
                self.hplaybar.Visible = 'on';
            end

            redraw(self);
        end
    end



    methods (Access = 'protected')

        function resize(self)
            resize@DataViewer(self);
            if ~isempty(self.api.resizeFcn)
                try
                    self.api.resizeFcn();
                catch
                    self.api.resizeFcn = [];
                end
            end
        end

        function playback(self)
            if ~isempty(self.api.playbackFcn)
                self.api.playbackFcn();
            end
        end

        function reset(self)
            self.hlisten_playbar.Enabled = false;
            stop(self.hplaybar);
            self.hplaybar.Min = 1;
            self.hplaybar.Max = 0;
            self.hplaybar.Visible = 'off';

            if ~isempty(self.api.deleteFcn)
                self.api.deleteFcn();
            end

            self.api = self.template;

            self.exportaxes = false;
            self.exportrect = [];
            self.isAllowExportImage = false;
            self.isAllowExportVideo = false;
        end
    end
end

