classdef DENSE3DPlugin < plugins.DENSEanalysisPlugin
    % DENSE3DPlugin - A DENSEanalysis plugin
    %
    %   A plugin for analyzing 3D DENSE data
    %
    % Copyright (c) 2016, Jonathan Suever

    properties
        hfig
        hpanel
        Handles
        hdense
    end

    methods

        function self = DENSE3DPlugin(varargin)
            import plugins.dense3D_plugin.*
            self@plugins.DENSEanalysisPlugin(varargin{:});

            self.hfig = findall(0, ...
                'Type', 'figure', ...
                'Name', 'DENSEanalysis', ...
                'tag', 'hfig');

            self.hdense = DENSE3D();

            self.Handles.datalistener = addlistener(self.hdense, ...
                'NewData', @(s,e)self.refresh());

            % Only load the GUI components if there's actually a GUI
            if ~isempty(self.hfig)
                self.initGUI();
            end
        end

        function refresh(self)
            hlist = self.Handles.hlist;
            newstrs = {self.hdense.Data.Description};
            val = min(numel(newstrs), get(hlist, 'value'));
            set(hlist, 'String', newstrs, 'Value', max([val,1]));
        end

        function initGUI(self)
            import plugins.dense3D_plugin.*
            h = guidata(self.hfig);

            self.hpanel = uipanel('Parent', self.hfig);
            h.hsidebar.addTab({'3D DENSE', 'Analysis'}, self.hpanel);

            self.Handles.TabIndex = h.hsidebar.NumberOfTabs;

            self.Handles.hpopup = uipanel( ...
                'Parent',   self.hpanel, ...
                'Position', [0 0 150 150]);

            self.Handles.hpanel = uipanel( ...
                'Parent',           self.hpanel, ...
                'BorderType',       'none', ...
                'BackgroundColor',  [0 0 0]);

            self.Handles.hpopup = uipanel( ...
                'Units',    'pixels', ...
                'Position', [0 0 200 200], ...
                'Parent',   self.hpanel);

            h.hpopup.addTab('3D Workspaces', self.Handles.hpopup);
            self.Handles.PopupIndex = h.hpopup.NumberOfTabs;

            h.hpopup.PanelHeights(self.Handles.PopupIndex) = 150;

            self.Handles.hlist = uicontrol( ...
                'Parent',   self.Handles.hpopup, ...
                'Style',    'Listbox', ...
                'String',   {}, ...
                'Units',    'normalized', ...
                'Position', [0.05 0.25 0.9 0.7]);

            self.Handles.himport = uicontrol( ...
                'Parent',   self.Handles.hpopup, ...
                'Style',    'pushbutton', ...
                'String',   'Add', ...
                'Units',    'Normalized', ...
                'Callback', @(s,e)self.hdense.addData(), ...
                'Position', [0.05 0.05 0.425 0.15]);

            self.Handles.hremove = uicontrol( ...
                'Parent',   self.Handles.hpopup, ...
                'Style',    'pushbutton', ...
                'String',   'Remove', ...
                'Units',    'Normalized', ...
                'Callback', @(s,e)self.removeCallback(), ...
                'Position', [0.525 0.05 0.425 0.15]);

            % Now add a popup for display options
            self.Handles.hresults = uipanel('Parent', self.hpanel);
            self.Handles.hoptions = uipanel('Parent', self.hpanel);

            h.hpopup.addTab('3D Options', self.Handles.hoptions);
            self.Handles.PopupIndex(end+1) = h.hpopup.NumberOfTabs;
            h.hpopup.PanelHeights(self.Handles.PopupIndex(end)) = 80;

            h.hpopup.addTab('3D Results', self.Handles.hresults);
            self.Handles.PopupIndex(end+1) = h.hpopup.NumberOfTabs;

            self.Handles.hviewer = DENSE3Dviewer( ...
                self.hdense, ...
                self.Handles.hpanel, ...
                self.Handles.hresults);

            bottoms = linspace(0.9, 0.05, 4);
            bottoms(1) = [];

            height = abs(diff(bottoms(1:2)));

            uicontrol( ...
                'Parent', self.Handles.hoptions, ...
                'Style', 'popupmenu', ...
                'String', {'Local Coords.', 'Cylindrical Coords.'}, ...
                'Units', 'normalized', ...
                'Enable', 'off', ...
                'Position', [0.05 bottoms(1) 0.9 height]);

            uicontrol( ...
                'Parent', self.Handles.hoptions, ...
                'Units', 'normalized', ...
                'Position', [0.05 bottoms(2) 0.9 height], ...
                'Style', 'checkbox', ...
                'String', 'Include Apex', ...
                'Enable', 'off', ...
                'Value', 1);

            uicontrol( ...
                'Parent', self.Handles.hoptions, ...
                'Units', 'normalized', ...
                'Position', [0.05 bottoms(3) 0.9 height], ...
                'Style', 'checkbox', ...
                'String', 'Flip Ventricle');

            h.hpopup.PanelHeights(self.Handles.PopupIndex(end)) = 350;

            [h.hpopup.Visible{self.Handles.PopupIndex}] = deal('off');

            % Add a listener to the tabs
            self.Handles.TabListener = addlistener(h.hsidebar, ...
                'SwitchTab', @(s,e)tab(s));

            function tab(src)
                if src.ActiveTab == self.Handles.TabIndex
                    self.Handles.state = h.hpopup.Visible;
                    h.hpopup.Visible = repmat({'off'}, size(h.hpopup.Visible));
                    [h.hpopup.Visible{self.Handles.PopupIndex}] = deal('on');
                else
                    if isfield(self.Handles, 'state')
                        h.hpopup.Visible = self.Handles.state;
                    end
                    [h.hpopup.Visible{self.Handles.PopupIndex}] = deal('off');
                end
            end
        end

        function removeCallback(self)
            val = get(self.Handles.hlist, 'Value');
            if val <= numel(self.hdense.Data)
                self.hdense.Data(val) = [];
            end
        end

        function delete(self)
            % Make sure all UI components are removed

            if ishghandle(self.hpanel)
                delete(self.hpanel)
            end

            delete(self.Handles.TabListener);

            if ishghandle(self.hfig)
                % Only worry about removing tabs etc. if the figure is not
                % going to be deleted anyhow
                if strcmpi(get(self.hfig, 'BeingDeleted'), 'off')
                    h = guidata(self.hfig);
                    h.hsidebar.removeTab(self.Handles.TabIndex);
                    h.hpopup.removeTab(self.Handles.hpopup);
                end

                % Call superclass destructor
                delete@plugins.DENSEanalysisPlugin(self);
            end
        end

        function validate(~, data)
            % validate - Check if the plugin can run.
            %
            %   Performs validation to ensure that the state of the program
            %   is correct to be able to run the plugin.
            %
            % USAGE:
            %   DENSE3DPlugin.validate(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all underlying
            %           data from the DENSEanalysis program.

            % Assert that image data base been loaded
            assert(~isempty(data.seq), ...
                'You must load imaging data into DENSEanalysis first.')

            % Make sure that we are dealing with 3D-encoded data
            is3D = ~cellfun(@(x)any(isnan(x)), {data.dns.PhaIndex});
            assert(any(is3D), ...
                '3D Analysis can only be run on 3D-encoded data.');
        end

        function run(self, data)
            % run - Method executed when user selects the plugin
            %
            % USAGE:
            %   DENSE3DPlugin.run(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all underlying
            %           data from the DENSEanalysis program.
        end
    end
end
