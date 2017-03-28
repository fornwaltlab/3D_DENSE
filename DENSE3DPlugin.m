classdef DENSE3DPlugin < plugins.DENSEanalysisPlugin
    % DENSE3DPlugin - A DENSEanalysis plugin
    %
    %   A plugin for analyzing 3D DENSE data
    %
    % Copyright (c) 2016, Jonathan Suever

    properties
        hpanel
        Handles
        hdense
    end

    methods

        function self = DENSE3DPlugin(varargin)
            import plugins.dense3D_plugin.*
            self@plugins.DENSEanalysisPlugin(varargin{:});

            self.hdense = DENSE3D();

            self.Handles.datalistener = addlistener(self.hdense, ...
                'NewData', @(s,e)self.refresh());

            % Only load the GUI components if there's actually a GUI
            if ~isempty(self.hfig)
                self.initGUI();
            end
        end

        function h = uimenu(self, varargin)
            h = gobjects(1,1);
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

            label = {'3D DENSE', 'Analysis'};
            [self.Handles.TabIndex, self.hpanel] = self.addTab(label);

            [~, self.Handles.hpopup] = self.addPopup('3D Workspaces', 150);
            [~, self.Handles.hoptions] = self.addPopup('3D Options', 80);
            [~, self.Handles.hresults] = self.addPopup('3D Results', 350);

            self.Handles.hpanel = uipanel( ...
                'Parent',           self.hpanel, ...
                'BorderType',       'none', ...
                'BackgroundColor',  [0 0 0]);

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

            %delete(self.Handles.TabListener);

            if ishghandle(self.hfig)
                % Only worry about removing tabs etc. if the figure is not
                % going to be deleted anyhow
                if strcmpi(get(self.hfig, 'BeingDeleted'), 'off')
                    h = guidata(self.hfig);
                    h.hpopup.removeTab(self.Handles.hpopup);
                    h.hpopup.removeTab(self.Handles.hresults);
                    h.hpopup.removeTab(self.Handles.hoptions);
                end

                % Call superclass destructor
                delete@plugins.DENSEanalysisPlugin(self);
            end
        end

        function validate(varargin)
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
        end

        function run(varargin)
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
