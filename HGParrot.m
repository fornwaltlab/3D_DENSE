classdef HGParrot < hgsetget & dynamicprops
    % HGParrot - Class for mimicking an internal graphics object
    %
    %   This class wraps an internal graphics object and makes it so that
    %   this object possesses the same properties and all set/get events
    %   are relayed to the underlying graphics object.
    %
    % USAGE:
    %   parrot = HGParrot(obj)
    %   parrot = HGParrot(classname)
    %
    % INPUTS:
    %   obj:        Handle, Handle to an existing graphics object that
    %               should be wrapped by this class.
    %
    %   classname:  String, Name of the internal graphics object that
    %               should be wrapped by this class.
    %
    % OUTPUTS:
    %   parrot:     Handle, Object which can be used to interact with the
    %               underlying graphics object.

    properties
        Type = 'HGParrot'   % Type which is displayed for the object
    end

    properties (Hidden)
        Handle      % Handle to the underlying graphics object
        Listener    % Delete listener to listen for object removal
    end

    methods
        function self = HGParrot(handle)
            % HGParrot - Constructor for the HGParrot object
            %
            % USAGE:
            %   parrot = HGParrot(obj)
            %   parrot = HGParrot(classname)
            %
            % INPUTS:
            %   obj:        Handle, Handle to an existing graphics object
            %               that should be wrapped by this class.
            %
            %   classname:  String, Name of the internal graphics object
            %               that should be wrapped by this class.
            %
            % OUTPUTS:
            %   parrot:     Handle, Object which can be used to interact
            %               with the underlying graphics object.

            if numel(handle) > 1 && ~ischar(handle)

                % Expand the array automatically to it's max size
                self(numel(handle)) = self;

                for k = 1:numel(handle)
                    self(k) = HGParrot(handle(k));
                end

                return
            end

            if ischar(handle)
                self.Handle = hgfeval(handle);
            else
                self.Handle = handle;
            end

            props = fieldnames(get(self.Handle));

            for k = 1:numel(props)
                % Skip properties that we have overloaded
                if isprop(self, props{k})
                    continue
                end

                prop = self.addprop(props{k});

                % Set/Get Methods
                prop.SetMethod = @(s,val)set(self.Handle, props{k}, val);
                prop.GetMethod = @(s,val)get(self.Handle, props{k});
            end

            % Link up a listener so that the object gets deleted properly
            self.Listener = addlistener(self.Handle, ...
                'ObjectBeingDestroyed', @(s,e)delete(self));
        end

        function delete(self)
            if ishghandle(self.Handle)
                delete(self.Handle)
            end
        end

        function res = double(self)
            % double - Converts this object to it's double equivalent
            %
            %   MATLAB graphics use a double handle typically rather than
            %   an object. This method returns that underlying double
            %   handle.
            %
            % USAGE:
            %   res = self.double()
            %
            % OUTPUTS:
            %   res:    Double, Handle to the underlying graphics object
            %           as a double.

            res = double(self.Handle);
        end

        function getdisp(self)
            disp(orderfields(get(self)));
        end

        function disp(self)
            % disp - Overloaded display of the object
            %
            %   The parrot object will simply be displayed similar to a
            %   graphics handle (as a double)
            %
            % USAGE:
            %   self.disp()

            if ~isvalid(self)
                disp@hgsetget(self)
            else
                disp(double([self.Handle]'));
            end
        end

        function res = ishghandle(self, varargin)
            % ishghandle - Overloaded function to behave like an hghandle
            %
            %   MATLAB uses this check periodically to determine if an
            %   object is an hghandle of a specific type. We relay these
            %   calls onto the underlying graphics object.
            %
            % USAGE:
            %   bool = self.ishghandle(type)
            %
            % INPUTS:
            %   type:   String, Type of handle that is requested. If no
            %           value is provided, then this tests if this is ANY
            %           type of graphics handle.
            %
            % OUTPUTS:
            %   bool:   Logical, Indicates if this is a graphics object of
            %           the desired type (true) or not (False)

            res = isvalid(self) && ishghandle(self.Handle, varargin{:});
        end

        function res = handle(self)
            % handle - Used to retrieve the underlying graphics handle
            %
            % USAGE:
            %   h = self.handle()
            %
            % OUTPUTS:
            %   h:  Handle, Handle to the underlying graphics object.

            res = self.Handle;
        end
    end

    methods (Access = 'protected')
        function addLinkedProperty(self, name, propname, hobj, validator)
            if ~exist('validator', 'var')
                validator = @(varargin)true;
            end

            if ~exist('propname', 'var')
                propname = name;
            end

            if ~exist('hobj', 'var')
                hobj = self.Handle;
            end

            prop = self.addprop(name);

            prop.SetMethod = @(s,val)setter(hobj, propname, val, validator);
            prop.GetMethod = @(s,val)get(hobj, propname);
        end
    end

    methods (Static, Hidden)
        function validateAxes(h, mfilename)
            types = {'axes', 'hggroup'};

            if ~ishghandle(h) || ~any(strcmpi(get(h, 'type'), types))
                error(sprintf('%s:InvalidParent', mfilename), ...
                    'Parent must be an axes or hggroup');
            end
        end

        function validateCallback(val, mfilename)
            msg = 'Callback must be a string, a function handle, or cell';
            if iscell(val) && ~isempty(val)
                if ~isa(val{1}, 'function_handle') || ~ischar(val{1})
                    error(sprintf('%s:InvalidCallback', mfilename), msg);
                end
            elseif ~isa(val, 'function_handle') && ~isempty(val) && ...
                   ~ischar(val)
                error(sprintf('%s:InvalidCallback', mfilename), msg);
            end
        end

        function validateCData(cdata, mfilename)
            if ~isnumeric(cdata) || ndims(cdata) > 3 || ...
                (ndims(cdata) == 3 && size(cdata, 3) ~= 3)
                error(sprintf('%s:InvalidCData', mfilename), ...
                    'CData must be a 2D or RGB matrix');
            end
        end
    end
end

function setter(hobj, propname, value, validator)
    if validator(value)
        set(hobj, propname, value)
    end
end

