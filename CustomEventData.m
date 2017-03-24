classdef CustomEventData < event.EventData & dynamicprops
    % CustomEventData - Object for creating events with custom fields
    %
    %   Typically, the user has to subclass event.EventData in order to
    %   create a custom data structure to return as an event. This class
    %   allows the user to create a data structure inline rather than
    %   creating a concrete class.
    %
    % USAGE:
    %   evnt = CustomEventData(S)
    %   evnt = CustomEventData(...)
    %
    % INPUTS:
    %   S:      Struct, Struct where each field contains data to be stored.
    %           The fields of the struct will be used as the names of the
    %           data fields in the event data.
    %
    %   ...:    Property/Value pairs, MATLAB-standard list of
    %           property(string) / value pairs

    methods
        function self = CustomEventData(varargin)
            ip = inputParser();
            ip.addParamValue('Source', [], @(x)true)
            ip.addParamValue('EventName', [], @(x)true)
            ip.KeepUnmatched = true;

            ip.parse(varargin{:})

            props = ip.Unmatched;

            fields = fieldnames(props);

            for k = 1:numel(fields)
                self.addprop(fields{k});
                self.(fields{k}) = props.(fields{k});
            end
        end
    end
end
