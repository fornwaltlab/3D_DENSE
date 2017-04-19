classdef LADualContour < ROIType
    methods
        function self = LADualContour()
            self@ROIType('hladual', 'Long Axis Dual Ventricle', 3, 'ladual', rand(16,16,3), false, true);
        end

        function [pos, iscls, iscrv, iscrn] = drawContour(self, hax, varargin)
            [outer, left, right] = getLAFull(hax, varargin{:});

            pos     = {outer, left, right};
            iscls   = {self.Closed};
            iscrv   = {self.Curved};
            iscrn   = {false};
        end
    end

    methods (Static)
        function tf = mask(X, Y, C)
            C = cat(1, C{:});
            tf = inpolygon(X, Y, C(:,1), C(:,2));
        end
    end
end
