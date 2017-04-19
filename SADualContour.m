classdef SADualContour < ROIType
    methods
        function self = SADualContour()
            self@ROIType('hsadual', 'Short Axis Dual Ventricle', 3, ...
                'sadual', rand(16,16,3), true, true);
        end

        function [pos, iscls, iscrv, iscrn] = drawContour(self, hax, varargin)
            [outer, left, right] = getSAFull(hax, varargin{:});

            pos     = {outer, left, right};
            iscls   = {self.Closed};
            iscrv   = {self.Curved};
            iscrn   = {false};
        end
    end

    methods (Static)
        function tf = mask(X, Y, C)
            [inep,onep] = inpolygon(X, Y, C{1}(:,1), C{1}(:,2));
            [inen1,onen1] = inpolygon(X, Y, C{2}(:,1), C{2}(:,2));
            [inen2,onen2] = inpolygon(X, Y, C{3}(:,1), C{3}(:,2));

            tf = (inep & ~inen1 & ~inen2) | onep | onen1 | onen2;
        end
    end
end
