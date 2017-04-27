classdef SADualContour < ROIType
    methods
        function self = SADualContour()

            cdata = [
                1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                1 1 1 2 3 3 3 3 3 3 3 3 2 1 1 1
                1 2 3 3 1 1 1 1 1 1 1 1 3 3 2 1
                2 3 1 1 2 3 3 3 3 2 3 3 2 1 3 2
                3 1 2 3 3 1 1 3 2 3 1 1 3 2 1 3
                3 1 3 1 1 1 3 2 3 1 1 1 1 3 1 3
                3 1 3 1 1 1 3 2 3 1 1 1 1 3 1 3
                3 1 2 3 3 1 1 3 2 3 1 1 3 2 1 3
                2 3 1 1 2 3 3 3 3 2 3 3 2 1 3 2
                1 2 3 3 1 1 1 1 1 1 1 1 3 3 2 1
                1 1 1 2 3 3 3 3 3 3 3 3 2 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

            cmap = [NaN NaN NaN;
                    0.7 0.7 0.7;
                    0.0 0.0 0.0];

            cdata = ind2rgb(cdata, cmap);

            self@ROIType('hsadual', 'Short Axis Dual Ventricle', 3, ...
                'sadual', cdata, true, true, {'r', 'g', 'y'});
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
