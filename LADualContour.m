classdef LADualContour < ROIType
    methods
        function self = LADualContour()

            cdata = [ ...
                1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1
                1 3 3 1 3 1 1 1 1 1 1 1 1 1 1 1
                1 3 1 2 3 1 1 3 1 1 1 1 1 1 1 1
                1 3 1 3 1 1 3 3 1 3 1 1 1 1 1 1
                1 3 1 3 1 1 3 1 1 3 1 1 1 1 1 1
                3 3 1 3 1 2 3 1 3 2 1 1 1 3 1 1
                3 1 2 3 1 3 1 1 3 1 1 1 2 3 1 3
                3 1 3 1 2 3 1 3 2 1 1 1 3 1 2 3
                3 1 3 1 3 1 1 3 1 1 1 2 3 1 3 1
                3 1 3 2 3 1 3 2 1 1 2 3 1 2 3 1
                3 1 3 3 1 2 3 1 1 2 3 1 1 3 1 1
                3 1 1 1 1 3 1 1 3 3 1 1 3 1 1 1
                2 3 1 1 1 3 2 3 3 1 1 3 1 1 1 1
                1 3 3 1 1 3 3 1 1 3 3 1 1 1 1 1
                1 2 3 1 1 1 1 2 3 2 1 1 1 1 1 1
                1 1 2 3 3 3 3 3 1 1 1 1 1 1 1 1];

            cmap = [NaN NaN NaN;
                    0.5 0.5 0.5;
                    0.0 0.0 0.0];

            cdata = ind2rgb(cdata, cmap);

            self@ROIType('hladual', 'Long Axis Dual Ventricle', 3, 'ladual', cdata, false, true);
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
