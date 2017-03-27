function B = bwr(m)
    % BWR - Blue White Red Colormap
    %
    %   Returns an M x 3 matrix containing colormap data. the colors begin with
    %   blue, progress through white, and end with red. BWR by itself is the
    %   same length as the current figure's colormap. If no figure exists,
    %   MATLAB will create one.
    %
    %   See also JET, HSV, PINK, FLAG, COLORMAP, RGBPLOT
    %
    %   Last Modified: 01-09-2014
    %   Modified By: Suever (suever@gmail.com)

    % Code shamelessly adapted from MATLAB
    if nargin < 1
        m = size(get(gcf,'colormap'),1);
    end

    red     = interp1([1 2 3], [0 1 1], linspace(1,3,m), 'linear');
    green   = interp1([1 2 3], [0 1 0], linspace(1,3,m), 'linear');
    blue    = interp1([1 2 3], [1 1 0], linspace(1,3,m), 'linear');

    B = [red(:), green(:), blue(:)];
end
