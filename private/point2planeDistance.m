function dists = point2planeDistance(points, pt, normal)
    % point2planeDistance - Computes the min dist between a point and plane
    %
    % USAGE:
    %   distance = point2planeDistance(point, planePoint, planeNormal)
    %
    % INPUTS:
    %   point:          [M x 3] Matrix, Points to project onto the 3D plane.
    %                   The columns are the x, y, z coordinates
    %   planePoint:     [1 x 3] Array, Coordinates of a point that lies
    %                   within the plane of interest
    %   planeNormal:    [1 x 3] Array, Normal vector to the plane of
    %                   interest
    %
    % OUTPUTS;
    %   distance:       [M x 1] Array, Distances between each input point
    %                   and the defined plane
    %
    % Last Modified: 06-10-2014
    % Modified By: Jonathan Suever (suever@gmail.com)
    dists = bsxfun(@minus, points, pt(:)');
    dists = dists * (normal(:) ./ norm(normal));
end
