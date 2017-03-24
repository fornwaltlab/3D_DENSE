function D = longitudinalParameterization(vertices, faces, apex)

    % Find the point that is closest to the apex to use as the apex index
    distances = sum(bsxfun(@minus, apex, vertices).^2, 2);
    [~, apexIndex] = min(distances);

    % Compute the distance between each point and the apex
    dist2apex = heat_geodesic(vertices, faces, apexIndex);

    % Find the open edge of the mesh to use as the base seed points
    basepoints = outline(faces);

    dist2base = heat_geodesic(vertices, faces, basepoints);

    % Normalize the distance between the base and apex
    D = dist2base ./ (dist2apex + dist2base);
end
