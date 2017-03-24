function msh = surfacemesh(points, varargin)

    m = Generate_Mesh(points, varargin{:});

    faces       = cat(2, m.tri_n1(:), m.tri_n2(:), m.tri_n3(:));
    vertices    = cat(2, m.node_x(:), m.node_y(:), m.node_z(:));

    verts = reshape(vertices(faces.', :), 3, [], 3);
    centroids = squeeze(mean(verts, 1));

    % Compute the surface normals using two edges
    e1 = squeeze(verts(2,:,:) - verts(1,:,:));
    e2 = squeeze(verts(3,:,:) - verts(1,:,:));
    normals = normr(cross(e1, e2));

    % Make sure that these normals point towards the centroid
    cdiff = bsxfun(@minus, centroids, mean(centroids, 1));

    direction = sign(dot(cdiff, normals, 2));

    % If the normals point away from the centroid (overall) then
    % flip them all
    if mode(direction) == 1
        normals = -normals;
    end

    msh = struct('faces',       faces, ...
                 'vertices',    vertices, ...
                 'normals',     normals, ...
                 'centroids',   centroids);
end
