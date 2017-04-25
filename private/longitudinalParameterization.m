function D = longitudinalParameterization(vertices, faces, apex)

    % Find the point that is closest to the apex to use as the apex index
    distances = sum(bsxfun(@minus, apex, vertices).^2, 2);
    [~, apexIndex] = min(distances);

    % Find the open edge of the mesh to use as the base seed points
    basepoints = unique(outline(faces));

    % Setup the biharmonic encoding so that we can quickly compute the
    % distance between vertices on the mesh
    B = biharmonic_embedding(vertices, faces, 4, 2);

    % Compute the distance from all vertices to the apex
    dist2apex = sqrt(sum(bsxfun(@minus, B(apexIndex,:), B).^2,2));

    % Now find the minimum distance to a base point
    D = sum(bsxfun(@minus, B(basepoints,:).', permute(B, [2 3 1])).^2, 1);
    dist2base = min(D, [], 2);

    % Normalize the distance between the base and apex
    D = dist2apex(:) ./ (dist2apex(:) + dist2base(:));
end
