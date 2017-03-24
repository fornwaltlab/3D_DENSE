function [normals, circum, long] = localcoordinates(normals, centroids, apex)
    % localcoordinates - Compute rotation matrix for local coordinates

    % Compute the longitudinal unit vectors
    vectorsLV = bsxfun(@minus, centroids, apex);

    % Compute dot product between these and unit normals
    dotp = dot(vectorsLV, normals, 2);

    % Calculate longitudinal vectors: Long = (V - (V . N) * N)
    long = normr(vectorsLV - bsxfun(@times, dotp, normals));

    % Now the circumferential is simply the cross product of the
    % longitudinal and radial unit vectors
    circum = cross(long, normals);

    % Now compute the rotation matrices which will be 3 x 3 x nFaces
    %rot(3,:,:) = long.';
    %rot(2,:,:) = circum.';
    %rot(1,:,:) = normals.';
end
