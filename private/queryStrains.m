function strains = queryStrains(mesh, query, apex, rbfs)
    % queryStrains - Queries an RBF at the specified location for strains
    %
    %   This function returns the cartesian, polar, and principal strains
    %   for any point supplied in the query argument.

    nFrames = numel(rbfs);
    nPoints = size(query, 1);

    assert(size(query, 2) == 3, 'Query points must be 3D coordinates.');

    N = normr(normals(mesh.vertices, mesh.faces));
    centers = barycenter(mesh.vertices, mesh.faces);

    [R, C, L] = localcoordinates(N, centers, apex);

    % Find the closest centers to each query point
    inds = dsearchn(centers, query);

    cfields = {'XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ'};
    pfields = {'RR', 'CR', 'LR', 'RC', 'CC', 'LC', 'RL', 'CL', 'LL'};
    ofields = {'p1', 'p2', 'p3'};
    allfields = [cfields, pfields, ofields];

    % Pre-allocate the output structure
    tmp = zeros(nPoints, nFrames + 1);
    values = cellfun(@(x)tmp, allfields, 'UniformOutput', false).';
    strains = cell2struct(values, allfields);

    rotation = [permute(R(inds,:), [3 2 1]);
                permute(C(inds,:), [3 2 1]);
                permute(L(inds,:), [3 2 1])];


    for frame = 2:(nFrames + 1)
        E = rbf2strain(rbfs(frame - 1), query);

        % Pre-allocate these
        Erot = E;

        D = nan(nPoints, 3);

        % Rotate the strain tensor
        for k = 1:size(E, 3)
            Erot(:,:,k) = rotation(:,:,k) * E(:,:,k) * rotation(:,:,k).';
            D(k,:) = eig(E(:,:,k), 'nobalance');
        end

        % Actually fill in the values in the output structure
        for k = 1:numel(cfields)
            [r,c] = ind2sub([3 3], k);
            strains.(cfields{k})(:,frame) = squeeze(E(r,c,:));
        end

        % Actually fill in the values in the output structure
        for k = 1:numel(pfields)
            [r,c] = ind2sub([3 3], k);
            strains.(pfields{k})(:,frame) = squeeze(Erot(r,c,:));
        end

        % Now the principal strains
        strains.p1(:,frame) = D(:,3);
        strains.p2(:,frame) = D(:,2);
        strains.p3(:,frame) = D(:,1);
    end
end

function [E, F] = rbf2strain(rbf, query)

    % Calculate the gradient tensor of the displacement field
    grad = rbfgradients(rbf, query);

    % Convert to the deformation gradient tensor
    F = bsxfun(@plus, grad, eye(3));

    % Pre-allocate E
    E = F;

    % Compute the Lagrangian Green strain tensor
    for k = 1:size(query, 1)
        E(:,:,k) = ((F(:,:,k).' *  F(:,:,k)) - eye(3)) / 2;
    end
end

function grad = rbfgradients(rbf, query)

    nPoints = size(query, 1);
    nNodes = size(rbf.Points, 1);

    % Find the distance between all points and the query points
    dx = bsxfun(@minus, query(:,1), rbf.Points(:,1).');
    dy = bsxfun(@minus, query(:,2), rbf.Points(:,2).');
    dz = bsxfun(@minus, query(:,3), rbf.Points(:,3).');

    % Compute the radius
    radius = sqrt(dx.^2 + dy.^2 + dz.^2);

    % Replace any zeros with episilon to avoid division by zero
    radius(radius == 0) = eps;

    ux = (dx ./ radius) * rbf.Weights(1:nNodes, :);
    uy = (dy ./ radius) * rbf.Weights(1:nNodes, :);
    uz = (dz ./ radius) * rbf.Weights(1:nNodes, :);

    grad = nan(3, 3, nPoints);

    grad(1,:,:) = ux.';
    grad(2,:,:) = uy.';
    grad(3,:,:) = uz.';

    grad = bsxfun(@plus, grad, rbf.Weights(nNodes+2:end,:));
end
