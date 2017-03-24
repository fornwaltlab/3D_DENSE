function mesh = Generate_Mesh(points, adjust, center)

    if ~exist('adjust', 'var')
        adjust = false;
    end

    if ~exist('center', 'var')
        center = mean(points, 1);
    end

    thisdir = fileparts(mfilename('fullpath'));

    % model directions filename
    MODEL_DIRECTIONS_FILENAME = fullfile(thisdir, 'model_directions_1000.mat');

    % mesh filename
    MESH_FILENAME = fullfile(thisdir, 'unit_sphere_triangularization_200_200_4000.mat');

    % surface smoothness parameter
    SURFACE_LAMBDA = .001; % .0001; % .00001 works well

    % load and store model directions
    model = load(MODEL_DIRECTIONS_FILENAME);
    surface_model_u = [ model.model_directions.x ;
                        model.model_directions.y ;
                        model.model_directions.z ];

    % load mesh file
    MESH_INDEX = 20; % 10; % between 1 and 20
    triang = load(MESH_FILENAME);
    unit_sphere_triangularization = triang.unit_sphere_triangularization(MESH_INDEX);



    % fit LV surface
    [a0, a] = Fit_Surface_To_Data(surface_model_u, ...
                                  points, ...
                                  SURFACE_LAMBDA, ...
                                  center(1), center(2), center(3));

    mesh = Surface_Mesh(a0, a, ...
                        surface_model_u, ...
                        unit_sphere_triangularization, ...
                        center(1), center(2), center(3));

    if adjust
        mesh = Adjust_Mesh(mesh, points);
    end
end


function [a0, a] = Fit_Surface_To_Data(u, data, lambda, XC, YC, ZC)

    v = [ data(:,1)-XC data(:,2)-YC data(:,3)-ZC ]';
    [a0, a] = Fit_Surface(u, v, lambda);

end % Fit_Surface_To_Data

function mesh_out = Surface_Mesh(a0, a, u, unit_sphere_triangularization, XC, YC, ZC)

    mesh_u = [ unit_sphere_triangularization.x ; unit_sphere_triangularization.y ; unit_sphere_triangularization.z ];
    s = Evaluate_Surface(a0, a, u, mesh_u);

    mesh_in.node_num = length(unit_sphere_triangularization.x);
    mesh_in.tri_num = size(unit_sphere_triangularization.tri,1);
    mesh_in.node_x = XC + s(1,:);
    mesh_in.node_y = YC + s(2,:);
    mesh_in.node_z = ZC + s(3,:);
    mesh_in.tri_n1 = unit_sphere_triangularization.tri(:,1)';
    mesh_in.tri_n2 = unit_sphere_triangularization.tri(:,2)';
    mesh_in.tri_n3 = unit_sphere_triangularization.tri(:,3)';

    % mesh_out
    mesh_out = mesh_in;

end % Surface_Mesh


%
% Function: Adjust_Mesh
% This function adjusts the mesh input mesh to better match the points.
%

function mesh = Adjust_Mesh(mesh, points)

    % method parameters
    initial_nor_R = 10; % [mm]
    final_nor_R = 5; % [mm]
    initial_tan_R = 5; % [mm]
    final_tan_R = 5; % [mm]
    initial_smoothness_weight = .9;
    final_smoothness_weight = .9;
    initial_points_weight = .1;
    final_points_weight = .2;
    iterations = 100;

    newnodes = Adjust_Mesh_Mex(mesh, points, iterations, ...
                               [initial_smoothness_weight, final_smoothness_weight], ...
                               [initial_points_weight, final_points_weight],...
                               [initial_nor_R, final_nor_R], ...
                               [initial_tan_R, final_tan_R]);

    mesh.node_x = newnodes(:,1)';
    mesh.node_y = newnodes(:,2)';
    mesh.node_z = newnodes(:,3)';
end % Adjust_Mesh


function nodes = Surface_Normals_at_Nodes(nodes, mesh)

    % process node by node
    for n = 1:mesh.node_num,

        % compute the center
        ind = nodes(n).neighbors(:,1)';
        mean_x = mean(mesh.node_x(ind));
        mean_y = mean(mesh.node_y(ind));
        mean_z = mean(mesh.node_z(ind));

        % node normal
        X = [ mesh.node_x(ind)'-mean_x mesh.node_y(ind)'-mean_y mesh.node_z(ind)'-mean_z ];

        % compute the eigenvalues of X'*X
        [V,D] = eig(X'*X);

        % the first eigenvalue should be the smallest - check
        m = min([ D(1,1) D(2,2) D(3,3) ]);
        if (m ~= D(1,1))
            error('Adjust_Mesh::Surface_Normals_at_Nodes: the first eignevalue is not the smallest, which is assumed to be the case!');
        end

        % node normal
        nodes(n).normal = V(:,1);

        % set the normal orientation
        n1 = mesh.tri_n1(nodes(n).tri_ind);
        n2 = mesh.tri_n2(nodes(n).tri_ind);
        n3 = mesh.tri_n3(nodes(n).tri_ind);
        e12 = [ mesh.node_x(n2)-mesh.node_x(n1) mesh.node_y(n2)-mesh.node_y(n1) mesh.node_z(n2)-mesh.node_z(n1) ];
        e13 = [ mesh.node_x(n3)-mesh.node_x(n1) mesh.node_y(n3)-mesh.node_y(n1) mesh.node_z(n3)-mesh.node_z(n1) ];
        if (dot(nodes(n).normal,cross(e12,e13)) > 0)
            nodes(n).normal = -nodes(n).normal;
        end

    end % n loop

end % Surface_Normals_at_Nodes

function mesh = Next(nodes, mesh, points, NOR_R, TAN_R, smoothness_weight, points_weight)

    % compute node normals
    nodes = Surface_Normals_at_Nodes(nodes, mesh);

    % node displacements from the smoothness
    nodes = Node_Displacements_from_Smoothness(nodes, mesh);

    % node displacements from points
    nodenew = Node_Displacements_from_Points_Mex(nodes, mesh, points', NOR_R, TAN_R);

    % update the mesh
    for n = 1:mesh.node_num,
        mesh.node_x(n) = mesh.node_x(n) + smoothness_weight*nodes(n).smooth_disp(1) + points_weight*nodes(n).point_disp(1);
        mesh.node_y(n) = mesh.node_y(n) + smoothness_weight*nodes(n).smooth_disp(2) + points_weight*nodes(n).point_disp(2);
        mesh.node_z(n) = mesh.node_z(n) + smoothness_weight*nodes(n).smooth_disp(3) + points_weight*nodes(n).point_disp(3);
    end

end % Next

function nodes = Node_Displacements_from_Smoothness(nodes, mesh)

    % loop over nodes
    for n = 1:mesh.node_num,

        % compute the center
        ind = nodes(n).neighbors(:,1)';
        ind2 = [ ind(2:end) ind(1) ];
        dx = mesh.node_x(ind2) - mesh.node_x(ind);
        dy = mesh.node_y(ind2) - mesh.node_y(ind);
        dz = mesh.node_z(ind2) - mesh.node_z(ind);
        d = sqrt( dx.^2 + dy.^2 + dz.^2 );
        mx = .5 * ( mesh.node_x(ind2) + mesh.node_x(ind) );
        my = .5 * ( mesh.node_y(ind2) + mesh.node_y(ind) );
        mz = .5 * ( mesh.node_z(ind2) + mesh.node_z(ind) );
        center_x = sum(mx.*d) / sum(d);
        center_y = sum(my.*d) / sum(d);
        center_z = sum(mz.*d) / sum(d);

        % compute the average normal distance to triangle sides
        d21x = mesh.node_x(ind2) - mesh.node_x(ind);
        d21y = mesh.node_y(ind2) - mesh.node_y(ind);
        d21z = mesh.node_z(ind2) - mesh.node_z(ind);
        d10x = mesh.node_x(ind) - center_x;
        d10y = mesh.node_y(ind) - center_y;
        d10z = mesh.node_z(ind) - center_z;
        d = sqrt( (d10x.^2+d10y.^2+d10z.^2).*(d21x.^2+d21y.^2+d21z.^2) - (d10x.*d21x+d10y.*d21y+d10z.*d21z).^2 ) ./ sqrt(d21x.^2+d21y.^2+d21z.^2);
        d_avg = mean(d);

        % node normal
        node_normal = nodes(n).normal;


        % compute the intersection of second triangle planes with the node
        % normal
        cp = [ 0 0 0 ]';
        l0 = [ center_x center_y center_z ]';
        M = size(nodes(n).neighbors,1);
        for m = 1:M,
            mn = m+1;
            if (m == M)
                mn = 1;
            end
            n1 = nodes(n).neighbors(m,1);
            n2 = nodes(n).neighbors(mn,1);
            n3 = nodes(n).neighbors(m,2);
            p1 = [ mesh.node_x(n1) mesh.node_y(n1) mesh.node_z(n1) ]';
            p2 = [ mesh.node_x(n2) mesh.node_y(n2) mesh.node_z(n2) ]';
            p3 = [ mesh.node_x(n3) mesh.node_y(n3) mesh.node_z(n3) ]';
            p = Line_Plane_Intersection(l0, node_normal, p1, p2, p3);
            cp = cp + p;
            % tmp
            % plot3([p1(1) p(1)], [p1(2) p(2)], [p1(3) p(3)], 'm');
            % plot3([p2(1) p(1)], [p2(2) p(2)], [p2(3) p(3)], 'm');
            % plot3([p3(1) p(1)], [p3(2) p(2)], [p3(3) p(3)], 'm');
        end % m loop
        cp = cp / M;

        % compute the new point
        H = norm(cp-l0);
        h = d_avg * tan(atan(H/d_avg)/3);
        np = l0 + h * (cp-l0)/H;

        % store the displacement due to smoothness
        nodes(n).smooth_disp = [ np(1)-mesh.node_x(n) np(2)-mesh.node_y(n) np(3)-mesh.node_z(n) ];

    end % n loop

end % Node_Displacements_from_Smoothness


%
% Function: Line_Plane_Intersection
% -----------------------------------------------------------------------
% This function computes the intersection of a line (defined by point l0
% and direction l) and a plane (defined by three points p1, p2, p3). The
% function returns the point of intersection p. If there is no intersection
% p = [].
%
function p = Line_Plane_Intersection(l0, l, p1, p2, p3)

    % plane normal
    n = cross(p2-p1,p3-p1);
    nn = norm(n);
    if (nn <= eps)
        error('Adjust_Mesh::Line_Plane_Intersection: The three points do not define a plane!');
    end
    n = n / nn;

    num = dot(l,n);
    if (abs(num) <= eps)
        p = [];
    else
        d = dot((p1-l0),n) / num;
        p = l0 + d*l;
    end

end % Line_Plane_Intersection


%
% Function: Node_Neighborhoods
% -----------------------------------------------------------------------
% It computes the neighbood for each node. Each node will have a two-column
% matrix "neighbors". The first column contains the index of the
% neighboring node, while the second column contines the index of the "star
% tip second neighbor". Typically a node has six neighbors and therefore
% typically "neighbors" is a 6x2 matrix, representing neighbor info around
% the node.
%
function nodes = Node_Neighborhoods(mesh)

    % preallocate and initialize the nodes array
    nodes(mesh.node_num).count = [];
    for n = 1:mesh.node_num,
        nodes(n).count = 0;
        nodes(n).point_disp = [ 0 0 0 ];
    end

    % loop over nodes
    fprintf('Adjust_Mesh: Computing node neighbors: ');
    previous_percent = -1;
    for n = 1:mesh.node_num,


        % print progress
        percent = 10*floor(10*n/mesh.node_num);
        if (percent ~= previous_percent)
            if (percent == 10)
                fprintf('\b\b');
            end
            if (percent > 10)
                fprintf('\b\b\b');
            end
            fprintf('%d%%', percent);
            previous_percent = percent;
        end

        % find a triangle that contains node n
        ind = find(mesh.tri_n1 == n, 1);
        if isempty(ind)
            ind = find(mesh.tri_n2 == n, 1);
            if isempty(ind)
                ind = find(mesh.tri_n3 == n, 1);
                if isempty(ind)
                    error('Adjust_Mesh::Node_Neighborhoods: Invalid situation!');
                else
                    n1 = mesh.tri_n1(ind);
                    n2 = mesh.tri_n2(ind);
                end
            else
                n1 = mesh.tri_n3(ind);
                n2 = mesh.tri_n1(ind);
            end
        else
            n1 = mesh.tri_n2(ind);
            n2 = mesh.tri_n3(ind);
        end

        % store the triangle index
        nodes(n).tri_ind = ind;

        % find the forth point
        [n3, ~] = Find_Triangle(mesh, n1, n2, n);

        % store the neighbors
        nodes(n).neighbors = [ n1 n3 ];

        % go around the node
        ind = 1;
        while(ind)

            % find the next triangle
            [n4, ~] = Find_Triangle(mesh, n, n2, n1);

            % rename the points
            n1 = n2;
            n2 = n4;

            % find the forth point
            [n3, ~] = Find_Triangle(mesh, n1, n2, n);

            % store the neighbors
            nodes(n).neighbors = [ nodes(n).neighbors ; n1 n3 ];

            % check if the circle is complete
            if (n4 == nodes(n).neighbors(1,1))
                ind = 0;
            end

        end % while loop

    end % n loop
    fprintf('\n');

end % Node_Neighborhoods


%
% Function: Find_Triangle
% -----------------------------------------------------------------------
% This function finds the triangle, i.e. the node n3 that together with
% nodes n1 and n2 forms a triangle, but it does not contain node n.
%
function [n3, ind] = Find_Triangle(mesh, n1, n2, n)

    % n1 - n2 - n
    ind = find( (mesh.tri_n1 == n1) & (mesh.tri_n2 == n2) & (mesh.tri_n3 ~= n) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n1-n2-n)!');
        end
        n3 = mesh.tri_n3(ind);
        return;
    end

    % n1 - n - n2
    ind = find( (mesh.tri_n1 == n1) & (mesh.tri_n2 ~= n) & (mesh.tri_n3 == n2) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n1-n-n2)!');
        end
        n3 = mesh.tri_n2(ind);
        return;
    end

    % n2 - n1 - n
    ind = find( (mesh.tri_n1 == n2) & (mesh.tri_n2 == n1) & (mesh.tri_n3 ~= n) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n2-n1-n)!');
        end
        n3 = mesh.tri_n3(ind);
        return;
    end

    % n2 - n - n1
    ind = find( (mesh.tri_n1 == n2) & (mesh.tri_n2 ~= n) & (mesh.tri_n3 == n1) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n2-n-n1)!');
        end
        n3 = mesh.tri_n2(ind);
        return;
    end

    % n - n1 - n2
    ind = find( (mesh.tri_n1 ~= n) & (mesh.tri_n2 == n1) & (mesh.tri_n3 == n2) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n-n1-n2)!');
        end
        n3 = mesh.tri_n1(ind);
        return;
    end

    % n - n2 - n1
    ind = find( (mesh.tri_n1 ~= n) & (mesh.tri_n2 == n2) & (mesh.tri_n3 == n1) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n-n2-n1)!');
        end
        n3 = mesh.tri_n1(ind);
        return;
    end

    error('Adjust_Mesh::Find_Triangle: Invalid situation (no triangle match)!');

end % Find_Triangle

function nodes = Node_Displacements_from_Points(nodes, mesh, points, NORMAL_R, TANGENT_R)
    % nodes:    Nodes of the mesh
    % mesh:     Connectivity of nodes establishing 3D mesh
    % points:   User-provided points in 3D
    % NORMAL_R: Distance from the mesh in the normal direction
    % TANGENT_R:Distance from the mesh in the tangent direction

    orig = nodes;
    % process node by node
    for n = 1:mesh.node_num,
        % node normal

        node_normal = nodes(n).normal;

        % compute the normal and tangent components
        dx = points(:,1) - mesh.node_x(n);
        dy = points(:,2) - mesh.node_y(n);
        dz = points(:,3) - mesh.node_z(n);
        nor = dx*node_normal(1) + dy*node_normal(2) + dz*node_normal(3);
        tan = sqrt( (dx-nor*node_normal(1)).^2 + (dy-nor*node_normal(2)).^2 + (dz-nor*node_normal(3)).^2 );

        eds = (nor/NORMAL_R).^2 + (tan/TANGENT_R).^2;
        ind = find(eds <= 1);

        nodes(n).point_disp = [0 0 0];
        nodes(n).point_disp_nc = 0; % point displacement normal component
        c = 0;
        nodes(n).point_disp_computed = 0;
        if (isempty(ind) == 0)
            for m = 1:length(ind),
                im = ind(m);
                w = 1 - eds(im);
                nm = nor(im);
                nodes(n).point_disp = nodes(n).point_disp + w * nm * node_normal';
                nodes(n).point_disp_nc = nodes(n).point_disp_nc + w * nm;
                c = c + w;
            end % m loop
            nodes(n).point_disp = nodes(n).point_disp / c;
            nodes(n).point_disp_nc = nodes(n).point_disp_nc / c;
            nodes(n).point_disp_computed = 1;
        end

    end % n loop


    %
    % Spread the point displacements using "Laplacian smoothing"
    %

    % point_disp_normal
    point_disp_nc = zeros(1, mesh.node_num);
    for n = 1:mesh.node_num,
        if (nodes(n).point_disp_computed == 1)
            point_disp_nc(n) = nodes(n).point_disp_nc;
        end
    end

    M = 100; % number of iterations (practice shows that this is a reasonable value)
    for m = 1:M,

        % tmp storage
        nc = zeros(1, mesh.node_num);

        % process node by node
        for n = 1:mesh.node_num,

            % if the node point displacement needs to be adjusted
            if (nodes(n).point_disp_computed == 0)

                % node neighbors
                ind = nodes(n).neighbors(:,1)';

                % the average of neighbors
                nc(n) = mean(point_disp_nc(ind));

            end % if ...

        end % n loop

        % copy the results
        for n = 1:mesh.node_num,
            if (nodes(n).point_disp_computed == 0)
                point_disp_nc(n) = nc(n);
            end
        end % n loop

    end % m loop

    % get the final point displacements
    for n = 1:mesh.node_num,

        nodes(n).point_disp = point_disp_nc(n) * nodes(n).normal';

    end % n loop

end % Node_Displacements_from_Points
