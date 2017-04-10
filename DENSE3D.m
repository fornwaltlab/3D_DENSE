classdef DENSE3D < hgsetget

    % TODO: Add APEX label to the slice that was determined automatically
    % to be the apex
    %
    % TODO: Add BASE label to the slice that was determined automatically
    % to be the base

    properties
        Apex
        Parameterization
        EndocardialMesh
        EndocardialMeshCut
        EpicardialMesh
        EpicardialMeshCut
        Interpolants
        LocalCoordinates
        RadialParams
        RBFParameters   % Parameters to use when fitting RBFs to the data
        Strains
        Flip = false
    end

    properties (Dependent)
        Data
    end

    properties (Hidden)
        data
    end

    events
        NewData
    end

    methods

        function set.Data(self, value)
            % Make sure that the data is unique such that there are no two
            % elements with the same SeriesInstanceUID

            if isempty(value)
                self.data = [];
                reset(self);
                return
            end

            seqs = arrayfun(@(x)x.SequenceInfo(1), value);
            if numel(seqs) > numel(unique({seqs.SeriesInstanceUID}))
                error(sprintf('%s:InvalidData', mfilename), ...
                    'There is duplicate data present')
            end

            seqs = arrayfun(@(x)x.SequenceInfo(1), value);
            [~, sortind] = sort([seqs.SliceLocation]);

            self.data = value(sortind);

            reset(self);
        end

        function reset(self)
            % reset - Resets all of the computed parameters because there
            % is either new data added or the ordering of the base/apex has
            % changed

            self.EndocardialMesh = [];
            self.EndocardialMeshCut = [];
            self.EpicardialMesh = [];
            self.EpicardialMeshCut = [];
            self.Interpolants = [];
            self.LocalCoordinates = [];
            self.RadialParams = [];
            self.RBFParameters = [];
            self.Strains = [];
            self.Parameterization = [];
            self.Apex = [];

            notify(self, 'NewData')
        end

        function set.Flip(self, value)
            self.Flip = value;
            reset(self);
        end

        function res = get.Data(self)
            if self.Flip
                res = flip(self.data);
            else
                res = self.data;
            end
        end

        function index = apicalSlice(self)
            index = numel(self.Data);
        end

        function movie(self)

            points = self.Interpolants(1).Points;

            figure

            axis equal

            plot3(points(:,1), points(:,2), points(:,3), '.');

            hold on

            endo = self.EndocardialMeshCut;

            inner = patch('Faces', endo.faces, 'Vertices', endo.vertices, 'FaceColor', 'w', 'FaceAlpha', 0.5);

            epi = self.EpicardialMeshCut;
            outer = patch('Faces', epi.faces, 'Vertices', epi.vertices, 'FaceColor', 'w', 'FaceAlpha', 0.5);

            axis equal
            axis manual

            for k = 1:1000

                frame = mod(k, numel(self.Interpolants)) + 1;

                delta = self.Interpolants(frame).query(endo.vertices);
                newpoints = endo.vertices + delta;
                set(inner, 'Vertices', newpoints);

                delta = self.Interpolants(frame).query(epi.vertices);
                newpoints = epi.vertices + delta;
                set(outer, 'Vertices', newpoints);
                pause
            end
        end


        function preview(self, cdata)
            figure;
            patch( ...
                'Faces', self.EndocardialMeshCut.faces, ...
                'Vertices', self.EndocardialMeshCut.vertices, ...
                'FaceColor', 'interp', ...
                'FaceVertexCData', cdata);

            hold on

            inds = dsearchn(self.EndocardialMeshCut.vertices, ...
                            self.EpicardialMeshCut.vertices);

            patch( ...
                'Faces', self.EpicardialMeshCut.faces, ...
                'Vertices', self.EpicardialMeshCut.vertices, ...
                'FaceColor', 'interp', ...
                'FaceVertexCData', cdata(inds));

            for k = 1:numel(self.Data)
                roi = self.Data(k).ROIInfo.RestingContour3D;
                plot3(roi{1}(:,1), roi{1}(:,2), roi{1}(:,3));
                plot3(roi{2}(:,1), roi{2}(:,2), roi{2}(:,3));
            end

            plot3(self.Apex(1), self.Apex(2), self.Apex(3), 'r*');

            axis equal

            hold on
        end

        function output = regionalStrains(self, nSegments)
            % Group strains based upon their parameterization

            % nSegments is: 16, 17 (default)
            if ~exist('nSegments', 'var');
                nSegments = 17;
            end

            % Compute the strains if necessary
            if isempty(self.Strains)
                self.computeStrains()
            end

            lsegments = linspace(0-eps, 1+eps, 5);
            lsegments(end) = [];

            csegments = linspace(0-eps, 1+eps, 25);
            csegments(end) = [];

            % Shift these circumferential segments

            % Break it into 4 segments longitudinally
            L = bsxfun(@lt, self.Strains.Parameterization.Longitudinal, lsegments);
            C = bsxfun(@gt, self.Strains.Parameterization.Circumferential, csegments);

            Lseg = min(4 - cumsum(L, 2), [], 2);
            Cseg = max(cumsum(C, 2), [], 2);

            fields = fieldnames(self.Strains);

            indices = accumarray([Cseg(:), Lseg(:)], 1:size(L,1), [], @(x){x(:).'});

            % This is going to be your standard model
            segments = {
                [indices{21:24, 4}]
                [indices{1:4,   4}]
                [indices{5:8,   4}]
                [indices{9:12,  4}]
                [indices{13:16, 4}]
                [indices{17:20, 4}]

                % Second longitudinal segment
                [indices{21:24, 3}]
                [indices{1:4,   3}]
                [indices{5:8,   3}]
                [indices{9:12,  3}]
                [indices{13:16, 3}]
                [indices{17:20, 3}]

                % Third longitudinal segment
                [indices{[20:24 1], 2}]
                [indices{2:7,    2}]
                [indices{8:13,   2}]
                [indices{14:19,  2}]

                % Apex
                [indices{:,1}]
            };

            % Omit the apex if needed
            if ismember(nSegments, [16, 18])
                segments(end) = [];
            end

            pts = self.Strains.Locations;
            segs = zeros(size(pts(:,1)));

            for k = 1:numel(segments)
                segs(segments{k}) = k;
            end

            % Create output struct to save everything in
            output = struct();
            output.Segmentation = segments;

            for k = 1:numel(fields)
                value = self.Strains.(fields{k});

                if isstruct(value)
                    continue
                end

                % Compute the mean value of whatever value this is within
                % this particular segment
                func = @(x)mean(self.Strains.(fields{k})(x,:), 1);
                tmp = cellfun(func, segments, 'UniformOutput', 0);
                output.(fields{k}) = cat(1, tmp{:});
            end

            import plugins.dense3D_plugin.*
            RD = RegionalDyssynchrony(permute(self.Strains.p2, [1 3 2]));
            delays = RD.computeRegionalDelays(mean(output.p2,1));
            output.DelayTimes = cellfun(@(x)mean(delays(x)), segments);

            % Compute CURE, RURE, and LURE

            if ismember(nSegments, [18, 19])
                lastinds = 13:18;
            else
                lastinds = 13:16;
            end

            output.CURE = [CURE(output.CC(1:6,:).'), ...
                           CURE(output.CC(7:12,:).'), ...
                           CURE(output.CC(lastinds,:).')].';

            output.RURE = [CURE(output.RR(1:6,:).'), ...
                           CURE(output.RR(7:12,:).'), ...
                           CURE(output.RR(lastinds,:).')].';

            output.LURE = [CURE(output.LL(1:6,:).'), ...
                           CURE(output.LL(7:12,:).'), ...
                           CURE(output.LL(lastinds,:).')].';

            output.CLShearAngle = rad2deg(self.torsion(output));
        end

        function rDelays = regionalDelays(self)
            regional = self.regionalStrains();
            strain = self.Strains.p2;

            RD = RegionalDyssynchrony(permute(strain, [1 3 2]));
            delays = RD.computeRegionalDelays(mean(regional.p2,1));

            rDelays = cellfun(@(x)mean(delays(x)), regional.Segmentation) * 100;
        end

        function computeStrains(self)

            mesh = self.EndocardialMeshCut;
            points = self.RadialParams.Points;
            self.Strains = queryStrains(mesh, points, self.Apex, ...
                                        self.Interpolants);

            % Figure out the radial/circumferential/longitudinal positions
            inds = dsearchn(mesh.vertices, points);

            self.Strains.Locations = points;

            p.Circumferential = self.Parameterization.Circumferential(inds);
            p.Longitudinal = self.Parameterization.Longitudinal(inds);
            p.Radial = self.RadialParams.Laplacian;

            self.Strains.Parameterization = p;
        end

        function coordinates(self)

            endo = self.EndocardialMeshCut;

            % Compute the normal for each face
            N = normr(normals(endo.vertices, endo.faces));

            % Compute the centroid of each face
            centers = barycenter(endo.vertices, endo.faces);

            [R,C,L] = localcoordinates(N, centers, self.Apex);

            self.LocalCoordinates = struct( ...
                'Points', centers, ...
                'Radial', R, ...
                'Circumferential', C, ...
                'Longitudinal', L);
        end

        function self = DENSE3D(matfiles)

            self.RBFParameters = struct( ...
                'Normalized', false, ...
                'Type', 'linear', ...
                'Constant', 0);

            if exist('matfiles', 'var')
                self.addData(matfiles)
            end
        end

        function interpolate(self, varargin)
            self.Interpolants = displacementSplines(self.Data, varargin{:});
        end

        function uipath = addData(self, data, description)

            uipath = '';

            if ~exist('data', 'var') || (ischar(data) && isdir(data))

                if exist('data', 'var')
                    startdir = data;
                else
                    startdir = pwd;
                end

                startdir = fullfile(startdir, '*.mat');
                [fname, pname] = uigetfile(startdir, 'MultiSelect', 'on');

                if isequal(fname, 0) || isequal(pname, 0)
                    return
                end

                uipath = pname;

                data = fullfile(pname, fname);
            end

            fields = {'ImageInfo', 'AnalysisInfo', 'SequenceInfo', ...
                      'ROIInfo', 'DENSEInfo'};

            if iscell(data)
                if exist('description', 'var')
                    cellfun(@self.addData, data, description, 'UniformOutput', false);
                else
                    cellfun(@self.addData, data, 'UniformOutput', false);
                end
                return
            elseif ischar(data)
                if ~exist('description', 'var')
                    [~, description] = fileparts(data);
                end
                data = load(data, fields{:}, '-mat');
            end

            % Delete any unnecessary fields
            data = rmfield(data, setdiff(fieldnames(data), fields));
            data.Description = description;

            for k = 1:numel(data)
                % Convert contours to 3D
                data(k) = self.two2three(data(k));
            end

            self.Data = cat(2, self.Data, data(:).');
        end

        function peaks = peakStrains(self)

            ops.RR = @(x)max(x, [], 2);
            ops.CC = @(x)min(x, [], 2);
            ops.LL = @(x)min(x, [], 2);

            regional = self.regionalStrains();

            fields = fieldnames(ops);

            peaks = struct();

            for k = 1:numel(fields)
                op = ops.(fields{k});
                peaks.(fields{k}) = op(regional.(fields{k}));
            end
        end

        function [base, baseindex] = basalSlice(self)
            % basalSlice - Returns information about the basal slice

            baseindex = 1;
            base = self.Data(baseindex).SequenceInfo(1);
            return

            seqs = arrayfun(@(x)x.SequenceInfo(1), [self.Data]);
            parallelids = [seqs.parallelid];
            issa = parallelids == mode(parallelids);

            % Get just the short-axis SequenceInfo values
            short_axis = seqs(issa);

            dists = nan(size(short_axis));

            for k = 1:numel(dists)
                IOP = short_axis(k).ImageOrientationPatient;
                IPP = short_axis(k).ImagePositionPatient;

                normal = cross(IOP(1:3), IOP(4:6));

                dists(k) = point2planeDistance(self.Apex, IPP, normal);
            end

            % Find the point that was the absolute furthest from the base
            [~, baseindex] = max(abs(dists));

            base = short_axis(baseindex);
        end

        function radialSample(self)
            endo = self.EndocardialMeshCut;
            epi = self.EpicardialMeshCut;

            % Rotate these so that the normal is facing up
            desired = [0 0 1];

            base = self.basalSlice();

            normal = cross(base.ImageOrientationPatient(1:3), ...
                           base.ImageOrientationPatient(4:6));

            theta = acos(dot(desired, normal));
            vec = cross(desired, normal);

            R = aa2mat(vec, theta);

            center = mean(epi.vertices, 1);

            % Vert a list of ALL vertices
            demeaned = bsxfun(@minus, epi.vertices, center);
            verts = bsxfun(@plus, demeaned * R, center);
            epi.vertices = verts;

            demeaned = bsxfun(@minus, endo.vertices, center);
            verts = bsxfun(@plus, demeaned * R, center);

            endo.vertices = verts;

            limits = [min(epi.vertices, [], 1); max(epi.vertices, [], 1)];

            spacing = 2.0;

            % Expand the limits by 3 * spacing
            limits(1,:) = limits(1,:) - 3 * spacing;
            limits(2,:) = limits(2,:) + 3 * spacing;

            X = limits(1,1):spacing:limits(2,1);
            Y = limits(1,2):spacing:limits(2,2);
            Z = limits(1,3):spacing:limits(2,3);

            endomask = inpolyhedron(endo.faces, endo.vertices, X, Y, Z, 'FlipNormals', true);
            epimask = inpolyhedron(epi.faces, epi.vertices, X, Y, Z, 'FlipNormals', true);

            result = ones(size(epimask), 'single');
            result(epimask) = 0.5;
            result(endomask) = 0;

            laplacian = relaxationLaplacian3D(result, result == 0.5);

            [xx,yy,zz] = meshgrid(X, Y, Z);

            ismyo = epimask & ~endomask;

            X = xx(ismyo);
            Y = yy(ismyo);
            Z = zz(ismyo);
            L = laplacian(ismyo);

            % Convert these back to their normal coordinates
            XYZ = [X, Y, Z];

            XYZ = bsxfun(@plus, bsxfun(@minus, XYZ, center) * R.', center);

            self.RadialParams.Points = XYZ;
            self.RadialParams.Laplacian = L;
        end

        function pruneMesh(self)

            % Identify the basal slice
            base = self.basalSlice();

            % Now cut the mesh using this basal plane information
            normal = cross(base.ImageOrientationPatient(1:3), ...
                           base.ImageOrientationPatient(4:6));

            % Check to make sure that the normal points towards the apex
            apex = self.Data(self.apicalSlice).SequenceInfo(1);

            dist = point2planeDistance(apex.ImagePositionPatient.', ...
                base.ImagePositionPatient, normal);

            if dist < 0
                normal = -normal;
            end

            [verts, faces] = half_space_intersect( ...
                self.EndocardialMesh.vertices, ...
                self.EndocardialMesh.faces, ...
                base.ImagePositionPatient, ...
                normal, 'Cap', false);

            self.EndocardialMeshCut = struct('faces', faces, ...
                'vertices', verts);

            [verts, faces] = half_space_intersect( ...
                self.EpicardialMesh.vertices, ...
                self.EpicardialMesh.faces, ...
                base.ImagePositionPatient, ...
                normal, 'Cap', false);

            self.EpicardialMeshCut = struct('faces', faces, ...
                'vertices', verts);
        end

        function parameterize(self)

            msh = self.EndocardialMeshCut;

            long = longitudinalParameterization(msh.vertices, msh.faces, self.Apex);

            self.Parameterization.Longitudinal = long;

            % Circumferential parameterization

            [base, baseindex] = self.basalSlice();
            analysis = self.Data(baseindex).AnalysisInfo;

            % Take into account that PositionB is the right-most so we need
            % to shift our parameterization by 1/6 later on
            insertion = tformfwd(base.tform, [analysis.PositionB, 0]);

            % For each isoline find the point that is closest to the
            % insertion
            nLines = 50;

            isolines = linspace(0, 1, nLines + 2);
            isolines([1 end]) = [];

            points = zeros(0, 3);
            params = zeros(0, 1);

            % Find the point in the mesh that corresponds to the apex
            [~, apexind] = ismember(self.Apex, msh.vertices, 'rows');

            for k = 1:numel(isolines)
                % Draw an isoline at this particular value
                [isoline, tris, ~, ~, SU] = slice_isolines(msh.vertices, ...
                    msh.faces, long(:), isolines(k), 'Manifold', true);

                % Assumes that the base is 0 and the apex is 1
                % Do a quick check here to ensure that this is actually the
                % case
                assert(SU(apexind) < 0.5, ...
                    'Expected the apex to be 1 and base to be 0')

                meshkeep = SU <= (isolines(k) + eps);

                trikeep = all(ismember(tris, find(meshkeep)), 2);
                tris = tris(trikeep, :);

                B = ordered_outline(tris);

                isoline = isoline(B,:);

                insertion_distances = sum(bsxfun(@minus, insertion, isoline).^2, 2);

                [~, startind] = min(insertion_distances);

                isoline = circshift(isoline, [1-startind, 0]);

                % Now parameterize by arc length
                dists = sqrt(sum(diff(isoline([1:end 1],:), [], 1).^2, 2));
                arclength = [0; cumsum(dists(:))];

                arclength = arclength ./ arclength(end);
                arclength(end) = [];

                % Now check whether the orientation of this sampling is
                % correct.

                % Compute the normal vector to the plane based upon the
                % contour orientation
                tmp = normr(bsxfun(@minus, isoline, mean(isoline, 1)));
                mean_norm = mean(cross(tmp(1:end-1,:), tmp(2:end,:)),1);
                D = point2planeDistance(self.Apex, isoline(1,:), mean_norm);

                % If the apex was not on the "correct" side, then the
                % orientation of the contour was the opposite of what we
                % expected and it needs to be flipped
                if D < 0
                    arclength = 1 - arclength;
                end

                points = cat(1, points, isoline);
                params = cat(1, params, arclength(:));
            end

            % Shift the parameterization by 1/6 to account for PositionB
            % being the EDGE of segment 1 and not the LV insertion
            params = mod(params - (1/6), 1);

            % Try to get Circ locations for each vertex in the initial mesh
            inds = dsearchn(points, msh.vertices);
            self.Parameterization.Circumferential = params(inds);
        end

        function apex = autoApex(self)
            % Sort from base to apex
            seqs = arrayfun(@(x)x.SequenceInfo(1), [self.Data]);
            parallelids = [seqs.parallelid];
            issa = parallelids == mode(parallelids);

            % Get just the short-axis SequenceInfo values
            short_axis = seqs(issa);

            maxInd = zeros(size(short_axis));

            % Now compute the normal for each of these slices
            for k = 1:numel(short_axis)

                vertices = self.EndocardialMesh.vertices;

                IOP = short_axis(k).ImageOrientationPatient;
                IPP = short_axis(k).ImagePositionPatient;

                normal = cross(IOP(1:3), IOP(4:6));

                dists = point2planeDistance(vertices, IPP, normal);

                [~, maxInd(k)] = min(dists);
            end

            apexIndex = unique(maxInd);
            assert(numel(apexIndex) == 1, 'More than one apex point identified');

            apex = vertices(apexIndex,:);
        end

        function generateMeshes(self)
            % Generate the surface meshes for the endo- and epicardial
            % contours
            rois = [self.Data.ROIInfo];
            rois = cat(1, rois.RestingContour3D);

            self.EpicardialMesh = surfacemesh(cat(1, rois{:,1}));
            self.EndocardialMesh = surfacemesh(cat(1, rois{:,2}));

            self.Apex = self.autoApex();

            self.pruneMesh();

            % Chop these off using the most basal slice
        end

        function data = two2three(self, data)
            % Converts all 2D contours to 3D contours

            roi = data.ROIInfo;
            seq = data.SequenceInfo(1);

            % Pre-allocate these to be the same size
            roi.RestingContour3D = roi.RestingContour;
            roi.Contour3D = roi.Contour;

            for k = 1:numel(roi.RestingContour3D)
                C = roi.RestingContour{k};
                C(:,3) = 0;
                roi.RestingContour3D{k} = tformfwd(seq.tform, C);
            end

            for k = 1:numel(roi.Contour)
                C = roi.Contour{k};
                C(:,3) = 0;
                roi.Contour3D{k} = tformfwd(seq.tform, C);
            end

            % Update the ROIInfo field
            data.ROIInfo = roi;
        end
    end

    methods (Static)
        function result = torsion(strains)
            % torsion - Compute the CL shear angle from the strain tensor
            result = -asin(2 * strains.CL ./ ...
                sqrt((1 + 2 * strains.CC) .* (1 + 2 * strains.LL)));
        end
    end
end

function splines = displacementSplines(data, progressfunc)

    if ~exist('progressfunc', 'var')
        progressfunc = @(x,y)fprintf('Processing frame %d of %d\n', x, y);
    end

    import plugins.dense3D_plugin.*

    analysis = [data.AnalysisInfo];
    images = [data.ImageInfo];

    % Figure out which frames we need to use
    analysisFrames = cat(2, analysis.FramesForAnalysis);
    frames = max(analysisFrames(1,:)):min(analysisFrames(2,:));

    % Pre-allocate the RBFs
    RBFs = cell(size(frames));

    Xdata = cat(4, images.Xunwrap);
    Ydata = cat(4, images.Yunwrap);
    Zdata = cat(4, images.Zunwrap);

    % Now select only the frames that we care about (the 3rd dimension)
    Xdata = Xdata(:,:,frames,:);
    Ydata = Ydata(:,:,frames,:);
    Zdata = Zdata(:,:,frames,:);

    % Now permute the dimensions so that we have [rows, cols, part, frame]
    Xdata = permute(Xdata, [1 2 4 3]);
    Ydata = permute(Ydata, [1 2 4 3]);
    Zdata = permute(Zdata, [1 2 4 3]);

    locations = find(~isnan(Xdata) & ~isnan(Ydata) & ~isnan(Zdata));

    % Now we know the row, column, partition, and frame of all data
    [row, col, part, frame] = ind2sub(size(Xdata), locations);

    % Create a nPartions x nFrames cell array of indices that we can use to
    % grab the displacements and positions for each of these
    indices = accumarray([part, frame], 1:numel(locations), [], @(x){x});

    for frame = 1:size(indices, 2)
        % Update the progress
        progressfunc(frame, numel(frames));

        XYZ = cell(size(indices, 1), 1);

        for slice = 1:size(indices, 1)
            seq = data(slice).SequenceInfo(1);

            % Convert the rows/columns to XYZ
            inds = indices{slice, frame};
            rc = [col(inds), row(inds), zeros(size(inds))];



            % Convert the displacements from image coordinates to world
            % coordinates
            locs = locations(inds);

            iminfo = data(slice).ImageInfo;

            dxIM = Xdata(locs) .* iminfo.Multipliers(1) .* seq.PixelSpacing(2);
            dyIM = Ydata(locs) .* iminfo.Multipliers(2) .* seq.PixelSpacing(1);
            dzIM = Zdata(locs) .* iminfo.Multipliers(3) .* seq.PixelSpacing(1);

            % Negate the displacements as necessary
            dzIM = -dzIM;

            xax = seq.ImageOrientationPatient(1:3);
            yax = seq.ImageOrientationPatient(4:6);
            zax = cross(xax, yax);

            dx = dxIM * xax(1) + dyIM * yax(1) + dzIM * zax(1);
            dy = dxIM * xax(2) + dyIM * yax(2) + dzIM * zax(2);
            dz = dxIM * xax(3) + dyIM * yax(3) + dzIM * zax(3);

            Xdata(locs) = dx;
            Ydata(locs) = dy;
            Zdata(locs) = dz;

            % Subtract the displacements from the points

            xyz_world = tformfwd(seq.tform, rc);


            XYZ{slice} = xyz_world - [dx, dy, dz];
        end

        allinds = cat(1, indices{:,frame});

        xyz = cat(1, XYZ{:});

        % Make sure that we don't have any duplicate values
        if ~isequal(size(unique(xyz, 'rows')), size(xyz))
            error('Non-unique displacement data')
        end

        I = locations(allinds);

        disps = single([Xdata(I), Ydata(I), Zdata(I)]);
        xyz = single(xyz);

        RBFs{frame} = RBFInterpolator(xyz, disps, 'Type', RBFInterpolator.LINEAR);
        disp frame
    end

    splines = cat(1, RBFs{:});
end


function mask = relaxationLaplacian3D(mask, innervals)

    % Laplacian filter
    kernel = ones(3,3,3,'single');
    kernel = kernel / numel(kernel);

    % Initialize convergence criteria to infinity
    convergence = Inf;

    % Pick 100 random points to test for convergence
    pos = find(innervals);
    opt = randperm(numel(pos));
    ndx = pos(opt(1:min(100, numel(pos))));

    while convergence > 0.001
        tmp = convn(mask, kernel, 'same');

        % Compute convergence criteria
        convergence = max(abs(tmp(ndx) - mask(ndx))) / min(mask(ndx));

        % Now update the mask with the relaxed values
        mask(innervals) = tmp(innervals);
    end
end

function M = aa2mat(ax, theta)
    s = sin(theta);
    c = cos(theta);
    t = 1 - c;

    ax = normr(ax);

    x = ax(1);
    y = ax(2);
    z = ax(3);

    M = [t*x*x + c,     t*x*y - s*z,    t*x*z + s*y;
         t*x*y + s*z,   t*y*y + c,      t*y*z - s*x;
         t*x*z - s*y,   t*y*z + s*x,    t*z*z + c];
end
