classdef DENSE3D < hgsetget

    properties
        Apex
        Parameterization
        EndocardialMesh
        EpicardialMesh
        Interpolants
        LocalCoordinates
        RadialParams
        RBFParameters   % Parameters to use when fitting RBFs to the data
        Strains
        Flip = false
    end

    properties (Dependent)
        AnalysisFrames
        Data
    end

    properties (Hidden)
        data
    end

    events
        NewData
    end

    methods

        function res = get.AnalysisFrames(self)
            AI = [self.Data.AnalysisInfo];
            frames = arrayfun(@(x)x.FramesForAnalysis(:).', AI, 'Uniform', 0);
            frames = cat(1, frames{:});
            res = [max(frames(:,1)), min(frames(:,2))];
        end

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

        function bool = isBiventricular(self)
            % isBiventricular - Check if this is biventricular data or not

            bool = false;

            if isempty(self.Data)
                return
            end

            rois = [self.Data.ROIInfo];
            types = {'SAFull', 'sadual', 'LAFull', 'ladual'};
            bool = any(ismember(types, {rois.ROIType}));
        end

        function reset(self)
            % reset - Resets all of the computed parameters because there
            % is either new data added or the ordering of the base/apex has
            % changed

            self.EndocardialMesh = [];
            self.EpicardialMesh = [];
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

        function output = save(self, varargin)
            parser = inputParser();
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});

            opts = structobj(parser.Unmatched);

            %--- Image Info ---%
            if getfield(opts, 'ImageInfo', false)
                output.ImageInfo = [self.Data.ImageInfo];
            end

            %--- ROI Info ---%
            if getfield(opts, 'ROIInfo', true)
                output.ROIInfo = [self.Data.ROIInfo];
            end

            %--- Analysis Info ---%
            if getfield(opts, 'AnalysisInfo', true)
                AI.InterpolationMethod = 'RadialBasisFunction';
                AI.RBFConstants = [self.Interpolants.Constant];
                AI.RBFNormalized = self.Interpolants(1).Normalized;
                AI.RBFType = self.Interpolants(1).Type;
                AI.FramesForAnalysis = self.AnalysisFrames;

                % Eventually we want to make these options
                AI.CoordinateSystem = 'local';
                AI.Apex = true;

                output.AnalysisInfo = AI;
            end

            %--- DENSE Info ---%
            if getfield(opts, 'DENSEInfo', true)
                output.DENSEInfo = [self.Data.DENSEInfo];
            end

            %--- Sequence Info ---%
            if getfield(opts, 'SequenceInfo', true)
                output.SequenceInfo = cat(3, self.Data.SequenceInfo);
            end

            %--- Displacement Info ---%
            if getfield(opts, 'DisplacementInfo', true)

                points = self.samplePoints();
                nPoints = size(points, 1);

                DI.X = points(:,1);
                DI.Y = points(:,2);
                DI.Z = points(:,3);

                frames = self.AnalysisFrames;

                % Now sample the splines at these locations
                DI.dX = zeros(nPoints, numel(frames));
                DI.dY = zeros(nPoints, numel(frames));

                for k = 1:numel(frames)
                    D = single(self.Interpolants(k).query(points));
                    DI.dX(:,k) = D(:,1);
                    DI.dY(:,k) = D(:,2);
                    DI.dZ(:,k) = D(:,3);
                end

                output.DisplacementInfo = DI;
            end

            %--- Strain Info ---%
            if getfield(opts, 'StrainInfo', true)

                for endo = 1:numel(self.EndocardialMesh)

                    mesh = self.EndocardialMesh(endo);
                    strains = self.Strains(endo);

                    inds = dsearchn(mesh.vertices, strains.Locations);

                    SI.X = single(mesh.vertices(:,1));
                    SI.Y = single(mesh.vertices(:,2));
                    SI.Z = single(mesh.vertices(:,3));

                    tmp = rmfield(strains, {'Locations', 'Parameterization'});

                    fields = fieldnames(tmp);

                    nFrames = numel(self.Interpolants);

                    for k = 1:numel(fields)
                        values = accumarray(inds, 1:numel(inds), [], ...
                            @(x){mean(tmp.(fields{k})(x,:),1)}, {nan(1, nFrames+1)});
                        SI.(fields{k}) = single(cat(1, values{:}));
                    end

                    output.StrainInfo(endo) = SI;
                end
            end

            %--- Regional Strain Info ---%
            if getfield(opts, 'RegionalStrainInfo', true)
                RSI = self.regionalStrains();
                RSI = rmfield(RSI, {'Segmentation', 'Locations'});

                output.RegionalStrainInfo = RSI;
            end

            output.AnalysisInstanceUID = dicomuid;
        end

        function movie(self)
            points = self.Interpolants(1).Points;

            figure

            axis equal

            plot3(points(:,1), points(:,2), points(:,3), '.');

            hold on

            endo = self.EndocardialMesh;

            for m = 1:numel(endo)
                inner(m) = patch(endo(m), 'FaceColor', 'w', 'FaceAlpha', 0.5); %#ok
            end

            epi = self.EpicardialMesh;

            outer = patch(epi, 'FaceColor', 'w', 'FaceAlpha', 0.5);

            axis equal
            axis manual

            for k = 1:1000

                frame = mod(k, numel(self.Interpolants)) + 1;

                for m = 1:numel(endo)
                    delta = self.Interpolants(frame).query(endo(m).vertices);
                    newpoints = endo(m).vertices + delta;
                    set(inner(m), 'Vertices', newpoints);
                end

                delta = self.Interpolants(frame).query(epi.vertices);
                newpoints = epi.vertices + delta;
                set(outer, 'Vertices', newpoints);
                pause
            end
        end


        function preview(self, cdata)
            figure;

            if ischar(cdata)
                if isfield(self.Strains, cdata)
                    cdata = @(x)self.Strains(x).(cdata)(:,10);
                elseif isfield(self.Parameterization, cdata)
                    cdata = @(x)self.Parameterization(x).(cdata);
                end
            end


            for k = 1:numel(self.EndocardialMesh)
                patch(self.EndocardialMesh(k), ...
                    'FaceColor', 'interp', ...
                    'FaceVertexCData', cdata(k));

                hold on

                apex = num2cell(self.Apex(k,:));
                plot3(apex{:}, 'r*');
            end

            axis equal
        end

        function out = regionalStrains(self, nSegments)
            % Group strains based upon their parameterization

            % nSegments is: 16, 17 (default)
            if ~exist('nSegments', 'var');
                nSegments = [17, 19];
            end

            % Compute the strains if necessary
            if isempty(self.Strains)
                self.computeStrains()
            end

            out = [];

            for index = 1:numel(self.Strains)
                strains = self.Strains(index);

                lsegments = linspace(0-eps, 1+eps, 5);
                lsegments(end) = [];

                csegments = linspace(0-eps, 1+eps, 25);
                csegments(end) = [];

                % Break it into 4 segments longitudinally
                L = bsxfun(@lt, strains.Parameterization.Longitudinal, lsegments);
                C = bsxfun(@gt, strains.Parameterization.Circumferential, csegments);

                Lseg = min(4 - cumsum(L, 2), [], 2);
                Cseg = max(cumsum(C, 2), [], 2);

                fields = fieldnames(strains);

                indices = accumarray([Cseg(:), Lseg(:)], 1:size(L,1), [], @(x){x(:).'});

                switch nSegments(index)
                    case {16, 17}
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
                    case {18, 19}
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
                            [indices{21:24, 2}]
                            [indices{1:4,   2}]
                            [indices{5:8,   2}]
                            [indices{9:12,  2}]
                            [indices{13:16, 2}]
                            [indices{17:20, 2}]

                            % Apex
                            [indices{:,1}]
                        };
                    otherwise
                        error(sprintf('%s:InvalidSegments', mfilename), ...
                            'Segment number must be 16, 17, 18, or 19');
                end

                % Omit the apex if needed
                if ismember(nSegments(index), [16, 18])
                    segments(end) = [];
                end

                pts = strains.Locations;
                segs = zeros(size(pts(:,1)));

                for k = 1:numel(segments)
                    segs(segments{k}) = k;
                end

                % Create output struct to save everything in
                output = struct();
                output.Segmentation = segments;

                for k = 1:numel(fields)
                    value = strains.(fields{k});

                    if isstruct(value)
                        continue
                    end

                    % Compute the mean value of whatever value this is within
                    % this particular segment
                    func = @(x)mean(strains.(fields{k})(x,:), 1);
                    tmp = cellfun(func, segments, 'UniformOutput', 0);
                    output.(fields{k}) = cat(1, tmp{:});
                end

                if ismember(nSegments(index), [18, 19])
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

                % Use the mean LV contraction as a reference
                if index == 1
                    reference = mean(output.p2, 1);
                end

                import plugins.dense3D_plugin.*

                RD = RegionalDyssynchrony(permute(strains.p2, [1 3 2]));
                delays = RD.computeRegionalDelays(reference);
                output.DelayTimes = cellfun(@(x)mean(delays(x)), segments);

                out = cat(2, out, output);
            end
        end

        function computeStrains(self)

            if isempty(self.Parameterization)
                self.parameterize();
            end

            points = self.samplePoints();

            % Determine the "ownership" of each point.
            [~, epidist] = dsearchn(self.EpicardialMesh.vertices, points);

            nEndo = numel(self.EndocardialMesh);

            [indices, endodist] = deal(nan(size(points, 1), nEndo));

            for k = 1:nEndo
                [indices(:,k), endodist(:,k)] = dsearchn(self.EndocardialMesh(k).vertices, points);
            end

            if self.isBiventricular()
                % Only need this in the biventricular case
                isSeptum = all(bsxfun(@lt, endodist, epidist), 2);
            else
                isSeptum = false(size(epidist));
            end

            % Find the closest endo to each point
            [~, endo_ownership] = min(endodist, [], 2);

            % Now re-assign the septum such that it doesn't belong to any
            % particular ventricle
            endo_ownership(isSeptum) = 0;

            self.Strains = [];

            for k = 1:numel(self.EndocardialMesh)
                % Use the septal datapoints and the ones that belong to
                % this ventricle
                touse = isSeptum | endo_ownership == k;

                strains = queryStrains(self.EndocardialMesh(k), ...
                    points(touse,:), self.Apex(k,:), self.Interpolants);

                inds = indices(touse, k);

                param.Circumferential = self.Parameterization(k).Circumferential(inds);
                param.Longitudinal = self.Parameterization(k).Longitudinal(inds);

                strains.Parameterization = param;
                strains.Locations = points(touse,:);

                self.Strains = cat(2, self.Strains, strains);
            end
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

            frames = self.AnalysisFrames;

            self.Interpolants = displacementSplines(self.Data, ...
                frames(1):frames(end), self.Flip, varargin{:});
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

            peaks = repmat(struct(), size(regional));

            for k = 1:numel(fields)
                op = ops.(fields{k});
                for n = 1:numel(regional)
                    peaks(n).(fields{k}) = op(regional(n).(fields{k}));
                end
            end
        end

        function [base, baseindex] = basalSlice(self)
            % basalSlice - Returns information about the basal slice
            baseindex = 1;
            base = self.Data(baseindex).SequenceInfo(1);
        end

        function R = rotationMatrix(self)
            % Rotate these so that the normal is facing up
            desired = [0 0 1];

            base = self.basalSlice();

            normal = cross(base.ImageOrientationPatient(1:3), ...
                           base.ImageOrientationPatient(4:6));

            R = aa2mat(cross(desired, normal), ...
                       acos(dot(desired, normal)));
        end

        function [coordinates, epimask, endomask, myomask] = samplePoints(self, spacing)

            if isempty(self.EpicardialMesh)
                self.generateMeshes();
            end

            if ~exist('spacing', 'var')
                spacing = 2.0;
            end

            epi = self.EpicardialMesh;

            R = self.rotationMatrix();
            center = mean(epi.vertices, 1);

            epi = self.rotateMesh(epi, R, center);

            % Figure out the extent of the meshes
            limits = [min(epi.vertices, [], 1); max(epi.vertices, [], 1)];
            limits(1,:) = limits(1,:) - 3 * spacing;
            limits(2,:) = limits(2,:) + 3 * spacing;

            X = limits(1,1):spacing:limits(2,1);
            Y = limits(1,2):spacing:limits(2,2);
            Z = limits(1,3):spacing:limits(2,3);

            % Pre-allocate
            endomask = cell(size(self.EndocardialMesh));

            for k = 1:numel(self.EndocardialMesh)
                endo = self.rotateMesh(self.EndocardialMesh(k), R, center);
                endomask{k} = inpolyhedron(endo.faces, endo.vertices, ...
                    X, Y, Z, 'FlipNormals', true);
            end

            epimask = inpolyhedron(epi.faces, epi.vertices, X, Y, Z, 'FlipNormals', true);

            [xx,yy,zz] = meshgrid(X, Y, Z);

            % Do an OR across all endo masks
            myomask = epimask & ~any(cat(4, endomask{:}), 4);
            coordinates = [xx(myomask), yy(myomask), zz(myomask)];

            % Transform these coordinates back to their normal locations
            demeaned = bsxfun(@minus, coordinates, center);
            coordinates = bsxfun(@plus, demeaned * R.', center);
        end

        function radialSample(self)

            if self.isBiventricular()
                error(sprintf('%s:NotImplemented', mfilename), ...
                    'Transmural strains not yet implemented for biventricular DENSE')
            end

            [XYZ, epimask, endomask, myomask] = self.samplePoints();

            % Rotate through each of the endocardial masks and compute the
            % laplacian for each. Then we'll take the minimum of all of
            % them

            laplacian = nan(size(epimask));

            for k = 1:numel(endomask)
                result = ones(size(epimask), 'single');
                result(myomask) = 0.5;
                result(endomask{k}) = 0;

                laplacian = min(laplacian, relaxationLaplacian3D(result, myomask));
            end

            L = laplacian(myomask);

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

            empty = cell(size(self.EndocardialMesh));
            meshes = struct('faces', empty, 'vertices', empty);

            for k = 1:numel(self.EndocardialMesh)
                [verts, faces] = half_space_intersect( ...
                    self.EndocardialMesh(k).vertices, ...
                    self.EndocardialMesh(k).faces, ...
                    base.ImagePositionPatient, ...
                    normal, 'Cap', false);

                meshes(k).vertices = verts;
                meshes(k).faces = faces;
            end

            self.EndocardialMesh = meshes;

            [verts, faces] = half_space_intersect( ...
                self.EpicardialMesh.vertices, ...
                self.EpicardialMesh.faces, ...
                base.ImagePositionPatient, ...
                normal, 'Cap', false);

            self.EpicardialMesh = struct('faces', faces, ...
                'vertices', verts);
        end

        function parameterize(self)
            % Compute the parameterization for each endocardial mesh

            if isempty(self.EndocardialMesh)
                self.generateMeshes();
            end
            func = @self.parameterizeEndocardialMesh;
            nEndo = numel(self.EndocardialMesh);
            self.Parameterization = arrayfun(func, 1:nEndo);
        end

        function insertion = anteriorInsertion(self)
            % If this is single ventricle
            if self.isBiventricular()
                % THING to produce the anterior insertion LINE

                epioutline = ordered_outline(self.EpicardialMesh.faces);
                epioutline = self.EpicardialMesh.vertices(epioutline,:);

                % Find the EPI point that is equidistance from each contour
                outlines = cell(size(self.EndocardialMesh));

                nEndo = numel(self.EndocardialMesh);

                [epidist, epidistind] = deal(nan(size(epioutline, 1), nEndo));

                for k = 1:nEndo
                    msh = self.EndocardialMesh(k);
                    OOL = ordered_outline(msh.faces);
                    OOL = msh.vertices(OOL,:);

                    % Compute the point-wise distance between this and the
                    % EPI
                    distances = bsxfun(@minus, epioutline.', permute(OOL, [2 3 1]));
                    distances = squeeze(sum(distances.^2, 1));

                    [epidist(:,k), epidistind(:,k)] = min(distances, [], 2);

                    outlines{k} = OOL;
                end

                % Repeat the last point
                epidist = epidist([1:end 1],:);

                % Find where the endos cross
                sgn = sign(diff(epidist, [], 2));

                % Find the inferior and anterior insertion points
                insert  = find(sgn == 1, 1, 'first');
                indices = num2cell(epidistind(insert,:));
                tmp     = cellfun(@(x,y)x(y,:), outlines, indices, 'uni', 0);
                one     = mean(cat(1, epioutline(insert,:), tmp{:}));

                insert  = find(sgn == 1, 1, 'last');
                indices = num2cell(epidistind(insert,:));
                tmp     = cellfun(@(x,y)x(y,:), outlines, indices, 'uni', 0);
                two     = mean(cat(1, epioutline(insert,:), tmp{:}));

                % Now we need to figure out which one of these is the
                % anterior one
                centers = cellfun(@(x)mean(x, 1), outlines, 'uni', 0);

                N = cross(centers{2} - centers{1}, ...
                          self.Apex(1,:) - centers{1});

                insertions = [one; two];

                % Positive distance is going to be inferior
                side = sign(point2planeDistance([one; two], centers{1}, N));

                assert(isequal(sort(side(:)), [-1; 1]), ...
                    'Insertion point error');

                anterior = insertions(side == 1, :);
                % inferior = insertions(side == -1, :);

                % Just return the anterior for now
                insertion = anterior;
            else
                [base, baseindex] = self.basalSlice();
                analysis = self.Data(baseindex).AnalysisInfo;

                % Rotate the insertion point by pi/3
                theta = (2 * pi) / analysis.Nmodel;

                % Rotate the other way if we need it to go clockwise
                if analysis.Clockwise
                    theta = -theta;
                end

                % Create a rotation matrix to shift the anterior insertion
                % to the correct location
                R = [cos(theta), -sin(theta);
                     sin(theta), cos(theta)];

                center = analysis.PositionA;
                insertion2D = ((analysis.PositionB - center) * R) + center;
                insertion = tformfwd(base.tform, [insertion2D, 0]);
            end
        end

        function params = parameterizeEndocardialMesh(self, index)
            msh = self.EndocardialMesh(index);
            apex = self.Apex(index,:);

            long = longitudinalParameterization(msh.vertices, msh.faces, apex);

            params.Longitudinal = long;

            % Circumferential parameterization

            % Take into account that PositionB is the right-most so we need
            % to shift our parameterization by 1/6 later on
            insertion = self.anteriorInsertion();

            % For each isoline find the point that is closest to the
            % insertion
            nLines = 50;

            isolines = linspace(0, 1, nLines + 2);
            isolines([1 end]) = [];

            points = zeros(0, 3);
            cparam = zeros(0, 1);


            % Find the point in the mesh that corresponds to the apex
            [~, apexind] = ismember(apex, msh.vertices, 'rows');

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
                D = point2planeDistance(apex, isoline(1,:), mean_norm);

                % If the apex was not on the "correct" side, then the
                % orientation of the contour was the opposite of what we
                % expected and it needs to be flipped
                if D < 0
                    arclength = 1 - arclength;
                end

                points = cat(1, points, isoline);
                cparam = cat(1, cparam, arclength(:));
            end

            % Shift the parameterization by 1/6 to account for PositionB
            % being the EDGE of segment 1 and not the LV insertion
            % XXX This is now handled by anteriorInsertion itself
            %cparam = mod(cparam - (1/6), 1);

            % Try to get Circ locations for each vertex in the initial mesh
            inds = dsearchn(points, msh.vertices);
            params.Circumferential = cparam(inds);
        end

        function apex = autoApex(self)
            apx = self.Data(self.apicalSlice).SequenceInfo(1);

            % Figure out what side the basal slice is on
            base = self.basalSlice();

            normal = cross(apx.ImageOrientationPatient(1:3), ...
                           apx.ImageOrientationPatient(4:6));

            aIPP = apx.ImagePositionPatient;
            bIPP = base.ImagePositionPatient;

            D = point2planeDistance(bIPP.', aIPP, normal);

            % Determine the function for finding the apical point
            if D > 0
                func = @min;
            else
                func = @max;
            end

            apex = zeros(numel(self.EndocardialMesh), 3);

            for k = 1:numel(self.EndocardialMesh)
                vertices = self.EndocardialMesh(k).vertices;
                dists = point2planeDistance(vertices, aIPP, normal);

                % Take the mean of all vertices that are equally far
                apexind = dists == func(dists);
                apex(k,:) = mean(vertices(apexind,:), 1);
            end

            % apex is an M x 3 matrix where it contains the apex of each of
            % the M endocardial meshes
        end

        function generateMeshes(self)
            % Generate the surface meshes for the endo- and epicardial
            % contours
            rois = [self.Data.ROIInfo];
            rois = cat(1, rois.RestingContour3D);

            self.EpicardialMesh = surfacemesh(cat(1, rois{:,1}));
            self.EndocardialMesh = surfacemesh(cat(1, rois{:,2}));

            if isBiventricular(self)
                self.EndocardialMesh(2) = surfacemesh(cat(1, rois{:,3}), 1);
            end

            self.Apex = self.autoApex();

            % Remove the mesh vertices that are above the basal slice
            self.pruneMesh();
        end

        function data = two2three(~, data)
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

        function msh = rotateMesh(msh, R, center)
            for k = 1:numel(msh)
                demeaned = bsxfun(@minus, msh(k).vertices, center);
                msh(k).vertices = bsxfun(@plus, demeaned * R, center);
            end
        end
    end
end

function splines = displacementSplines(data, frames, flip, progressfunc)

    if ~exist('progressfunc', 'var')
        progressfunc = @(x,y)fprintf('Processing frame %d of %d\n', x, y);
    end

    import plugins.dense3D_plugin.*

    images = [data.ImageInfo];

    % Pre-allocate the RBFs
    RBFs = cell(size(frames));

    Xdata = cat(4, images.Xunwrap);
    Ydata = cat(4, images.Yunwrap);
    Zdata = cat(4, images.Zunwrap);

    if flip
        Zdata = -Zdata;
    end

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
    % Convert axis/angle notation to a rotation matrix
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
