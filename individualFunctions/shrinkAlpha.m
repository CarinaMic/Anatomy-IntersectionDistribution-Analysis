%% Shrink Wrap with alphaShape (without user-input)

% Input:    inter: struct to be updated (output container)
%           baseIntersect: label used for logging
%           typeName: 'alpha' | 'refinedAlpha' | 'inside' | 'grid'
%           intersectName: target field under inter.(typeName).*
%           aRadius: initial alpha radius (will be increased if open)
%           numDefect: number of defect vertices in the combined set (first block of indices)
%           numVertices: number of reference vertices in the combined set (follow the defect block)
%           inputPoints: base points used for shrink wrapping
%           addPoints: additional points appended to inputPoints (can be empty)

% Output:   inter: updated struct with inter.(typeName)

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [inter] = shrinkAlpha(inter, baseIntersect, typeName, intersectName, aRadius, numDefect, numVertices, ...
    inputPoints, addPoints)

% Combine input points and additional points
allVerticesPoints = [inputPoints; addPoints];

% Initialize
isClosed = false;
maxIterations = 150;
iterationCount = 0;
while ~isClosed && iterationCount < maxIterations
    try
        % Create alphaShape and get the boundary facets
        shrinkWrap = alphaShape(allVerticesPoints, aRadius);
        shrinkFaces = boundaryFacets(shrinkWrap);
        % Check if the mesh is closed
        edgesList = sort([shrinkFaces(:,1) shrinkFaces(:,2);
            shrinkFaces(:,2) shrinkFaces(:,3);
            shrinkFaces(:,3) shrinkFaces(:,1)], 2);
        [~, ~, edgeOccurrences] = unique(edgesList, 'rows', 'stable');
        edgeCounts = accumarray(edgeOccurrences, 1);
        % Check if any edge is not shared by exactly two triangles
        isClosed = all(edgeCounts == 2);
        % If the mesh is closed, break out of the loop (no increase of alpha radius)
        if isClosed
            break;
        end
    catch
        % If an error occurs, increase aRadius and try again
        isClosed = false;
    end

    % Increment alpha radius and iteration count
    aRadius = aRadius + 1;
    iterationCount = iterationCount + 1;
end

% Store results
inter.(typeName).(intersectName).radius = aRadius;
%obj.(typeName).(intersectName).alpha = shrinkWrap;
inter.(typeName).(intersectName).faces = shrinkFaces;
inter.(typeName).(intersectName).allVerticesPoints = allVerticesPoints; % all vertices
% Calculate used vertices and defect count
idxUsed = unique(shrinkFaces(:));
idxDefect    = idxUsed(idxUsed <= numDefect);
idxVertices  = idxUsed(idxUsed > numDefect & idxUsed <= numDefect+numVertices);
idxPoints    = idxUsed(idxUsed > numDefect+numVertices);
inter.(typeName).(intersectName).numOutputDefect   = numel(idxDefect);
inter.(typeName).(intersectName).numOutputVertices = numel(idxVertices);
inter.(typeName).(intersectName).numOutputPoints   = numel(idxPoints);
inter.(typeName).(intersectName).usedDefectCount   = numel(idxDefect) / numDefect;
inter.(typeName).(intersectName).usedCount = numel(idxUsed) / length(allVerticesPoints);

% Volume (divergence theorem)
volumeAlpha = 0; % Initialize volume
% Loop through each face
for i = 1:size(shrinkFaces, 1)
    % Get the vertices of the face
    v1 = allVerticesPoints(shrinkFaces(i, 1), :);
    v2 = allVerticesPoints(shrinkFaces(i, 2), :);
    v3 = allVerticesPoints(shrinkFaces(i, 3), :);
    % Calculate the signed volume of the tetrahedron
    tetraVolume = dot(v1, cross(v2, v3)) / 6;
    % Add to total volume
    volumeAlpha = volumeAlpha + tetraVolume;
end
% The volume should be positive
inter.(typeName).(intersectName).volume = abs(volumeAlpha);

disp(['alphaShape calculated: intersect type "', (typeName),'", pelvis defect categorie ', ...
    (baseIntersect), ' and ' ,(intersectName)])

end