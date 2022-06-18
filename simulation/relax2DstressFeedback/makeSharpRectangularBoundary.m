function [vNew, cNew] = ...
    makeSharpRectangularBoundary(v, c, xLim, yLim, cellSize, ...
    collapseEdges, pickCornerPoints, customPoints)
%MAKESHARPRECTANGULARBOUNDARY Shifts vertices so that a polygonal tiling of
%a simply connected region of the plane occupies a rectangular domain
%
%   INPUT PARAMETERS:
%
%       - v:        #Vx2 list of vertex coordinates
%       - c:        #Cx1 cell containing a CCW ordered list of vertices
%                   definining polygonal faces
%       - xLim:     1x2 horizontal limits of the rectangular domain
%       - yLim:     1x2 vertical limits of the rectangular domain
%       - cellSize: Length scale of polygons
%       - collapseEdges: If true, collinear points within the same cells
%                        will be collapsed to share a single edge
%       - pickCornerPoints: If true, the user will be prompted to pick the
%                           four corner points
%       - customPoints:     Allow the use to supply the corner vertex IDs
%                           directly
%
%   OUTPUT PARMETERS:
%
%       - vNew:     #V'x2 list of new vertex coordinates
%       - cNew:     #C'x2 cell array of face connectivity lists
%
%   by Dillon Cislo 05/20/2022

% Input Processing --------------------------------------------------------
if (nargin < 5), collapseEdges = true; end
if (nargin < 6), pickCornerPoints = false; end
if (nargin < 7), customPoints = []; end

maxXVal = xLim(2);  minXVal = xLim(1);
assert(xLim(1) < xLim(2), 'Invalid horizontal domain limits');

maxYVal = yLim(2); minYVal = yLim(1);
assert(yLim(1) < yLim(2), 'Invalid vertical domain limits');

% Generate the CCW Oriented Polygon Defining the Tissue Boundary ----------
% NOTE: This method assumes that all cells are already CCW oriented
% internally and that no boundary antibonds are maintained iwthin the
% cell lattice graph. It also assumes that the boundary polygon is
% simple and non-intersecting

% The IDs of the vertices on the boundary polygon
bdyVertIDx = find(histcounts([c{:}], (0:size(v,1)).'+0.5).' < 3);
locVertIDx = (1:numel(bdyVertIDx)).';

% Extract the edges connecting boundary vertices
bdyEdges = cellfun(@(x) [x.', circshift(x.', [-1 0])], ...
    c, 'UniformOutput', false);
bdyEdges = cell2mat(bdyEdges);
bdyEdges = bdyEdges(all(ismember(bdyEdges, bdyVertIDx), 2), :);
bdyEdges = changem(bdyEdges, locVertIDx, bdyVertIDx);

% Construct a graph representation of the boundary polygon
A = sparse( bdyEdges(:), [bdyEdges(:,2); bdyEdges(:,1)], ...
    1, numel(locVertIDx), numel(locVertIDx) );
if any(sum(A,2) ~= 2)
    error('Vertices belong to more than two edges')
end
G = graph(A);

% Sort the edge list
newEdges = dfsearch(G, 1, {'edgetonew'});
newEdges = [ newEdges; newEdges(end,2), newEdges(1) ];
[~, sortIDx] = ismember( sort(newEdges, 2), sort(bdyEdges, 2), 'rows');
bdyEdges = bdyEdges(sortIDx, :);
assert(isequal(sort(bdyEdges, 2), sort(newEdges, 2)), 'Bad sort');

if (bdyEdges(1,2) ~= bdyEdges(2,1))
    bdyEdges = flipud(bdyEdges);
end

bdyEdges = changem(bdyEdges, bdyVertIDx, locVertIDx);
assert(size(bdyEdges,1) == numel(bdyVertIDx), ...
    ['Number of boundary vertices should ' ...
    'match number of boundary edges']);

% Divide the Polygonal Boundary Into Four Edges ---------------------------
%
%   c4 <--------------- c3
%   |        e3         ^
%   | e4             e2 |
%   |                   |
%  \/        e1         |
%   c1 ---------------> c2

if ~isempty(customPoints)

    v1i = customPoints(1);
    v2i = customPoints(2);
    v3i = customPoints(3);
    v4i = customPoints(4);

    shift = find(bdyEdges(:,1) == v1i)-1;
    bdyEdges = circshift(bdyEdges, [-shift 0]);

    e1IDx = false(size(bdyEdges, 1), 1);
    e1IDx(find(bdyEdges(:,1) == v1i):find(bdyEdges(:,2) == v2i)) = true;

    e2IDx = false(size(bdyEdges, 1), 1);
    e2IDx(find(bdyEdges(:,1) == v2i):find(bdyEdges(:,2) == v3i)) = true;
    
    e3IDx = false(size(bdyEdges, 1), 1);
    e3IDx(find(bdyEdges(:,1) == v3i):find(bdyEdges(:,2) == v4i)) = true;
    
    e4IDx = false(size(bdyEdges, 1), 1);
    e4IDx(find(bdyEdges(:,1) == v4i):end) = true;

    assert(sum([e1IDx; e2IDx; e3IDx; e4IDx]) == size(bdyEdges,1), ...
        'Bad edge orientation segregation');

elseif pickCornerPoints

    maxFaceSize = max(cellfun(@(x) numel(x), c));
    voronoiFace = nan(size(c,1), maxFaceSize);
    for i = 1:size(c,1)
        voronoiFace(i, 1:numel(c{i})) = c{i};
    end

    fig = figure('units', 'normalized', ...
    'outerposition', [0.5 0 0.5 1],  'color', [1 1 1]);

    patch('Faces', voronoiFace, ...
        'Vertices', v, 'FaceColor', 0.8 * ones(1,3), ...
        'EdgeColor', 'k', 'LineWidth', 2);

    hold on
    scatter(v(bdyVertIDx, 1), v(bdyVertIDx, 2), 'filled', 'c');
    hold off

    axis equal

    disp('Please select the bottom left corner');
    v1 = [0 0]; [v1(1), v1(2)] = ginput(1);

    disp('Please select the bottom right corner');
    v2 = [0 0]; [v2(1), v2(2)] = ginput(1);

    disp('Please select the top right corner');
    v3 = [0 0]; [v3(1), v3(2)] = ginput(1);

    disp('Please select the top left corner');
    v4 = [0 0]; [v4(1), v4(2)] = ginput(1);

    v1i = bdyVertIDx(knnsearch(v(bdyVertIDx, :), v1));
    v2i = bdyVertIDx(knnsearch(v(bdyVertIDx, :), v2));
    v3i = bdyVertIDx(knnsearch(v(bdyVertIDx, :), v3));
    v4i = bdyVertIDx(knnsearch(v(bdyVertIDx, :), v4));

    disp([v1i v2i v3i v4i]);

    shift = find(bdyEdges(:,1) == v1i)-1;
    bdyEdges = circshift(bdyEdges, [-shift 0]);

    e1IDx = false(size(bdyEdges, 1), 1);
    e1IDx(find(bdyEdges(:,1) == v1i):find(bdyEdges(:,2) == v2i)) = true;

    e2IDx = false(size(bdyEdges, 1), 1);
    e2IDx(find(bdyEdges(:,1) == v2i):find(bdyEdges(:,2) == v3i)) = true;
    
    e3IDx = false(size(bdyEdges, 1), 1);
    e3IDx(find(bdyEdges(:,1) == v3i):find(bdyEdges(:,2) == v4i)) = true;
    
    e4IDx = false(size(bdyEdges, 1), 1);
    e4IDx(find(bdyEdges(:,1) == v4i):end) = true;

    assert(sum([e1IDx; e2IDx; e3IDx; e4IDx]) == size(bdyEdges,1), ...
        'Bad edge orientation segregation');

    close(fig);

else

    % Segregate edges based on their orientation
    edgeAngles = atan2(v(bdyEdges(:,2), 2)-v(bdyEdges(:,1), 2), ...
        v(bdyEdges(:,2), 1)-v(bdyEdges(:,1), 1));

    e1IDx = -pi/4 <= edgeAngles & edgeAngles < pi/4;
    e2IDx = pi/4 <= edgeAngles & edgeAngles < 3*pi/4;
    e3IDx = 3*pi/4 <= edgeAngles | edgeAngles < -3*pi/4;
    e4IDx = -3*pi/4 <= edgeAngles & edgeAngles < -pi/4;

    assert(sum([e1IDx; e2IDx; e3IDx; e4IDx]) == numel(edgeAngles), ...
        'Bad edge orientation segregation');

    % WARNING: THIS MAY NOT BE ROBUST!!

    % Try to fill simple holes
    e1IDx(circshift(e1IDx, [-1 0]) & circshift(e1IDx, [1 0])) = true;
    eIDx = [e1IDx, e2IDx, e3IDx, e4IDx];

    % bwcc = bwconncomp(e1IDx);
    % [~, maxID] = max(cellfun(@numel, bwcc.PixelIdxList));
    % shift = bwcc.PixelIdxList{maxID}(1) - 1;

    if e1IDx(1)
        shift = find(~e1IDx, 1, 'last');
    else
        shift = find(e1IDx, 1, 'first')-1;
    end

    % eIDx0 = eIDx;
    eIDx = circshift(eIDx, [-shift 0]);

    tmpIDx= (1:find(eIDx(:,1), 1, 'last')).';
    eIDx(tmpIDx, :) = repmat([1 0 0 0], numel(tmpIDx), 1);

    tmpIDx = (find(eIDx(:,2), 1, 'first'):find(eIDx(:,2), 1, 'last')).';
    eIDx(tmpIDx, :) = repmat([0 1 0 0], numel(tmpIDx), 1);

    tmpIDx = (find(eIDx(:,3), 1, 'first'):find(eIDx(:,3), 1, 'last')).';
    eIDx(tmpIDx, :) = repmat([0 0 1 0], numel(tmpIDx), 1);

    tmpIDx = (find(eIDx(:,4), 1, 'first'):find(eIDx(:,4), 1, 'last')).';
    eIDx(tmpIDx, :) = repmat([0 0 0 4], numel(tmpIDx), 1);

    % eIDx = circshift(eIDx, [shift 0]);
    bdyEdges = circshift(bdyEdges, [-shift 0]);
    e1IDx = eIDx(:,1); e2IDx = eIDx(:,2);
    e3IDx = eIDx(:,3); e4IDx = eIDx(:,4);

end

% Pin Boundary Vertices to the Rectangular Boundary -----------------------

% Segregate vertices according to edge association
v1IDx = bdyEdges(e1IDx, :); v1IDx = [v1IDx(:,1); v1IDx(end, 2)];
v2IDx = bdyEdges(e2IDx, :); v2IDx = [v2IDx(:,1); v2IDx(end, 2)];
v3IDx = bdyEdges(e3IDx, :); v3IDx = [v3IDx(:,1); v3IDx(end, 2)];
v4IDx = bdyEdges(e4IDx, :); v4IDx = [v4IDx(:,1); v4IDx(end, 2)];

v(v1IDx, 2) = minYVal;
v(v2IDx, 1) = maxXVal;
v(v3IDx, 2) = maxYVal;
v(v4IDx, 1) = minXVal;

% Try to force interior vertices inside the rectangular boundary
v(v(:,1) < minXVal) = minXVal + 0.25 * cellSize;
v(v(:,1) > maxXVal) = maxXVal - 0.25 * cellSize;
v(v(:,2) < minYVal) = minYVal + 0.25 * cellSize;
v(v(:,2) > maxYVal) = maxYVal - 0.25 * cellSize;

% Jiggle boundary vertices to have the correct ordering
for i = 2:(numel(v1IDx)-1)
    if ( (v(v1IDx(i), 1) <= v(v1IDx(i-1), 1)) || ...
            (v(v1IDx(i+1), 1) <= v(v1IDx(i), 1)) )
        v(v1IDx(i), 1) = ( v(v1IDx(i+1), 1) + v(v1IDx(i-1), 1) ) / 2;
    end
end

for i = 2:(numel(v2IDx)-1)
    if ( (v(v2IDx(i), 2) <= v(v2IDx(i-1), 2)) || ...
            (v(v2IDx(i+1), 2) <= v(v2IDx(i), 2)) )
        v(v2IDx(i), 2) = ( v(v2IDx(i+1), 2) + v(v2IDx(i-1), 2) ) / 2;
    end
end

for i = 2:(numel(v3IDx)-1)
    if ( (v(v3IDx(i), 1) >= v(v3IDx(i-1), 1)) || ...
            (v(v3IDx(i+1), 1) >= v(v3IDx(i), 1)) )
        v(v3IDx(i), 1) = ( v(v3IDx(i+1), 1) + v(v3IDx(i-1), 1) ) / 2;
    end
end

for i = 2:(numel(v4IDx)-1)
    if ( (v(v4IDx(i), 2) >= v(v4IDx(i-1), 2)) || ...
            (v(v4IDx(i+1), 2) >= v(v4IDx(i), 2)) )
        v(v4IDx(i), 2) = ( v(v4IDx(i+1), 2) + v(v4IDx(i-1), 2) ) / 2;
    end
end

% Collapse Collinear Boundary Edges Within a Single Cell ------------------

if collapseEdges

    % Number of cells containing each vertex
    vcounts = histcounts([c{:}], (1:(size(v,1)+1))-0.5);
    allSingleV = find(vcounts == 1);

    for cid = 1:numel(c)

        singleV = ismember(c{cid}, allSingleV);
        if ~any(singleV), continue; end
        singleV = find(singleV);

        vid = 1:numel(c{cid});
        pid = circshift(vid, [0 1]);
        nid = circshift(vid, [0 -1]);

        badV = false(1, numel(c{cid}));
        for i = 1:numel(singleV)

            v1 = v(c{cid}(nid(singleV(i))), :) - ...
                v(c{cid}(vid(singleV(i))), :);
            v2 = v(c{cid}(pid(singleV(i))), :) - ...
                v(c{cid}(vid(singleV(i))), :);

            l1 = sqrt(sum(v1.^2, 2));
            l2 = sqrt(sum(v2.^2, 2));

            ang = 2*atan2( v1(1)*v2(2)-v1(2)*v2(1), l1*l2 + dot(v1, v2));
            if (abs(ang) < deg2rad(10)), badV(singleV(i)) = true; end

        end

        % THIS MAY NOT BE ROBUST
        c{cid}(badV) = [];

    end

end

incVIDx = unique([c{:}]); % The Voronoi vertices remaining
oldVIDx = (1:size(v,1)).'; % Original vertex IDs
rmVIDx = ~ismember(oldVIDx, incVIDx); % IDs of unused vertices
v(rmVIDx, :) = []; % Prune vertex list

% Prune individual cell vertex lists and renumber vertices
oldVIDx(rmVIDx) = [];
c = cellfun(@(x) changem(x, (1:size(v,1)).', oldVIDx), ...
    c, 'Uni', false);

assert(all(ismember((1:size(v,1)).', unique([c{:}]))), ...
    'Voronoi vertex inclusion mismatch');

vNew = v; cNew = c;

end