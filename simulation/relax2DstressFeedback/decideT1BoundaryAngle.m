function t1idx = decideT1BoundaryAngle(g, angMin)
% DECIDET1BOUNDARYANGLE Determine if any edges should undergo a T1
% transition. Boundary edges undergo a T1 transition if their external
% angle is less than a user-supplied minimum length.
%
%   2       1             3
%    \     /              |
%  c2 \   / c1  ===>      |    
%      \ /                |
%       0                 0

if isempty(angMin)
    t1idx = [];
    return;
end

%--------------------------------------------------------------------------
% Generate the CCW Oriented Polygon Defining the Tissue Boundary
%--------------------------------------------------------------------------
% NOTE: This method assumes that all cells are already CCW oriented
% internally and that no boundary antibonds are maintained within the cell
% lattice graph. It also assumes that the boundary polygon is simple and
% non-intersecting

% bdyCell = struct();

% The IDs of the boundary bonds in the global boundary list
assert(~any(g.bonds(:,3) == 0), 'Cell lattice contains boundary bonds');
bdyBondIDx = find(g.bonds(:,4) == 0);
% bdyBondIDx = find( (g.bonds(:,3) ~= 0) & (g.bonds(:,4) == 0) );

% Local IDs for the boundary bonds
% locBondIDx = (1:numel(bdyBondIDx)).';

% The (unsorted) edge connectivity list
bonds = g.bonds(bdyBondIDx, 1:2);

% The cells attached to each bond
% bondCells = g.bonds(bdyBondIDx, 3:4);

% The vertices defining the boundary edges
bdyVertIDx = unique(bonds(:));
bdyVerts = g.verts(bdyVertIDx, :);

% Local IDs for the boundary vertices
locVertIDx = (1:numel(bdyVertIDx)).';

% Update bonds to be defined in terms of local IDs
bonds = changem(bonds, locVertIDx, bdyVertIDx);

% Construct a graph representation of the boundary polygon
A = sparse( bonds(:), [bonds(:,2); bonds(:,1)], ...
    1, numel(locVertIDx), numel(locVertIDx) );
% if any(~ismember(sum(A,2), [0, 2]))
if any(sum(A,2) ~= 2)
    error('Vertices belong to more than two edges')
end
G = graph(A);

% Sort the edge list
newBonds = dfsearch(G, 1, {'edgetonew'});
newBonds = [ newBonds; newBonds(end,2), newBonds(1) ];

[~, sortIDx] = ismember( sort(newBonds, 2), sort(bonds, 2), 'rows');
bonds = bonds(sortIDx, :);
bdyBondIDx = bdyBondIDx(sortIDx, :);
% bondCells = bondCells(sortIDx, :);
assert(isequal(sort(bonds, 2), sort(newBonds, 2)), 'Bad sort');

% Calculate the area of the boundary polygon to determine if the polygon is
% CCW or CW
tissueArea = cross( bdyVerts(bonds(:,1), :), bdyVerts(bonds(:,2), :), 2);
tissueArea = dot(tissueArea, repmat([0 0 1], size(tissueArea,1), 1), 2);
tissueArea = sum(tissueArea) / 2;
if (tissueArea < 0)
    warning('Bond order in cells does not appear to be CCW');
    bonds = bonds(:, [2 1]);
    % bondCells = bondCells(:, [2 1]);
end

if ~isequal(bonds(:,1), circshift(bonds(:,2), [1 0]))
    bonds = flipud(bonds);
    bdyBondIDx = flipud(bdyBondIDx);
    if ~isequal(bonds(:,1), circshift(bonds(:,2), [1 0]))
        error('Invalid boundary polygon construction');
    end
end

%--------------------------------------------------------------------------
% Determine Which Boundary Edges are Degenerate
%--------------------------------------------------------------------------

% Boundary edge vectors
edgeVecs = bdyVerts(bonds(:,2), :) - bdyVerts(bonds(:,1), :);

% Outward pointing boundary edge unit normals
edgeNormals = cross(edgeVecs, repmat([0 0 1], size(edgeVecs, 1), 1), 2);
edgeNormals = edgeNormals ./ sqrt(sum(edgeNormals.^2, 2));

% Calculate bending angles between subsequent edge unit normals
% Edges are parallel if bendAngles(i) == 0
% Vertex is convex if bendAngles(i) > 0
% Vertex is concave if bendAngles(i) < 0
en2 = circshift(edgeNormals, [-1 0]);
bendAngles = cross(edgeNormals, en2, 2);
bendAngles = bendAngles(:,3);
% bendAngles = dot(bendAngles, repmat([0 0 1], size(bendAngles, 1), 1), 2);
bendAngles = 2 * atan2( bendAngles, 1 + dot(edgeNormals, en2, 2));

allDegenerateEdges = (bendAngles < (-pi+angMin));
allDegenerateEdges = allDegenerateEdges | ((pi-angMin) < bendAngles);
if ~any(allDegenerateEdges)
    t1idx = [];
    return;
end

%--------------------------------------------------------------------------
% Determine Which Edges Have a Valid Valence Structure
%--------------------------------------------------------------------------
t1idx = []; %-ones(sum(allDegenerateAngles), 1);

wrapN = @(x, N) (1 + mod(x-1, N)); % Circular indexing function
% mergeID1 = find(allDegenerateEdges);
% mergeID2 = wrapN(mergeID1+1, size(edgeNormals, 1));
mergeID1 = find(bendAngles < (-pi+angMin));
mergeID1 = [mergeID1; ...
    wrapN(find((pi-angMin) < bendAngles), size(edgeNormals, 1)) ];
mergeID2 = wrapN(mergeID1+1, size(edgeNormals, 1));

while ~isempty(mergeID1)

    % Pop the next pair to merge from the queue
    e1i = mergeID1(1);
    e2i = mergeID2(1);
    mergeID1(1) = [];
    mergeID2(1) = [];

    % Avoid problems that may arise from merging subsequent pairs of edges
    rmIDx = mergeID1 == e2i;
    mergeID1(rmIDx) = [];
    mergeID2(rmIDx) = [];

    % Local vertex IDs defining the bonds that will be merged
    v1i = bonds(e1i, 1);
    v0i = bonds(e1i, 2);
    v2i = bonds(e2i, 2);

    assert(v1i ~= v2i, 'v1 == v2');
    assert(v1i ~= v0i, 'v1 == v0');
    assert(v0i ~= v2i, 'v2 == v0');

    % Check vertex valence
    v1in = g.bonds(:,2) == bdyVertIDx(v1i);
    v1out = g.bonds(:,1) == bdyVertIDx(v1i);
    v2in = g.bonds(:,2) == bdyVertIDx(v2i);
    v2out = g.bonds(:,1) == bdyVertIDx(v2i);

    badValence = (sum(v1in) == 1) && (sum(v1out) == 1) && ...
        (sum(v2in) == 1) &&(sum(v2out) == 1);

    if ~badValence
        t1idx = [t1idx; bdyBondIDx(e1i)];
    end

end