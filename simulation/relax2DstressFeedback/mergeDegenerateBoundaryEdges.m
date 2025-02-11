function gnew = mergeDegenerateBoundaryEdges(g, angMin)
%MERGEDEGENERATEBOUNDARYEDGES Two successive boundary edges whose external
%angle falls below a user specified threshold are merged into a single
%edge, joining the cells containing the original two edges.
%
%   NOTE: To fit with the rest of the cell division pipeline this function
%   will only merge the angle with the smallest minimum angle below the
%   threshold. It can be called multiple times until all degenerate edges
%   are merged. This should avoid issues with edge re-labelling
%
%   2       1         2-------1 
%    \     /              c2 /
%  c2 \   / c1  ===>        /  c1  
%      \ /                 /
%       0                 0
%
%   INPUT PARAMETERS:
%
%       - g:        The original cell lattice graph
%       - angMin:   The minimum angle below which edges are merged
%
%   OUTPUT PARAMETERS:
%
%       - gnew:     The updated cell lattice graph

%--------------------------------------------------------------------------
% Generate the CCW Oriented Polygon Defining the Tissue Boundary
%--------------------------------------------------------------------------
% NOTE: This method assumes that all cells are already CCW oriented
% internally and that no boundary antibonds are maintained iwthin the cell
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
bondCells = g.bonds(bdyBondIDx, 3:4);

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
bondCells = bondCells(sortIDx, :);
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
    bondCells = flipud(bondCells);
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

allDegenerateEdges = bendAngles < (-pi+angMin);
if ~any(allDegenerateEdges)
    gnew = g;
    return;
end

%--------------------------------------------------------------------------
% Merge the Two Edges With the Mininum Bending Angle
%--------------------------------------------------------------------------

wrapN = @(x, N) (1 + mod(x-1, N)); % Circular indexing function
mergeID1 = find(allDegenerateEdges);
mergeID2 = wrapN(mergeID1+1, size(edgeNormals, 1));

% disp('Merging Edge Pairs: ')
% disp([mergeID1, mergeID2]);

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
    v0in = g.bonds(:,2) == bdyVertIDx(v0i);
    v0out = g.bonds(:,1) == bdyVertIDx(v0i);
    v1in = g.bonds(:,2) == bdyVertIDx(v1i);
    v1out = g.bonds(:,1) == bdyVertIDx(v1i);
    v2in = g.bonds(:,2) == bdyVertIDx(v2i);
    v2out = g.bonds(:,1) == bdyVertIDx(v2i);

    goodValence = (sum(v1in) == 1) && (sum(v1out) == 1) && ...
        (sum(v2in) == 1) && (sum(v2out) == 1);

    if ~goodValence
        warning(['Invalid vertex valence structure. '...
            'Edge merger should be addressed by T1 instead']);
        continue;
    else
        warning('Merger actually happening');
    end

    % Cells attached to the bonds that will be merged
    c1i = bondCells(e1i, 1);
    c2i = bondCells(e2i, 1);

    if ((c1i == c2i) && (sum(v0in) == 1) && (sum(v0out) == 1))

        % In this case, the merger is actually an edge collapse.
        % We must remove the middle vertex
        
        % Remove the collapsed vertex
        bdyVerts(v0i, :) = nan(1,3);
        g.verts(bdyVertIDx(v0i), :) = nan(1,3);

        % Update the edge that will remain
        g.bonds(bdyBondIDx(e1i), :) = ...
            [bdyVertIDx(v1i), bdyVertIDx(v2i), c1i, 0];

        % Remove the collapsed edges
        g.bonds(bdyBondIDx(e2i), :) = nan(1,4);

        % Update the cell containing the bonds
        rmIDx = g.cells{c1i} == bdyBondIDx(e2i);
        g.cells{c1i}(rmIDx) = [];

        % Ensure no remaining edge contains the deleted vertex
        assert(~any(any(g.bonds(:,1:2) == bdyVertIDx(v0i))), ...
            'Failed edge collapse connectivity update');

    else

        % e1 is becomes a bulk bond pointing v1-->v0
        g.bonds(bdyBondIDx(e1i), :) = ...
            [ bdyVertIDx(v1i), bdyVertIDx(v0i), c1i, c2i ];

        % e2 becomes a bulk antibond pointing v0-->v1
        g.bonds(bdyBondIDx(e2i), :) = ...
            [ bdyVertIDx(v0i), bdyVertIDx(v1i), c2i, c1i ];

        % We must also add a new boundary bond v1-->v2
        g.bonds = cat(1, g.bonds, ...
            [ bdyVertIDx(v1i), bdyVertIDx(v2i), c2i, 0 ]);

        % We must also update the cells containing the bonds
        % NOTE: Only c2 should change - we simply insert
        % the new edge after e2 in the edge ordering
        e2Loc = find(g.cells{c2i} == bdyBondIDx(e2i));
        if e2Loc == numel(g.cells{c2i})
            g.cells{c2i} = [g.cells{c2i}, size(g.bonds, 1)];
        else
            g.cells{c2i} = [ g.cells{c2i}(1:e2Loc), ...
                size(g.bonds,1), g.cells{c2i}((e2Loc+1):end) ];
        end

    end

end

% Handle Vertex/Edge Deletions --------------------------------------------

% The IDs of the vertices to remove
rmIDx = find(any(isnan(g.verts), 2));

if ~isempty(rmIDx)

    vertInBonds = unique(g.bonds(:, 1:2));
    vertInBonds(isnan(vertInBonds)) = [];
    assert( ~any(ismember(rmIDx, vertInBonds)), ...
        'Vertex marked for deletion is still part of an extant bond' );

    oldVertIDx = (1:size(g.verts,1)).';
    oldVertIDx(rmIDx) = [];

    newVertIDx = (1:(size(g.verts,1)-numel(rmIDx))).';

    g.bonds(:, 1:2) = changem(g.bonds(:, 1:2), newVertIDx, oldVertIDx);

    g.verts(rmIDx, :) = [];

end

% Handle Bond Deletions ---------------------------------------------------

% The IDs of the bonds to remove
rmIDx = find(any(isnan(g.bonds), 2));

if ~isempty(rmIDx)

    assert( ~any(ismember(rmIDx, unique([g.cells{:}]))), ...
        'Bond marked for deletion is still part of extant cell' );

    oldBondIDx = (1:size(g.bonds,1)).';
    oldBondIDx(rmIDx) = [];

    newBondIDx = (1:(size(g.bonds,1)-numel(rmIDx))).';

    g.cells = cellfun(@(x) changem(x, newBondIDx, oldBondIDx), ...
        g.cells, 'UniformOutput', false);

    g.bonds(rmIDx, :) = [];

end

gnew = g;

end



