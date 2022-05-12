function gnew = removeDegenerateEdges(g)
%REMOVEDEGENERATEEDGES Remove degenerate edges that have collapsed within a
%single cell
%
%   INPUT PARAMETERS:
%
%       - g:    The original cell lattice graph
%
%   OUTPUT PARAMETERS:
%
%       - gnew: The updated cell lattice graph

% WARNING: This method may destroy the validity of the graph if multiple
% degenerate edges are removed at once. I have not checked it is robust.

% Determine which bonds/antibonds must be removed
badBond = false(size(g.bonds,1), 1);
for bi = 1:size(g.bonds, 1)

    % Skip if we already know the bond is bad
    if badBond(bi), continue; end

    % Unpack the current bond
    bond = g.bonds(bi,:);
    % v1i = bond(1); v2i = bond(2);
    c1i = bond(3); c2i = bond(4);

    % A degenerate edge must be contained entirely within a single cell
    if (c1i ~= c2i), continue; end
    
    % Find the antibond of the current bond
    abi = find((g.bonds(:,1) == bond(2)) & (g.bonds(:,2) == bond(1)));
    if isempty(abi), continue; end
    assert(numel(abi) == 1, ...
        'Bond has more than one antibond. Graph structure is invalid');
    antibond = g.bonds(abi, :);

    assert((antibond(3) == c2i) && (antibond(4) == c1i), ...
        'Bond and antibond do not agree on shared cells');

    % Check that no additional cells share the edge
    hasBadBond = cellfun(@(x) (any(x == bi) || any(x == abi)), g.cells);
    hasBadBond = find(hasBadBond);
    assert((numel(hasBadBond) == 1) && (hasBadBond == c1i), ...
        'Degenerate bond is somehow attached to other cells');

    % Mark the bond and antibond for removal
    badBond(bi) = true; badBond(abi) = true;
    g.bonds(bi, :) = nan(1,4);
    g.bonds(abi, :) = nan(1,4);

    % Remove the bonds from the cell list
    g.cells{c1i}(g.cells{c1i} == bi) = [];
    g.cells{c1i}(g.cells{c1i} == abi) = [];

end

if ~any(badBond)
    gnew = g;
    return;
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

% Handle Vertex Deletions -------------------------------------------------
% Remove vertices that are no longer referenced by any bond

% The IDs of the vertices to remove
rmIDx = find(~ismember((1:size(g.verts,1)).', unique(g.bonds(:, 1:2))));

if ~isempty(rmIDx)

    oldVertIDx = (1:size(g.verts,1)).';
    oldVertIDx(rmIDx) = [];

    newVertIDx = (1:(size(g.verts,1)-numel(rmIDx))).';

    g.bonds(:, 1:2) = changem(g.bonds(:, 1:2), newVertIDx, oldVertIDx);

    g.verts(rmIDx, :) = [];

else

    error('How can a vertex not be removed here?');

end

gnew = g;

end

    
    



