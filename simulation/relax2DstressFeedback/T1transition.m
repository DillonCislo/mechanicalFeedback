function [gnew, newBondIDx, oldBondIDx] = T1transition(g, bi)
    % Perform a T1 transition on a cell lattice. I'm assuming that bonds in
    % polygonal faces are oriented CCW. NOTE: THE FOLLOWING PROCEDURE FOR
    % RE-CONSTRUCTING THE CONNECTIVTY OF THE ALTERED LATTICE ASSUMES ALL
    % VERTICES HAVE VALENCE == 3
    % 
    % gnew = T1transition(g, bi)
    %
    % bi:   bond index
    % newBondIDx:   The updated bond IDs in the event of bond deletion
    % oldBondIDx:   The original bond IDs in the event of bond deletion
    %
    %    0                   0               0       0
    %     \       c2        /                 \     /
    %      \               /                   \   /
    %       \     abi     /                      v1 
    %          --------->                        |
    %  c4   v2 <-------- v1   c3 ========>  c4   |   c3
    %       /      bi     \                      v2
    %      /               \                    /  \
    %     /        c1       \                  /    \
    %    0                   0                0      0
    %

    %----------------------------------------------------------------------
    % Determine Connectivity of T1 Candidate Bond
    %----------------------------------------------------------------------

    % Assign here to avoid bad assignment errors for atypical termination
    oldBondIDx = [];
    newBondIDx = [];

    % Assess Validity of c1 and c2 ----------------------------------------

    v1i = g.bonds(bi,1); % The ID of the first vertex in the bond
    v2i = g.bonds(bi,2); % The ID of the second vertex in the bond
    c1i = g.bonds(bi,3); % The ID of the cell of the left of the bond
    c2i = g.bonds(bi,4); % The ID of the cell to the right of the bond

    % Turn off ALL boundary T1s
    % if ~(c1i && c2i)
    %     warning('no T1s are allowed on the boundary, returning');
    %     gnew = g;
    %     return
    % end

    % Check for detached edges
    if ~(c1i || c2i)
        warning(['Input edge is not attached to any cell. ' ...
            'Graph is invalid']);
        gnew = g;
        return;
    end

    % Antibond index (CELLS CAN SHARE MORE THAN ONE EDGE!)
    abi = find(g.bonds(:,2) == g.bonds(bi,1) & ...
        g.bonds(:,1) == g.bonds(bi,2));

    assert((numel(abi) == 1) || (numel(abi) == 0), ...
        'Input bond has multiple antibonds');

    % Switch halfedges so that c1 is a bulk cell
    if (~c1i && c2i)

        if isempty(abi)
            warning(['Boundary halfedge does not have bulk partner. ' ...
                'Graph is invalid']);
            gnew = g;
            return;
        end

        tmp = v1i; v1i = v2i; v2i = tmp;
        tmp = c1i; c1i = c2i; c2i = tmp;
        tmp = bi; bi = abi; abi = tmp;

    end

    if ~c2i
        assert(isempty(abi), 'Boundary halfedge exists in bond list?');
    end
    
    % Assess Validity of c3 and c4 ----------------------------------------

    % The bond that terminates on v1 and is a part of c1
    bintov1 = g.bonds(:,2) == v1i & g.bonds(:,3) == c1i;
    c3i = g.bonds(bintov1, 4); 

    % The bond that begins on v2 and is a part of c1
    bintov2 = g.bonds(:,1) == v2i & g.bonds(:,3) == c1i;
    c4i = g.bonds(bintov2, 4);
    
    % How can these ever be empty?
    % If they were it would mean c1 was missing halfedges...
    if isempty(c3i), c3i = 0; end    
    if isempty(c4i), c4i = 0; end

    if ~(c3i || c4i)
        warning(['c3 and c4 are both empty. ' ...
            'Disallowed T1 transition would change tissue topology']);
        gnew = g;
        return;
    end

    allCI = [c1i c2i c3i c4i];
    validCells = numel(unique(allCI)) == (sum(allCI ~= 0)+any(allCI == 0));
    if ~validCells
        warning('Disallowed T1 transition would create degenerate edge');
        gnew = g;
        return;
    end
    
    % disp(['T1 on cells ' num2str([c1i c2i c3i c4i])]);

    %----------------------------------------------------------------------
    % Update Vertex Positionts
    %----------------------------------------------------------------------

    % Rotate bond vertices CW around the bond midpoint
    % See bond re-ordering description to see why it has to be CW and not
    % CCW assuming that the bonds in each cell are all ordered CCW
    v1 = g.verts(v1i, :);
    v2 = g.verts(v2i, :);
    mp = v1 + (v2-v1)/2;
    v1 = cross(v1-mp, [0 0 1], 2) + mp;
    v2 = cross(v2-mp, [0 0 1], 2) + mp;
    g.verts(v1i, :) = v1;
    g.verts(v2i, :) = v2;

    %----------------------------------------------------------------------
    % Update Bond/Face Connectivity Lists
    %----------------------------------------------------------------------
    % NOTE: I'm basically assuming here that no boundary halfedges are ever
    % maintained within the mesh. This is potentially problematic ...
    
    % Handle c1 -----------------------------------------------------------
    % For c1 all we have to do is delete the primary bond and update the
    % subsequent bond to now begin on v1

    % bond, previous and next within the cell
	cbi = g.cells{c1i} == bi;
    bin = g.cells{c1i}(circshift(cbi, [0 1]));
    bip = g.cells{c1i}(circshift(cbi, [0 -1]));

    % Delete bond v1-->v2 from c1
    g.cells{c1i}(g.cells{c1i} == bi) = [];
    
    % The subsequent bond now begins on v1
    g.bonds(bin,1) = g.bonds(bip,2);

    % Handle c4 -----------------------------------------------------------
    % c4 inherits a new bond from the primary cell. If c4 is not the
    % exterior, we update the bond that previously terminated on v2 to now
    % terminate on v1, we add the primary bond from v1-->v2 after that, and
    % then update the next bond now begin on v2. If c4 is the exterior, we
    % must remove the primary bond from the bond list entirely
    % NOTE: This handling does NOT depend at all on c2 and can be performed
    % prior to c2 handling

    if c4i

        % The bond that terminates on v2 and is a part of c4
        b4 = find(g.bonds(:,2) == v2i & g.bonds(:,3) == c4i);

        % The location of this bond in the c4 bond list
        cb4 = find(g.cells{c4i} == b4);

        % This bond now terminates on v1
        g.bonds(g.cells{c4i}(cb4),2) = v1i;

        % Update the subsequent bond to begin on v2 (I think it
        % already should... ) and then add the bond from v1-->v2 in
        % between them
        if cb4 < numel(g.cells{c4i})
            g.bonds(g.cells{c4i}(cb4+1),1) = v2i;
            g.cells{c4i} = ...
                [g.cells{c4i}(1:cb4) bi g.cells{c4i}(cb4+1:end)];
        else
            g.bonds(g.cells{c4i}(1),1) = v2i;
            g.cells{c4i} = [g.cells{c4i}(1:cb4) bi];
        end

        g.bonds(bi,3:4) = [c4i c3i];

    else

        g.bonds(bi,:) = [0 0 0 0];

    end

    % Handle c2 and c3 ----------------------------------------------------

    if c2i

        % Handle Interior c2 ----------------------------------------------
        % If c2 is not the exterior, we must delete the primary antibond
        % v2-->v1, and update the subsequent bond to now begin on v2

        % antibond, previous and next within the cell
        cabi = g.cells{c2i} == abi;
        abin = g.cells{c2i}(circshift(cabi, [0 1]));
        abip = g.cells{c2i}(circshift(cabi, [0 -1]));

        % Delete antibond v2-->v1 from c2
        g.cells{c2i}(g.cells{c2i} == abi) = [];

        % The subsequent bond now begins on v2
        g.bonds(abin,1) = g.bonds(abip,2);

        % Handle c3 -------------------------------------------------------
        % c3 inherits a new bond from c2. If c3 is not the exterior, we
        % update the bond that previously terminated on v1 to now terminate
        % on v2, we add the primary antibond v2-->v1 after that, and then
        % update the next bond to begin on v1. If c3 is the exterior, we
        % must remove the primary antibond from the bond list entirely

        if c3i

            % The bond that terminates on v1 and is a part of c3
            b3 = find(g.bonds(:,2) == v1i & g.bonds(:,3) == c3i);

            % The location of this bond in the c3 bond list
            cb3 = find(g.cells{c3i} == b3);
            
            % This bond now terminates on v2
            g.bonds(g.cells{c3i}(cb3), 2) = v2i;

            % Update the subsequent bond to begin on v1 (I think it already
            % should... ) and then add the antibond from v2-->v1 in between
            % them
            if cb3 < numel(g.cells{c3i})
                g.bonds(g.cells{c3i}(cb3+1),1) = v1i;
                g.cells{c3i} = ...
                    [g.cells{c3i}(1:cb3) abi g.cells{c3i}(cb3+1:end)];
            else
                g.bonds(g.cells{c3i}(1),1) = v1i;
                g.cells{c3i} = [g.cells{c3i}(1:cb3) abi];
            end

            g.bonds(abi,3:4) = [c3i c4i];

        else

            g.bonds(abi,:) = [0 0 0 0];

        end

    else

        % Handle c3 -------------------------------------------------------
        % If c3 is not exterior, it still needs a new bond, but it cannot
        % now inherit that new bond from c2. Instead we must make an
        % entirely new bond and append it to the bond list. We
        % update the bond that previously terminated on v1 to now terminate
        % on v2, insert the new antibond v2-->v1 after that, and then
        % update the next bond to begin on v1. If c3 is the exterior, no
        % handling is necessary

        if c3i

            abi = size(g.bonds,1)+1;
            newBond = [v2i, v1i, c3i, c4i];

            % The bond that terminates on v1 and is a part of c3
            b3 = find(g.bonds(:,2) == v1i & g.bonds(:,3) == c3i);

            % The location of this bond in the c3 bond list
            cb3 = find(g.cells{c3i} == b3);
            
            % This bond now terminates on v2
            g.bonds(g.cells{c3i}(cb3), 2) = v2i;

            % Update the subsequent bond to begin on v1 (I think it already
            % should... ) and then add the antibond from v2-->v1 in between
            % them
            if cb3 < numel(g.cells{c3i})
                g.bonds(g.cells{c3i}(cb3+1),1) = v1i;
                g.cells{c3i} = ...
                    [g.cells{c3i}(1:cb3) abi g.cells{c3i}(cb3+1:end)];
            else
                g.bonds(g.cells{c3i}(1),1) = v1i;
                g.cells{c3i} = [g.cells{c3i}(1:cb3) abi];
            end

            g.bonds = [g.bonds; newBond];

        end

    end

    % Handle Bond Deletions -----------------------------------------------

    % The IDs of the bonds to remove
    rmIDx = find(all(g.bonds == 0, 2));

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