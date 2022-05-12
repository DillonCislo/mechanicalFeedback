function t1idx = decideT1(g, Lmin)
% DECIDET1 Determine if any edges should undergo a T1 transition. Edges
% undergo a T1 transition if they are shorter than a user-supplied minimum
% length. We arbitrarily choose only a single bond in the (bond, anti-bond)
% pair to be the T1 output ID

if (Lmin <= 0)
    t1idx = [];
    return;
end

% Determine duplicate bonds
bonds = sort(g.bonds(:, 1:2), 2);
dupl = find_duplicate_rows(bonds);

duplIDx = struct2cell(dupl);
duplIDx = duplIDx(2,:).';
assert(all(cellfun(@numel, duplIDx) == 2), 'Non-manifold edge found');
duplIDx = [duplIDx{:}].';

antiBondIDx = repmat((1:size(g.bonds,1)).', 1, 2);
antiBondIDx(duplIDx, 1) = [duplIDx(:,2); duplIDx(:,1)];

% Calculate all bond lengths
bondL = inf(size(bonds,1), size(g.verts, 2));
bulkBonds = (bonds(:,1) ~= 0) & (bonds(:,2) ~= 0);
bondL(bulkBonds, :) = ...
    g.verts(bonds(bulkBonds,2), :)-g.verts(bonds(bulkBonds,1), :);
bondL = sqrt(sum(bondL.^2, 2));

t1idx = (bondL < Lmin) & (antiBondIDx(:,1) <= antiBondIDx(:,2));
t1idx = find(t1idx);

end

function [ dupl, C ] = find_duplicate_rows( A )
%FIND_DUPLICATE_ROWS Finds duplicate rows in a matrix

[C, ia, ~] = unique( A, 'rows' );

if size(A,1) == size(C,1)
    % disp('There are no duplicate rows!');
    dupl = {};
    return;
end

rep_idx = setdiff(1:size(A,1), ia);
rep_val = unique( A(rep_idx,:), 'rows');

dupl_val = cell( size(rep_val,1), 1 );
dupl_idx = cell( size(rep_val,1), 1 );

for i = 1:size(rep_val,1)
    
    dupl_val{i} = rep_val(i,:);
    
    diffA = A - repmat( rep_val(i,:), size(A,1), 1 );
    dupl_idx{i} = find( sqrt(sum(diffA .* conj(diffA), 2)) < eps );
    
end

dupl = struct( 'val', dupl_val, 'idx', dupl_idx );

end
