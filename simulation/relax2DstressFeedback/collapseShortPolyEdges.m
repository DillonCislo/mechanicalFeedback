function [newE, newV] = collapseShortPolyEdges(E, V, Lmin)
%COLLAPSESHORTPOLYEDGES Collapse edges of a polygon that are shorter than a
%user specified threshhold. Edges are replaced by a single new vertex at
%the average location of the two old vertices that defined the collapsed
%edge. NOTE: This function will also simply sort the edge list if you set
%the threshold to 0
%
%   INPUT PARAMETERS:
%
%       E:          #Ex2 list of vertex IDs defining polygon edges
%       V:          #VxD list of vertex coordinates
%       Lmin:       The minimum edge length
%
%   OUTPUT PARAMETERS:
%
%       newE:       #(E')x2 list of new vertex IDs defining edges
%       newV:       #(V')xD list of new vertex coordinates
%
%   by Dillon Cislo 07/06/2020

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Validate input properties
validateattributes(V, {'numeric'}, {'2d', 'real', 'finite'} );
validateattributes(Lmin, {'numeric'}, ...
    {'scalar', 'nonnegative', 'finite', 'real'});
validateattributes(E, {'numeric'}, ...
    {'2d', 'integer', 'real', 'finite', 'positive', '<=', size(V,1)} );

numE = size(E,1); % The number of polygon edges
numV = size(V,1); % The number of polygon vertices

if (numE ~= numV)
    error('Number of edges does not equal number of vertices');
end

% Construct a graph representation of the input polygon -------------------

% A vertex adjacency matrix
A = sparse( E(:), [E(:,2), E(:,1)], 1, numV, numV );

if any(sum(A,2) ~= 2)
    error('Vertices belong to more than two edges');
end

G = graph(A);
G.Nodes.V = V;

%--------------------------------------------------------------------------
% COLLAPSE EDGES
%--------------------------------------------------------------------------

while true
    
    newE = G.Edges.EndNodes; % The current edge list
    newV = G.Nodes.V; % The current vertex list
    
    % Calculate the edge lengths
    L = newV(newE(:,2), :) - newV(newE(:,1), :);
    L = sqrt(sum(L.^2, 2));
    
    % Determine if any edges are shorter than the threshold
    shortE = L < Lmin;
    if ~any(shortE), break; end
    
    % The edge to remove (we work one edge at a time)
    rmEdgeID = find(shortE, 1);
    
    % The IDs of the nodes that will be removed
    rmNodeIDx = newE(rmEdgeID, :);
    
    % Average the positions of the vertices defining the edge
    addV = [ newV(rmNodeIDx(1), :); newV(rmNodeIDx(2), :) ];
    addV = mean( addV, 1 );
    
    % The IDs of the new vertex
    addVID = numnodes(G) + 1;
    
    % Add the new nodes to the graph
    nodeProps = table( addV, 'VariableNames', {'V'} );
    G = addnode(G, nodeProps);
   
    % Find the neighbors of the vertices that will be removed
    neighborIDx = [ neighbors(G, rmNodeIDx(1)); ...
        neighbors(G, rmNodeIDx(2)) ];
    neighborIDx(ismember(neighborIDx, rmNodeIDx)) = [];
    
    % Add connecting edges betwen the collapsed edge neighbors and the new
    % vertex
    G = addedge(G, [addVID; addVID], neighborIDx, [1; 1]);
    
    % Remove the old vertices
    G = rmnode(G, rmNodeIDx(:));
    
end

%--------------------------------------------------------------------------
% SORT FINAL EDGE LIST
%--------------------------------------------------------------------------

% The final vertex list
newV = G.Nodes.V;

% The final sorted edge list
newE = dfsearch(G, 1, {'edgetonew'});
newE = [ newE; newE(end,2), newE(1) ];

end