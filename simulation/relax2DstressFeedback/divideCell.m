function newg = divideCell(gcur, ci, divaxis, intersectMethod)
% gcur: input cell lattice graph structure
% ci: cell ID of the cell about to divide
% divaxis: 1x2 vector along division axis
%
% newg: output cell lattice graph structure including daughter cells

if (nargin < 4), intersectMethod = 'edgeAverage'; end
if isempty(intersectMethod), intersectMethod = 'edgeAverage'; end

% Update Graph Topology ---------------------------------------------------

offset = 0;
[faces, bonds, verts, nfi, nbi, ibi] = ...
    divideFace( gcur.cells, gcur.bonds, gcur.verts, ...
    ci, divaxis, offset, intersectMethod );

newg = gcur;
newg.cells = faces;
newg.bonds = bonds;
newg.verts = verts;

% Update Cell Area Parameters ---------------------------------------------

% OPTION 1: New preferred area is inherited from mother
% newg.A0 = cat(1, newg.A0, gcur.A0(ci));

% OPTION 2: (Use for E ~ A0/A model)
% New preferred area is ~a quarter of the mother
% gcur.A0(ci) =  gcur.A0(ci)/(2*sqrt(2));

% OPTION 3: (Use for E ~ (A-A0)^2 model)
% New preferred area is half the area of the mother
gcur.A0(ci) =  gcur.A0(ci)/2;

newg.A0 = cat(1,gcur.A0, gcur.A0(ci));

% For some reason Idse recommends doubling the bulk modulus for the
% E ~ (A-A0)^2 model. Why? Ignoring this for now...
% gcur.kA0(ci) = 2*gcur.kA0(ci);
newg.kA0 = cat(1, newg.kA0, gcur.kA0(ci));

% Update All Other Cell Perimeter Parametrs -------------------------------

% New preferred perimeter is inherited from mother
newg.p0 = cat(1,newg.p0, gcur.p0(ci));

% Perimeter tension is inherited from mother
% gcur.pT0(ci) = gcur.pT0(ci);
newg.pT0 = cat(1, gcur.pT0, gcur.pT0(ci));

% Stress is inherited from mother
% This prevents stress response in dividing cells
newg.stress = cat(1, newg.stress, gcur.stress(ci));

% New tension is inherited from mother:
% Intersected bonds will have the same tension on each side
% The partition itself will have the mean of the bonds it intersects
Told0 = newg.T0(ibi);
newg.T0 = cat(1, newg.T0, ones([6 1]));
newg.T0(nbi(1:2)) = mean(Told0);
newg.T0(nbi(3:4)) = Told0(1);
newg.T0(nbi(5:6)) = Told0(2);

% New preferred edge length is inherited in the same way as the tension
l0old = newg.l0(ibi);
newg.l0 = cat(1, newg.l0, ones([6 1]));
newg.l0(nbi(1:2)) = mean(l0old);
newg.l0(nbi(3:4)) = l0old(1);
newg.l0(nbi(5:6)) = l0old(2);

% Inherit clonal identity;
newg.clones = cat(1, newg.clones, gcur.clones(ci));

end