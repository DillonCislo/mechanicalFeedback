clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));
addpath(genpath('/Users/idse/repos/epitheliumAnalysis'));

% import data locations file
dataLocations;

bi = 7;
fi = 2;

disp(description{bi})
disp(bareFnames{bi}{fi})

dataDir = fullfile(dataDir{bi},['analysis_' bareFnames{bi}{fi}]);
bareFnames = bareFnames{bi};

t = 1;

%%
% load apical surface

% SOI
SOIdir = fullfile(dataDir, [bareFnames{fi} '_apicalSOI']);
apicalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

% clone segmentation
fname = [fullfile(dataDir,'tripleLayerSeg', bareFnames{fi}), '_cloneMaskApical.tif'];
cloneMaskApical = imread(fname);
cloneBdryApical = bwboundaries(cloneMaskApical');

fname = [fullfile(dataDir, 'tripleLayerSeg', bareFnames{fi}), '_apicalL.tif'];
cloneApicalL = imread(fname);

% membrane segmentation
fname = [fullfile(dataDir, 'tripleLayerSeg', bareFnames{fi}), '_apicalSeg.tif'];
apicalSeg = imread(fname);

% vertex model
fname = fullfile(dataDir, 'tripleLayerSeg', 'apicalVertexModel');
S = load(fname);
cellLayer = S.cellLayer;
clear S;

%% visualize

clf
bdryBondIdx = cellLayer.getStateBoundary(t,'clone',1);
cloneBdryBonds = cellLayer.bonds{t}(bdryBondIdx(:,1));

areas = cellLayer.getCellGeometry(t,'area');
options = struct('colorTable', areas, 'transparent', 0);
cellLayer.visualize(t, options);

for i = 1:numel(cloneBdryBonds)
    V = cellLayer.vertices{t}(cloneBdryBonds(i).vertInd,1:2);
    hold on 
    plot(V(:,1),V(:,2),'-k','LineWidth',3);
    hold off
end
axis equal

%% fake shifted clones

cellCM = cellLayer.getCellGeometry(t, 'centroid');
CMsub = round(cellCM)';
linind = sub2ind(size(cloneMaskApical), CMsub(:,2), CMsub(:,1));

cloneMaskApicalShift = false(size(cloneMaskApical));
shift = 100;
cloneMaskApicalShift(:, (shift+1):end) = cloneMaskApical(:,1:(end-shift));
fakecloneL = bwlabel(cloneMaskApicalShift);
infakeclone = uint16(fakecloneL(linind));

%% visualize

inclone = cat(1,cellLayer.cells{t}.state);

figure
apicalarm = apicalSOI.data.patches{1}.apply{3};
imshow(apicalarm,[])
lw = 1;
hold on

options = struct('edgeColor', 'b', 'colorTable', inclone, 'transparent', 0);
cellLayer.visualize(t, options);

shift = -300;
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1) + shift, cloneBdryApical{i}(:,2), '-r','LineWidth',lw)
end

%% myosin or jub with clone outline

field = apicalSOI.getField('data');
pbs = field(t).getPatch('xy_index').apply;
for i = 1:3, pbs{i} = mat2gray(pbs{i}); end;
discProperImage = mat2gray(pbs{1});

h = figure;
imshow(discProperImage, [], 'InitialMagnification', 50);
hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-r')
end
hold off

%% plot distributions

myoOrJub = 'jub';

figure,
% area stats
stats = {};
cellArea = cellLayer.getCellGeometry(t, 'area');
stats{1}.area = makeDist(cellArea(inclone > 0), 30);
stats{2}.area = makeDist(cellArea(infakeclone > 0), stats{1}.area.bins);

row = 1;
nrows = 2;
ncols = 2;
plotInfo = struct('plotTitle', 'cell area');
fieldname = 'area';
subplotDist(nrows, ncols, row, plotInfo, stats, fieldname, {'clone','WT'})
                
% jub stats
ci = 1;
stats{1}.jub = makeDist(pbs{ci}(cloneMaskApical), 30);
stats{2}.jub = makeDist(pbs{ci}(cloneMaskApicalShift), stats{1}.jub.bins);

row = 2;
nrows = 2;
ncols = 2;
plotInfo = struct('plotTitle', 'Jub intensity');
fieldname = myoOrJub;
subplotDist(nrows, ncols, row, plotInfo, stats, fieldname, {'clone','WT'})
  
fname = [fullfile(dataDir, bareFnames{fi}), '_' myoOrJub '_stats.png'];
saveas(gcf, fname);


%%
%----------------------------------------------------------------------
% NOT FUNCTIONAL BELOW
%----------------------------------------------------------------------

% start for determining graph connected components for cleanup

dualLayer = cellLayer.dualize;
dualLayer.visualize(1)

% looks funny because vertices are not ordered
% cellLayer2 = dualLayer.dualize;
% cellLayer2.visualize;
% this also doesn't work

%%
CL = dualLayer;
NCells = numel(CL.cells{1});
NVertices = size(CL.vertices{1},1);

VorCells = {};
for i = 1:NCells
    tmp = cat(1,CL.bonds{t}(CL.cells{t}(i).bondInd).vertInd);
    VorCells{i} = tmp(:,1);b
end

S=sparse(NVertices,NVertices);
for i = 1 : NCells
    for j = 1 : (length(VorCells{i})-1)             
        v1 = VorCells{i}(j);
        v2 = VorCells{i}(j+1);
        S(v1,v2) = 1;
        S(v2,v1) = 1;
    end
end

[labels ~] = graph_connected_components(S);
n = hist(labels)

%% write function cellLayer.removeCells, cellLayer.removeVerts

CL = cellLayer;
CL.removeCells(t, 1:500);
figure,
CL.visualize(1);

% REMOVE UNUSED VERTS

