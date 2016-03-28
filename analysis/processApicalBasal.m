clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));
addpath(genpath('/Users/idse/repos/epitheliumAnalysis'));

% import data locations file
dataLocations;

bi = 11; %bi = 1;
fi = 2; %fi = 3;

% we don't store channelsUsed so here I assume it's everything except DNA 
setlabels = labels{bi};
setlabels(strcmp(setlabels,'DNA')) = [];
cloneC = find(strcmp(setlabels,'clone'));
memC = find(strcmp(setlabels,'arm'));

disp(description{bi})
disp(bareFnames{bi}{fi})

dataDir = dataDir{bi};
bareFname = bareFnames{bi}{fi};

tidx = 1;

tldir = fullfile(dataDir, ['analysis_' bareFname], 'tripleLayerSeg');
if ~exist(tldir)
   mkdir(tldir) 
end

resultsDir = fullfile(dataDir, ['analysis_' bareFname]);

%% 
%----------------------------------------------------------------------
% APICAL
%----------------------------------------------------------------------

% load apical surface

SOIdir = fullfile(dataDir, ['analysis_' bareFname], [bareFname '_apicalSOI']);
apicalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

%%
% apical clone mask
%-------------------

clthresh = 20;

% apical
field = apicalSOI.getField('data');
clonesApical = field.patches{1}.apply{cloneC};

cloneMaskApical = bwareaopen(clonesApical > clthresh, 15);
cloneMaskApical = imclose(cloneMaskApical,strel('disk',5));
%cloneMaskApical = imfill(cloneMaskApical,'holes');
cloneMaskApical = ~bwareaopen(~cloneMaskApical, 400);
cloneMaskApical = bwareaopen(cloneMaskApical, 500);

% % special for anticlones
% cloneMaskApical = imopen(~cloneMaskApical & clonesApical > 0, strel('disk',10));
% cloneMaskApical = imfill(cloneMaskApical,'holes');

% boundary
cloneBdryApical = bwboundaries(cloneMaskApical');

% visualize
imshow(clonesApical,[])
lw = 1;
hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-g','LineWidth',lw)
end
hold off

%%
imshow(clonesApical,[])

%%
% save clones

imshow(clonesApical,[])
lw = 1;
hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-g','LineWidth',lw)
end
hold off

fname = [fullfile(dataDir, ['analysis_' bareFname], 'tripleLayerSeg', bareFname), '_apicalClones.tif'];
saveas(gcf, fname);

fname = [fullfile(dataDir, ['analysis_' bareFname], 'tripleLayerSeg', bareFname), '_cloneMaskApical.tif'];
imwrite(cloneMaskApical, fname);

close;

%% load membrane segmentation 

SOIdir = fullfile(dataDir, ['analysis_' bareFname], [bareFname '_apicalSOI']);
segfile = fullfile(SOIdir, 'fields', 'data_MIP', 'xy_index', 'xy', 'cmp_1_3_T0000_seg.tif');
%segfile = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_apicalSurfaceForSeg_seg.tif']);
apicalmem = imread(segfile)==2;

% keep only largest connected component
CC = bwconncomp(apicalmem);
stats = regionprops(CC,'Area');
maxidx = find([stats.Area] == max([stats.Area]));
memclean = false(size(apicalmem));
memclean(CC.PixelIdxList{maxidx}) = true;

fname = [fullfile(dataDir, ['analysis_' bareFname], 'tripleLayerSeg', bareFname), '_apicalSeg.tif'];
imwrite(memclean, fname);

imshow(memclean)

%%
% identify enormous cells as background
CC = bwconncomp(~memclean);
stats = regionprops(CC,'area');
memcleaner = int8(memclean);
bgregions = find([stats.Area]>4000);
for bgi = 1:numel(bgregions)
    memcleaner(CC.PixelIdxList{bgregions(bgi)}) = -1;
end

% remove border cells
border = ~imclearborder(memcleaner == 0) & ~memcleaner;
memcleaner(border) = -1;

imshow(memcleaner,[])

%% turn into vertex model

% create a CellLayer object called cellLayer
cellStateVars = {'clone'};
bondStateVars = {};
nTimePts = 1;
cellLayer = CellLayer(nTimePts, cellStateVars, bondStateVars);

% initialize timepoint from Lattmin
% FOR SOME REASON minVertexDist == 0 WILL MESS UP BONDLABELS
% I SHOULD FIGURE OUT WHY
t = 1;
options = struct('closeSize', 1, 'areaOpenSize',30, 'minVertexDist', 1,...
                    'bgCells', 0, 'trim', 0);
tic
cellLayer.initTime(t, 'image', memcleaner, options);
cellLayer.calcCellGeometry(t);
toc

% % HERE SOMETHING GOES WRONG
% cellLayer.bonds{1}(cellLayer.cells{1}(37).bondInd).cellInd

%% visualize

figure,

field = apicalSOI.getField('data');
pbs = field(tidx).getPatch('xy_index').apply{memC};

imshow(pbs, [], 'InitialMagnification', 50);
fname = [fullfile(tldir, bareFname), '_apicalVertexModelOff.tif'];
saveas(gcf, fname);

hold on
options = struct('transparent','all','edgeColor','green','lineWidth',1);
cellLayer.visualize(1,options)
hold off

fname = [fullfile(tldir, bareFname), '_apicalVertexModel.tif'];
saveas(gcf, fname);

%% determine clones

cellCM = cellLayer.getCellGeometry(t, 'centroid');
CMsub = round(cellCM)';
cloneL = bwlabel(cloneMaskApical);
linind = sub2ind(size(cloneMaskApical), CMsub(:,2), CMsub(:,1));
%inclone = cloneMaskApical(linind);
inclone = uint16(cloneL(linind));
cellLayer.setCellState(t, inclone);

%% visualize

clf
imshow(pbs);

areacutoff = 600;
areas = uint32(cellLayer.getCellGeometry(t,'area'));
bdryBondIdx = cellLayer.getStateBoundary(t,'clone',1);
cloneBdryBonds = cellLayer.bonds{t}(bdryBondIdx(:,1));
areas(areas > areacutoff) = areacutoff;
colormap jet;

hold on 
options = struct('colorTable', areas, 'transparent', 0, 'colorbar', true);
cellLayer.visualize(t, options);

for i = 1:numel(cloneBdryBonds)
    V = cellLayer.vertices{t}(cloneBdryBonds(i).vertInd,1:2);
    plot(V(:,1),V(:,2),'-k','LineWidth',3);
end
hold off
axis equal

fname = [fullfile(tldir, bareFname), '_areaMap.tif'];
saveas(gcf, fname);

%% save cellLayer

fname = fullfile(tldir, 'apicalVertexModel');
save(fname, 'cellLayer');

%% identify clones to keep, only if no basal surface

% HOW: 
% - click on clone, press enter
% - repeat until no desired clones are left
% - then click outside clones and press enter

% select special clones apically
bw1 = bwlabel(cloneMaskApical);
newlabels1 = zeros(size(cloneMaskApical),'uint16');
label = 1;
clickedobject = 1;

while clickedobject

    apicalPB = apicalSOI.data.patches{1}.apply{memC};
    imshow(cat(3,mat2gray(apicalPB),bw1,0*bw1));
    [x1,y1] = getpts(gcf);
    clickedobject = bw1(round(y1(1)),round(x1(1)));

    % a dumb way to get the label for the selected clones
    % (first 2 arguments of bwselect are xdata ydata)
    cl1 = bwselect(bw1,x1,y1);
    newlabels1(cl1) = label;
    bw1(cl1) = 0;
    label = label + 1;
end

cloneApicalL = newlabels1;
fname = [fullfile(resultsDir, 'tripleLayerSeg', bareFname), '_apicalL.tif'];
imwrite(cloneApicalL, fname);
close;

%%
%----------------------------------------------------------------------
% BASAL
%----------------------------------------------------------------------

% load the basal surface
SOIdir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_basalSOI']);
basalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

%% basal clone mask

clthresh = 250;

% basal
field = basalSOI.getField('data_SIP');
clonesBasal = field.patches{1}.apply{cloneC};
cloneMaskBasal = bwareaopen(clonesBasal > clthresh, 15);
cloneMaskBasal = imclose(cloneMaskBasal,strel('disk',5));
%cloneMaskBasal = imfill(cloneMaskBasal,'holes');
cloneMaskBasal = ~bwareaopen(~cloneMaskBasal,400);
cloneMaskBasal = bwareaopen(cloneMaskBasal,100);
cloneMaskBasal = imopen(cloneMaskBasal,strel('disk',4));

% % special for anticlones
% cloneMaskBasal = imopen(~cloneMaskBasal & clonesBasal > 0, strel('disk',10));
% cloneMaskBasal = imfill(cloneMaskBasal,'holes');
% imshow(cloneMaskBasal)

cloneBdryBasal = bwboundaries(cloneMaskBasal');

imshow(clonesBasal,[])
lw = 1;
hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-r','LineWidth',lw)
end
for i = 1:numel(cloneBdryBasal)
    plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2), '-g','LineWidth',lw)
end
hold off

% save

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], 'tripleLayerSeg', bareFnames{fi}), '_basalClones.tif'];
saveas(gcf, fname);

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], 'tripleLayerSeg', bareFnames{fi}), '_cloneMaskBasal.tif'];
imwrite(cloneMaskBasal, fname);

%% match apical and basal

% MAKE IT SO THAT APICAL KEEPS THE SAME LABEL, 
% SO WE CAN MATCH THINGS UP LATER

bw1 = cloneMaskApical;
bw2 = cloneMaskBasal;

% click clone in green then enter (multiple clicks means merge for
% landmark), then click clone in red and enter
% to exit the loop, click outside the clones for both points

label = 1;

newlabels1 = zeros(size(cloneMaskApical),'uint16');
newlabels2 = newlabels1;

imshow(cat(3,cloneMaskApical, 0*cloneMaskApical, cloneMaskBasal));

% click outside clone to break out of loop
clickedobject = 1;
while clickedobject
    
    imshow(cat(3,bw1, 0*bw2, bw2));
    [x1,y1] = getpts(gcf);
    [x2,y2] = getpts(gcf); 
    
    clickedobject = bw1(round(y1(1)),round(x1(1)));
    
    % a dumb way to get the label for the selected clones
    % (first 2 arguments of bwselect are xdata ydata)
    cl1 = bwselect(bw1,x1,y1);
    cl2 = bwselect(bw2,x2,y2);

    newlabels1(cl1) = label;
    newlabels2(cl2) = label;
    
    bw1(cl1) = 0;
    bw2(cl2) = 0;
    
    label = label + 1;
end

cloneApicalL = newlabels1;
cloneBasalL = newlabels2;

cloneBdryBasal = bwboundaries((cloneBasalL > 0)');
cloneBdryApical = bwboundaries((cloneApicalL > 0)');

%% save matched apical and basal

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], 'tripleLayerSeg', bareFnames{fi}), '_basalL.tif'];
imwrite(cloneBasalL, fname)

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], 'tripleLayerSeg', bareFnames{fi}), '_apicalL.tif'];
imwrite(cloneApicalL, fname)





%%
% ATTEMPT BASAL SEGMENTATION

% load membrane segmentation 

segfile = [fullfile(dataDir, ['analysis_' bareFnames{fi}], 'basalSegTests', 'actinAVG10_20_seg.tif')];
basalmem = imread(segfile)==2;

% keep only largest connected component
CC = bwconncomp(basalmem);
stats = regionprops(CC,'Area');
maxidx = find([stats.Area] == max([stats.Area]));
memclean = false(size(basalmem));
memclean(CC.PixelIdxList{maxidx}) = true;

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], 'tripleLayerSeg', bareFnames{fi}), '_basalSeg.tif'];
imwrite(memclean, fname);

imshow(memclean)

%%
% identify enormous cells as background
CC = bwconncomp(~memclean);
stats = regionprops(CC,'area');
memcleaner = int8(memclean);
bgregions = find([stats.Area]>4000);
for bgi = 1:numel(bgregions)
    memcleaner(CC.PixelIdxList{bgregions(bgi)}) = -1;
end

% remove border cells
border = ~imclearborder(memcleaner == 0) & ~memcleaner;
memcleaner(border) = -1;

imshow(memcleaner,[])

%% turn into vertex model

% create a CellLayer object called cellLayer
cellStateVars = {'clone'};
bondStateVars = {};
nTimePts = 1;
cellLayer = CellLayer(nTimePts, cellStateVars, bondStateVars);

% initialize timepoint from Lattmin
% FOR SOME REASON minVertexDist == 0 WILL MESS UP BONDLABELS
% I SHOULD FIGURE OUT WHY
t = 1;
options = struct('closeSize', 0, 'areaOpenSize',0, 'minVertexDist', 1,...
                    'bgCells', 0, 'trim', 0);
tic
cellLayer.initTime(t, 'image', memcleaner, options);
cellLayer.calcCellGeometry(t);
toc

% % HERE SOMETHING GOES WRONG
% cellLayer.bonds{1}(cellLayer.cells{1}(37).bondInd).cellInd