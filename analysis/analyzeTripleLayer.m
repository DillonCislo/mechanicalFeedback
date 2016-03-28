clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));
addpath(genpath('/Users/idse/repos/epitheliumAnalysis'));

% import data locations file
dataLocations;

bi = 7; %bi = 1;
fi = 1; %fi = 3;

disp(description{bi})
disp(bareFnames{bi}{fi})

% we don't store channelsUsed so here I assume it's everything except DNA 
setlabels = labels{7};
setlabels(strcmp(setlabels,'DNA')) = [];
cloneC = find(strcmp(setlabels,'clone'));
memC = find(strcmp(setlabels,'arm'));

dataDir = fullfile(dataDir{bi}, ['analysis_' bareFnames{bi}{fi}]);
bareFnames = bareFnames{bi};

tidx = 1;

%%
%-------------------
% load surfaces
%-------------------

% load apical surface

SOIdir = fullfile(dataDir, [bareFnames{fi} '_apicalSOI']);
apicalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

fname = [fullfile(dataDir,'tripleLayerSeg', bareFnames{fi}), '_cloneMaskApical.tif'];
cloneMaskApical = imread(fname);
cloneBdryApical = bwboundaries(cloneMaskApical');

fname = [fullfile(dataDir, 'tripleLayerSeg', bareFnames{fi}), '_apicalL.tif'];
cloneApicalL = imread(fname);

fname = [fullfile(dataDir, 'tripleLayerSeg', bareFnames{fi}), '_apicalSeg.tif'];
apicalSeg = imread(fname);

fname = fullfile(dataDir, 'tripleLayerSeg', 'apicalVertexModel');
S = load(fname);
cellLayer = S.cellLayer;
clear S;

% load the basal surface
SOIdir = fullfile(dataDir, [bareFnames{fi} '_basalSOI']);
basalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

fname = [fullfile(dataDir, 'tripleLayerSeg', bareFnames{fi}), '_cloneMaskBasal.tif'];
cloneMaskBasal = imread(fname);
cloneBdryBasal = bwboundaries(cloneMaskBasal');

fname = [fullfile(dataDir, 'tripleLayerSeg', bareFnames{fi}), '_basalL.tif'];
cloneBasalL = imread(fname);

%%
%--------------------
% nuclei 
%--------------------

% load nuclear centroids 
nuclearStats = load(fullfile(dataDir, 'tripleLayerSeg', 'nuclearStats'), 'nuclearStats');
nuclearStats = nuclearStats.nuclearStats;

nucX = nuclearStats.centroid(:,1);
nucY = nuclearStats.centroid(:,2);
nucZ = nuclearStats.centroid(:,3);

% get the nuclei between apical and basal surface
basalZ = basalSOI.embedding.patches{1}.apply{3};
apicalZ = apicalSOI.embedding.patches{1}.apply{3};

nucI = round(nuclearStats.centroid(:,2));
nucJ = round(nuclearStats.centroid(:,1));

nucLinInd = sub2ind(size(basalZ),nucI,nucJ);

goodNucIdx = nucZ > apicalZ(nucLinInd) & nucZ < basalZ(nucLinInd);

% height coloring
idx = goodNucIdx;
scatter3(nucX(idx),nucY(idx),nucZ(idx),20,nucZ(idx),'o')
view([0 0 1]);
axis equal
colormap jet

%% assign nuclei to specific clones

% nuclei
nucCloneIdx = (nuclearStats.cloneMarker > 4)';
idx = goodNucIdx & nucCloneIdx;

% nuclei spread out more than the apical surface of the clone, so we have
% to expand the apical mask
% watershed so that clones can touch in a sensible way
cloneL = cloneApicalL;
D = bwdist(cloneL > 0);
cloneLexpand = watershed(D).*uint8(D < 70);
%imshow(D,[])

% assign nuclei to clone
cloneNucSub = round(nuclearStats.centroid(idx,:));
cloneNucInd = sub2ind([size(D,1),size(D,2)],...
                                    cloneNucSub(:,2), cloneNucSub(:,1));
nucCloneLabel = cloneLexpand(cloneNucInd);

% determine number of nuclei for each label
nClones = max(cloneL(:));
[nnuclei,~] = hist(nucCloneLabel,double(1:nClones));

% remove clonal nuclei that cannot be assigned to apical clone area 
nucCloneLabelFull = zeros(size(nucX),'uint8');
nucCloneLabelFull(idx) = nucCloneLabel;
idx = nucCloneLabelFull > 0;

scatter3(nucX(idx),nucY(idx),nucZ(idx),20,nucCloneLabelFull(idx),'o');
view([0 0 1]);
colormap jet
axis equal

%% visualize

figure
apicalarm = apicalSOI.data.patches{1}.apply{3};
imshow(apicalarm,[])
lw = 2;
hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-r','LineWidth',lw)
end
for i = 1:numel(cloneBdryBasal)
    plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2), '-b','LineWidth',lw)
end

scatter3(nucX(idx),nucY(idx),nucZ(idx),50,'.g');
view([0 0 1]);
axis equal
hold off

fname = [fullfile(dataDir, bareFnames{fi}), '_cloneNuclei.tif'];
saveas(gcf, fname);

%% myosin with clone outline

lw = 2;
field = apicalSOI.getField('data');
pbs = field(tidx).getPatch('xy_index').apply;
for i = 1:3, pbs{i} = mat2gray(pbs{i}); end;
apicalPB = cat(3,pbs{[3 1 2]});
apicalPB = mat2gray(pbs{1});

h = figure;
imshow(apicalPB, [], 'InitialMagnification', 50);
hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-r','LineWidth',lw)
end
for i = 1:numel(cloneBdryBasal)
    plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2), '-b','LineWidth',lw)
end
hold off

fname = [fullfile(dataDir, bareFnames{fi}), '_apicalMyo.tif'];
saveas(h, fname);

%% height map with clone outline

height = basalZ - apicalZ;
zmin = min(height(:));
zmax = max(height(:));
%height(isnan(height)) = 0;
%imshow(height,[zmin zmax])

h = figure
imshow(height, [zmin zmax])
colormap jet
shading interp
axis equal
view([0 0 1]);
colorbar

hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-k', 'LineWidth',2)
    plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2), '-k', 'LineWidth',1)
end
hold off

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_heightMap.tif'];
saveas(h, fname);

%% density map

[H,X,Y] = hist2d(nuclearStats.centroid(:,1:2), 64,64, [1 size(clonesApical,1)], [1 size(clonesApical,2)]);
% H2 = imfilter(H,fspecial('gaussian',[15 15], 2));
% H2 = imresize(H2, [size(clonesApical,1) size(clonesApical,2)]);
H2 = imfilter(imresize(H, [size(clonesApical,1) size(clonesApical,2)]),fspecial('gaussian',[80 80], 20));

figure, imshow(H2,[])

colorbar
hold on
for i =1:nClones
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2),'LineStyle','-','LineWidth',2,'Color', 'red');
    plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2),'LineStyle','-','LineWidth',2,'Color', 'blue');
end
hold off

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_density.tif'];
saveas(gcf, fname);

%% cell volume map

figure,
imshow(height./H2,[])
colorbar
lw = 1;
hold on
for i =1:nClones
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2),'LineStyle','-','LineWidth',lw,'Color', 'red');
    plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2),'LineStyle','-','LineWidth',lw,'Color', 'blue');
end
hold off

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_volume.tif'];
saveas(gcf, fname);

%% clone volume

statsA = regionprops(cloneApicalL,'Area');
statsB = regionprops(cloneBasalL,'Area');

%% apical and basal density in clone

i = 2;
statsA(i)
statsB(i)
nnuclei(i)
nnuclei(i)./statsA(i).Area
nnuclei(i)./statsB(i).Area

%% total

txtfile = [fullfile(dataDir, bareFnames{fi}), '_stats.txt'];
if exist(txtfile)
    delete(txtfile);
end
diary(txtfile)

apicalArea = cat(1,statsA.Area);
basalArea = cat(1,statsB.Area);
areaRatio = apicalArea./basalArea;
basalDensity = (nnuclei'./basalArea);
apicalDensity = (nnuclei'./apicalArea);

totBasalDensity = sum(nnuclei)./sum(basalArea)
totApicalDensity = sum(nnuclei)./sum(apicalArea)
totAreaRatio = sum(apicalArea)./sum(basalArea)

avgAreaRatio = [mean(areaRatio) std(areaRatio)];
avgBasalDensity = [mean(basalDensity) std(basalDensity)];
avgApicalDensity = [mean(apicalDensity) std(apicalDensity)];

areaRatio
basalDensity 
apicalDensity
basalArea

% NAIVE DENSITIES

WTMaskApical =  pbs{1} > 0 & ~cloneMaskApical;
linNucInd = sub2ind(size(cloneMaskApical),round(nucY), round(nucX));

% naive densities
disp('- naive WT -');
sum(WTMaskApical(linNucInd))./sum(WTMaskApical(:))

disp('- naive clone -');
sum(cloneMaskApical(linNucInd))./sum(cloneMaskApical(:))

disp('- per clone naive apical -');
for i = 1:max(cloneApicalL(:))
    mask = cloneApicalL == i;
    sum(mask(linNucInd))./sum(mask(:))
end

disp('- per clone naive basal -');
for i = 1:max(cloneBasalL(:))
    mask = cloneBasalL == i;
    sum(mask(linNucInd))./sum(mask(:))
end

diary off

%% cell count apically 

napical = hist(double(cat(1,cellLayer.cells{t}.state)));

%% count nuclei in some polygonal mask

apicalPB = cat(3,pbs{[3 1 2]});
[~, xi, yi] = roipoly(apicalPB);
IN = inpolygon(nucX,nucY,xi,yi);
A = polyarea(xi,yi);
roughDensity = sum(IN)/A


