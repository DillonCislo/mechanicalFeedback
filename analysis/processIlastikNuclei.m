clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));

% import data locations file
dataLocations; 

bi = 6;
fi = 1;

disp(description{bi})
disp(bareFnames{bi}{fi})

dataDir = dataDir{bi};
bareFnames = bareFnames{bi};

tldir = fullfile(dataDir, ['analysis_' bareFnames{fi}], 'tripleLayerSeg');
if ~exist(tldir)
   mkdir(tldir) 
end

%% read the data

i = 1;
fname = [fullfile(dataDir, bareFnames{i}), '.RGB_Seg.tiff'];

disp(bareFnames{i});

info = imfinfo(fname);
im = zeros([info(1).Height, info(1).Width numel(info)], 'uint8');

for zi = 1:size(im,3)
    im(:,:,zi) = imread(fname, zi);
end

im = im == 1;

% %% crop 
%im2 = im(230:780,:,10:350);
im2 = im;

%% imopen; will take about 4 minutes on a full disc
tic
im2 = imopen(im2, strel3D('sphere',5));
toc

%% mask the columnar cells within the folds

% SOIdir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_apicalSOI']);
% apicalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);
% 
% SOIdir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_basalSOI']);
% basalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

% [~,~,Z] = meshgrid(1:size(im2,1),1:size(im2,2),1:size(im2,3));
% 
% diskmask = false(size(im2));
% Zapical = apicalSOI.embedding.patches{1}.apply{3};
% Zbasal = basalSOI.embedding.patches{1}.apply{3};
% for i = 1:size(diskmask,3)
%     diskmask(:,:,i) = i > Zapical & i < Zbasal;
% end
% imshow(squeeze(diskmask(:,400,:)),[]);
%
%im2 = im2 & diskmask;

%% save segmentation

fname = [fullfile(dataDir, bareFnames{i}), '.RGB_Seg_Clean.tiff'];

imwrite(mat2gray(im2(:,:,1) > 0), fname); 
for i = 2:size(im,3)
    imwrite(mat2gray(im2(:,:,i) > 0), fname, 'WriteMode', 'append'); 
end

%% color label
% 
% cmap = hsv(15);
% cmap = repmat(cmap, [ceil(CC.NumObjects/15) 1]);
% rgbL = label2rgb3d(L,cmap,[0 0 0],'shuffle'); 
% 
% imwrite(0.5*squeeze(rgbL(:,:,1,:)), 'colorL.tif'); 
% for i = 2:size(im,3)
%     imwrite(0.5*squeeze(rgbL(:,:,i,:)), 'colorL.tif', 'WriteMode', 'append'); 
% end

%% read clones

fname = [fullfile(dataDir, bareFnames{fi}), '_scaleSS.tif'];
info = imfinfo(fname);

nChannels = 4;
nSlices = size(info,1)/nChannels;
cloneChannel = 3;
clones = zeros([info(1).Height, info(1).Width nSlices], 'uint8');

for zi = 1:size(clones,3)
    clones(:,:,zi) = imread(fname, 3 + (zi-1)*4);
end

%%
% connected components and label matrix
tic
CC = bwconncomp(im2);
toc
tic
L = labelmatrix(CC);
toc

%% nuclear centroids and clone marker level

clI = zeros([1 CC.NumObjects]);

for i = 1:CC.NumObjects
    clI(i) = mean(clones(CC.PixelIdxList{i}));
end

stats = regionprops(CC, 'centroid');
CM = cat(1,stats.Centroid);

nuclearStats = struct('centroid', CM, 'cloneMarker', clI);
save(fullfile(dataDir, ['analysis_' bareFnames{fi}], 'tripleLayerSeg', 'nuclearStats'), 'nuclearStats');

%% plot nuclear centroids 

clf
zmin = 15;
zmax = nSlices - 15;

apicalGarbage = find(CM(:,3) < zmin);
idx = CM(:,3) > zmin & CM(:,3) < zmax;

% height coloring
scatter3(CM(idx,1),CM(idx,2),CM(idx,3),20,CM(idx,3),'o')

view([0 0 1]);
colormap jet
axis equal






