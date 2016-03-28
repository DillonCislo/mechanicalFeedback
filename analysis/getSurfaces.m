clear all; close all;

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(scriptPath));

dataLocations; % particular to this case

bi = 10;
fi = 1;

disp(description{bi})
disp(bareFnames{bi}{fi})

dataDir = dataDir{bi};
bareFnames = bareFnames{bi};
detectionChannel = detectionChannel{bi};

projectDir = fullfile(dataDir, 'projectFiles');

xp = project.Experiment(projectDir, dataDir);

if ~exist(fullfile(dataDir, ['analysis_' bareFnames{fi}]),'dir')
    mkdir(fullfile(dataDir, ['analysis_' bareFnames{fi}]));
end

%%
% Next we set the metadata pertaining to the raw data files in the structure
% fileMeta. ImSAnE assumes that timepoints are saved individually and that 
% filenames in a timeseries are identical up to an integer specifying the 
% timepoint. Therefore we have
%
% * filenameFormat:               Filename, %u in the position of the integer
% * timePoints = [t1, t2, ..] :   List of times available. In this example 
% we have just a single time point 0.
% * stackResolution :             Stack resolution in micron.

fileMeta = struct();
fileMeta.dataDir = dataDir;
fileMeta.filenameFormat = [bareFnames{fi} '_scaleSS.tif'];
fileMeta.timePoints = [0];
fileMeta.stackResolution = [1 1 1]; 

%% 
% In the structure expMeta we set general parameters for the surface
% analysis we will do. 
%
% * channelsUsed:   Which channels do we need.
% * channelColor:   Assign color to the channels, RGB = [1 2 3].
%                   In this example the first channel is E-cadherin, and 
%                   the second is actin. We want these in green and red,
%                   respectively.
% * dynamicSurface: Does the surface shape change with time?  
%                   For a single time point this is false. True is not yet
%                   supported.
% * jitterCorrection:   Not needed here.
% * detectorType:       Which type of surface detector will be used.
%                       We will look at only one side of the wing, which is
%                       a planar surface so we use planarDetector.
% * fitterType:         Which type of fitter will be used.
%                       We fit planar surfaces using Thin Plate Spline:
%                       tpsFitter.

expMeta = struct();
expMeta.description = 'Bantam clones';
expMeta.channelsUsed = [1 2 3];
expMeta.channelColor = [2 3 1];
expMeta.dynamicSurface = false;
expMeta.jitterCorrection = false;
expMeta.detectorType = 'surfaceDetection.planarEdgeDetector';
expMeta.fitterType = 'surfaceFitting.tpsFitter';

xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);

%%
% Finally we call initNew(), which reads the stack size from the first 
% available time point, then initializes fitter and detector and creates 
% fitOptions and detectOptions based on their defaults.

xp.initNew();

%% Load a time point from the data
%
% Now that we have set up the project, we can load a time point.
% loadTime sets xp.currentTime, loads the stack into xp.stack 
% and resets detector and fitter with the default options for that time.
%
% We rescale to unit aspect ratio to make detection work better later and
% visualize in the right proportions.

xp.loadTime(0);
xp.rescaleStackToUnitAspect();

%% 
% xp.stack is not an array but a Stack object.
% The easy way to look at a slice through the data is using getSlice.

imshow(xp.stack.getSlice('x', 100), []);

%% Mask data
%
% Detection speed depends on stack size. We can speed it up by some basic
% masking. Masking out unwanted information can also improve detection and
% fitting.
% Rather than creating 3D masks, we create 2D masks that are applied to
% each cross section. In this case we simply cut out the bottom half of the
% stack because we want to detect the top of the wing. We are doing this in
% cross section normal to the y direction. 

zab = 270; % z-apical-basal separator

projectionMask = cell([3 1]);
projectionMask{2} = true([xp.stack.imageSize(1) xp.stack.imageSize(3)]);
projectionMask{2}(:, 1:zab) = false;

xp.stack.setProjectionMask(projectionMask);

imshow(xp.stack.getSlice('y', 500), []);

% %% speckle removal?
% 
% ims = xp.stack.image.apply;
% tic
% nospeckle = ims{1};
% for zi = zab:xp.stack.imageSize(3)
%     fprintf('.');
%     nospeckle(:,:,zi) = medfilt2(squeeze(ims{1}(:,:,zi)),[15 15]);
% end
% fprintf('\n');
% toc
% imshow(imadjust(squeeze(ims{1}(:,890,:))),[])

%% Detect the surface
%
% planarDetector.detectSurface detects the surface as the position of the 
% maximal Gaussian z-derivative in some direction, i.e. the position of the
% largest intensity jump along some direction and smoothened over some
% scale.
%
% A number of detection options directly affect detection:
%
% * sigma :     Width of the Gaussian z-derivative.
% * channels :  Channels (summed) to use for detection.
% * zdir :      Dimension corresponding to z, minus flips direction.
% Flipping the direction can sometimes improve detection.
%
% Then there are options which filter the result and can be modified
% without redetecting:
%
% * maxIthresh:     Throw out points with MIP dimmer than this.
% * summedIthresh:  Throw out points with SIP dimmer than this.
% * sigZoutliers:   Remove height outliers after all other masks.
% * scaleZoutliers: Spatial scale of outlier removal.
%
% scaleZoutliers is the linear size of a region over which the
% distribution of height is computed, sigZoutliers is then a cutoff in
% units of standard deviation of this distribution to remove misdetected
% points far above or below the other points in the region.

detectOptions = struct(  'sigma', 6, 'channels', [1], 'zdir', 3,...
                        'maxIthresh', 0.02, 'summedIthresh', 0,...
                        'sigZoutliers', 2, 'scaleZoutliers', 50); 

% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.

xp.setDetectOptions(detectOptions);
xp.detectSurface();

imshow(xp.detector.surfaceMatrix, [zab xp.stack.imageSize(3)],...
                                            'InitialMagnification', 40);
                                        
%%

inspectOptions= struct('dimension', 'x', 'value', 290, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% 
% One can then find better filter parameters without redetecting the
% surface by changing the second block of options in detectOptions and 
% calling resetMask and applyMasks. 

startmask = true(size(xp.detector.mask));

xp.detector.resetMask();
xp.detector.setManualMask(startmask);

myDetectOpts = struct(  'sigma', detectOptions.sigma, 'channels', detectOptions.channels,...
                        'zdir', detectOptions.zdir,...
                        'maxIthresh', 0.05, 'summedIthresh', 0,...
                        'sigZoutliers',1, 'scaleZoutliers', 80); 

xp.detector.setOptions(myDetectOpts);    
xp.detector.applyMasks();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [zab xp.stack.imageSize(3)],...
                                            'InitialMagnification', 40);

%% temporary chunck filter

badchunks = imclose(imfilter(xp.detector.surfaceMatrix,fspecial('log',9))<1,strel('disk',0));

% im = imfilter(xp.detector.surfaceMatrix,fspecial('log',6));
% % for i = 1%:2
% %     im = imdilate(im,strel('disk',13));
% % end
% imshow(stdfilt(im, ones(21)),[]);

CC = bwconncomp(badchunks);
stats = regionprops(CC,'area');
biggest = find([stats.Area] == max([stats.Area]));
badchunks(CC.PixelIdxList{biggest}) = false;
badchunks = imclose(badchunks,strel('disk',1));
badchunks = imdilate(badchunks,strel('disk',1));

figure, imshow(badchunks)

xp.detector.setManualMask(~badchunks);
xp.detector.applyMasks();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [zab xp.stack.imageSize(3)]);

%%
inspectOptions= struct('dimension', 'x', 'value', 300, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%%
% IMPORTANT VALUE TO PICK HERE

basalMaskParam = struct('mindepth', 310, 'closeSize', 10, 'openSize',50,...
                        'areaOpen', 1000);

% No.5 
% basalMaskParam = struct('mindepth', 333, 'closeSize', 20, 'openSize', 30)

% masking between the folds
bla = xp.detector.mask.*xp.detector.surfaceMatrix;
betweenfolds = (xp.detector.mask.*xp.detector.surfaceMatrix > basalMaskParam.mindepth);

% keep largest connected component
CC = bwconncomp(betweenfolds);
stats = regionprops(CC,'Area');
maxidx = find([stats.Area] == max([stats.Area]));
betweenfolds = false(size(betweenfolds));
betweenfolds(CC.PixelIdxList{maxidx}) = true;

imshow(betweenfolds)

% clean up
vla = imclose(betweenfolds, strel('disk',basalMaskParam.closeSize));
%vla = imfill(vla,'holes');
vla = ~bwareaopen(~vla, basalMaskParam.areaOpen);
vla = imopen(vla,strel('disk',basalMaskParam.openSize));

% biggest piece again
CC = bwconncomp(vla);
stats = regionprops(CC,'Area');
maxidx = find([stats.Area] == max([stats.Area]));
vla = false(size(vla));
vla(CC.PixelIdxList{maxidx}) = true;

imshow(vla)

%%
% % final outliers 1%
% blA = xp.detector.surfaceMatrix.*vla;
% imshow(blA,[300 xp.stack.imageSize(3)])
% [zdist,bins] = hist(blA(blA>0),100);
% zdist = zdist./sum(zdist);
% zcum = cumsum(zdist(:));
% zthresh = bins(find(zcum < 0.01,1,'last'))
% vla = vla & (blA > zthresh);

xp.detector.setManualMask(vla & ~badchunks);
xp.detector.applyMasks();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [zab xp.stack.imageSize(3)],...
                                            'InitialMagnification', 40);

                                        
%%
% We can also inspect a point cloud cross section over the data with
% detector.inspectQuality. In the pointCloud option, 'c' specifies the 
% color cyan.

inspectOptions= struct('dimension', 'x', 'value', 335, 'pointCloud', 'c');
xp.detector.inspectQuality(inspectOptions, xp.stack);

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_sectionx' num2str(inspectOptions.value) '.tif'];
saveas(gcf, fname);

%% save detection parameters and result

fname = fullfile(dataDir, ['analysis_' bareFnames{fi}], ['basalDetectOptions_' bareFnames{fi}]);
save(fname, 'myDetectOpts', 'basalMaskParam');

fname = fullfile(dataDir, ['analysis_' bareFnames{fi}], ['basalSurfaceMatrix_' bareFnames{fi} '.tif']);
imwrite(uint8(xp.detector.mask.*xp.detector.surfaceMatrix), fname);

%% 
% Or we can look at the point cloud in 3d, with some subsampling factor.
ssfactor = 50;
xp.detector.pointCloud.inspect(ssfactor);

%% Fit the surface for the disc proper cells
%
% By detecting the largest intensity jump along z for each x,y in the
% E-cad channel and filtering out local outliers we have found the apical
% surface of the disc proper cells. We can now fit a smooth surface
% representation to that.
%
% tpsFitter fits the pointcloud using a thin plate spline fit. It has the
% following options:
%
% * gridSize:     Size of grid on which to generate fitted surface
%               default [50 50], full size takes long.
% * smoothing:    TPS smoothing parameter (default 1000).

fitmask = imclose(vla,strel('disk',50));
fitOptions = struct('smoothing', 100, 'gridSize', [100 100],'fitMask',fitmask);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%xp.fitter.inspectTPS;

%%
% We can visualize the result on a cross section with
% fitter.inspectQuality.
figure,
inspectOptions= struct('dimension', 'x', 'value', 305, 'pointCloud', 'c');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_sectionx' num2str(inspectOptions.value) '.tif'];
saveas(gcf, fname);

%%
% We now generate the Surface Of Interest. The charts to be generated are 
% specified in xp.fitter.charts. In this case there is only one, called
% 'xy'. 

basalshift = -20;
fname = fullfile(dataDir, ['analysis_' bareFnames{fi}], 'basalshift');
save(fname,'basalshift');

xp.zEvolve(basalshift);
xp.generateSOI();

%% Pull back the data to the surface
% 
% We pull back the data to the SOI using pullbackStack.
onionOpts = struct('nLayers', 41, 'layerDistance', 1, 'sigma', 1,...
                    'makeIP', 'both', 'zEvolve', true);
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime);%onionOpts

%%
% To look at the pullback, we call the data field of the SOI at the right
% time, and get a particular patch from that with getPatch. A patch is a 
% part of a surface. In this case, there is only one called xy_index.
% Then we get the data in some patch in a particular coordinate system with
% getTransform. In this case there is only one coordinate system: xy.
% What we get is an object not only holding the image data but also
% metadata and methods to manipulate it. The actual data is obtained by
% calling the method apply. This returns a cell array with entries for each
% channel.

% xp.tIdx converts the time into an index in a list of time points
tidx = xp.tIdx(0);

discProperPatch = xp.SOI.data(tidx).getPatch('xy_index');
%discProperImage = discProperPatch.getTransform('xy').apply{channel};
pbs = discProperPatch.getTransform('xy').apply;
%pbs{1} = pbs{1}*0;
%pbs{2} = pbs{2}*0;

for i = 1:3, pbs{i} = imadjust(mat2gray(pbs{i})); end;
discProperImage = cat(3,pbs{[3 1 2]});
discProperImage = pbs{1};
figure,
imshow(discProperImage, [], 'InitialMagnification', 50);
fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_basalSurface.tif'];
imwrite(discProperImage, fname);

%% Save the result
%
% Finally we save the SOI using SOI.save. We set the following options:
%
% * dir:            The directory to save the SOI to.
% * imwriteOptions: Pullbacks are saved to image files using imwrite, we
% can pass options to change file format, compression etc. For example we
% could change this option to
% imwriteOptions = {'jp2', 'Mode', 'lossless'}; 
% * make8bit:       Often absolute intensities don't matter and 8 bit offers
% a large enough dynamic range. This options rescales the lookup table and
% converts to 8 bit before saving.

imwriteOptions = {'tif', 'Compression', 'deflate'};
savedir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_basalSOI']);

options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         true);
xp.SOI.save(options)

%% flat stacks

savedir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_basalSOI']);
xp.SOI.multilayer2stack(1, savedir); %dataDir);

%% 
%--------------------------------------------------------
% Get apical surface
%--------------------------------------------------------

% copy basal SOI before we overwrite it
savedir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_basalSOI']);
basalSOI = surfaceAnalysis.SurfaceOfInterest(savedir);

%%
zab = 160;

% change detection mask
projectionMask = cell([3 1]);
projectionMask{1} = true([xp.stack.imageSize(1) xp.stack.imageSize(3)]);
projectionMask{1}(:, zab:end) = false;
xp.stack.setProjectionMask(projectionMask);

imshow(xp.stack.getSlice('x',500));

%%
detectOptions = struct(  'sigma', 3.5, 'channels', 2, 'zdir', -3,...
                        'maxIthresh', 0.02, 'summedIthresh', 0,...
                        'sigZoutliers', 2, 'scaleZoutliers', 3); 

% Calling detectSurface runs the surface detector and creates the point
% cloud in detector.pointCloud.

xp.setDetectOptions(detectOptions);
xp.detectSurface();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [1 zab],...
                                            'InitialMagnification', 40);
                                        
%%
xp.detector.resetMask();
myDetectOpts = struct(  'sigma', detectOptions.sigma, 'channels', detectOptions.channels,...
                        'zdir', -3,...
                        'maxIthresh', 0.1, 'summedIthresh', 0,...
                        'sigZoutliers', 2, 'scaleZoutliers', 50); 

xp.detector.setOptions(myDetectOpts);    

% apply the basal mask too
apicalMask = fitmask;
xp.detector.setManualMask(apicalMask);
xp.detector.applyMasks();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [1 50],...
                                            'InitialMagnification', 40);
             
%% temporary chunck filter

badchunks = imclose(imfilter(xp.detector.surfaceMatrix,fspecial('log',11))<1,strel('disk',0));

% im = imfilter(xp.detector.surfaceMatrix,fspecial('log',6));
% % for i = 1%:2
% %     im = imdilate(im,strel('disk',13));
% % end
% imshow(stdfilt(im, ones(21)),[]);

CC = bwconncomp(badchunks);
stats = regionprops(CC,'area');
biggest = find([stats.Area] == max([stats.Area]));
badchunks(CC.PixelIdxList{biggest}) = false;
badchunks = imclose(badchunks, strel('disk',5));

%figure, imshow(~badchunks)

xp.detector.setManualMask(~badchunks.*fitmask);
xp.detector.applyMasks();

imshow(xp.detector.mask.*xp.detector.surfaceMatrix, [1 1]);

%%
% Looking at the fold-masked point cloud in a cross section we see that the
% mask works well.
inspectOptions= struct('dimension', 'x', 'value', 850, 'pointCloud', 'r');
xp.detector.inspectQuality(inspectOptions, xp.stack);

%% save detection parameters and result

fname = fullfile(dataDir, ['analysis_' bareFnames{fi}], ['apicalDetectOptions_' bareFnames{fi}]);
save(fname, 'myDetectOpts');

fname = fullfile(dataDir, ['analysis_' bareFnames{fi}], ['apicalSurfaceMatrix_' bareFnames{fi} '.tif']);
imwrite(uint8(xp.detector.mask.*xp.detector.surfaceMatrix), fname);

%%
% Fitting to this masked point cloud with we set smoothing higher because
% we are fitting a smoother surface (the peripodial cells don't fold).

fitOptions = struct('smoothing', 100, 'gridSize', [100 100],'fitMask',apicalMask);
xp.setFitOptions(fitOptions);
xp.fitSurface();

%%
apicalshift = 6;
xp.zEvolve(apicalshift);
xp.generateSOI();

%% inspect apical and basal surfaces together

projectionMask = cell([3 1]);
projectionMask{1} = true([xp.stack.imageSize(1) xp.stack.imageSize(3)]);
xp.stack.setProjectionMask(projectionMask);

h = figure,
BembGrids = basalSOI.embedding.patches{1}.apply;
inspectOptions= struct('dimension', 'x', 'value', 730, 'pointCloud', 'c');
xp.fitter.inspectQuality(inspectOptions, xp.detector, xp.stack);
hold on
plot(squeeze(BembGrids{3}(:,inspectOptions.value,:)),'-r','LineWidth',2)
hold off
%%
fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_surfaceXsect' num2str(inspectOptions.value) '.tif'];
saveas(h,fname);

%%
fname = fullfile(dataDir, ['analysis_' bareFnames{fi}], 'apicalshift');
save(fname,'apicalshift');

%% Pull back the data to the surface
% 
% We pull back the data to the SOI using pullbackStack.
onionOpts = struct('nLayers', abs(apicalshift)*2-3, 'layerDistance', 1, 'sigma', 1, 'makeIP', 'both', 'zEvolve', true);
xp.SOI.pullbackStack(xp.stack, xp.currentROI, xp.currentTime, onionOpts);

%%
t=1;
discProperPatch = xp.SOI.data(t).getPatch('xy_index');
%discProperImage = discProperPatch.getTransform('xy').apply{channel};
pbs = discProperPatch.getTransform('xy').apply;
for i = 1:3, pbs{i} = mat2gray(pbs{i}); end;
[~,colorInv] = sort(xp.expMeta.channelColor);
discProperImage = cat(3,pbs{colorInv});

h = figure;
imshow(discProperImage, [], 'InitialMagnification', 50);

fname = [fullfile(dataDir, ['analysis_' bareFnames{fi}], bareFnames{fi}), '_apicalSurface.tif'];
saveas(h, fname);

%%
% Finally we save it to a different directory from before, because it is a
% different surface.

imwriteOptions = {'tif', 'Compression', 'deflate'};
savedir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_apicalSOI']);

options = struct(   'dir',              savedir,...
                    'imwriteOptions',   {imwriteOptions},...
                    'make8bit',         true);
xp.SOI.save(options)

%% flat stacks

savedir = fullfile(dataDir, ['analysis_' bareFnames{fi}], [bareFnames{fi} '_apicalSOI']);
xp.SOI.multilayer2stack(1, savedir); %dataDir);

