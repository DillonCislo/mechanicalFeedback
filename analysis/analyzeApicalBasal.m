clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));
addpath(genpath('/Users/idse/repos/epitheliumAnalysis'));

% import data locations file
dataLocationsYW20150903;
results = {};

biRange = [1 7 12:14];%1:5%1:9%; %bi = 1;%7

% bi : batch index
% fi : file index

%%
for bi = biRange

bi
masterDataDir = dataDir{bi};

for fi = 1:numel(resolutions{bi});

fi
bareFname = bareFnames{bi}{fi};

disp(description{bi})
disp(bareFname)

% we don't store channelsUsed so here I assume it's everything except DNA 
setlabels = labels{bi};
setlabels(strcmp(setlabels,'DNA')) = [];
cloneC = find(strcmp(setlabels,'clone'));
memC = find(strcmp(setlabels,'arm'));
specialC = find(~(strcmp(setlabels,'arm') | strcmp(setlabels,'clone')));

resultsDir = fullfile(masterDataDir, ['analysis_' bareFname]);
if ~exist(fullfile(resultsDir, 'results'),'dir')
    mkdir(fullfile(resultsDir, 'results'));
end

t = 1;

%-------------------
% load surfaces
%-------------------

% load apical surface

SOIdir = fullfile(resultsDir, [bareFname '_apicalSOI']);
apicalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

fname = [fullfile(resultsDir,'tripleLayerSeg', bareFname), '_cloneMaskApical.tif'];
cloneMaskApical = imread(fname);
cloneBdryApical = bwboundaries(cloneMaskApical');

fname = [fullfile(resultsDir, 'tripleLayerSeg', bareFname), '_apicalSeg.tif'];
apicalSeg = imread(fname);

fname = fullfile(resultsDir, 'tripleLayerSeg', 'apicalVertexModel');
S = load(fname);
cellLayer = S.cellLayer;
clear S;

% load the basal surface

fname = [fullfile(resultsDir, 'tripleLayerSeg', bareFname), '_apicalL.tif'];

% existence of clone matching is a proxy for existence of basal surface
if exist(fullfile(resultsDir, [bareFname '_basalSOI']),'dir')
    
    hasbasal = true;
    
    cloneApicalL = imread(fname);

    SOIdir = fullfile(resultsDir, [bareFname '_basalSOI']);
    basalSOI = surfaceAnalysis.SurfaceOfInterest(SOIdir);

    fname = [fullfile(resultsDir, 'tripleLayerSeg', bareFname), '_cloneMaskBasal.tif'];
    cloneMaskBasal = imread(fname);
    cloneBdryBasal = bwboundaries(cloneMaskBasal');

    fname = [fullfile(resultsDir, 'tripleLayerSeg', bareFname), '_basalL.tif'];
    cloneBasalL = imread(fname);
    
else
    
    hasbasal = false;
    basalSOI = [];
    cloneBasalL = [];
    
    % select special clones apically
    bw1 = bwlabel(cloneMaskApical);
    newlabels1 = zeros(size(cloneMaskApical),'uint16');
    label = 1;
    clickedobject = 1;
    
    if ~exist(fname, 'file')
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
        close
    else
        cloneApicalL = imread(fname);
    end
end

wing = struct('apicalSOI', apicalSOI, 'basalSOI', basalSOI,...
              'cloneApicalL', cloneApicalL, 'cloneBasalL', cloneBasalL,...
              'apicalSegmentation', cellLayer, 'resolution', resolutions{bi}{fi},...
              'channelLabels', {setlabels}, 'specialChannel', specialC,...
              'rawDataName', bareFname);

% % add proper areas to wing structure
% 
% % proper cell area
% if isempty(wing.apicalSOI.g.patches)
%     wing.apicalSOI.NCalcInducedMetric(t);
% end
% 
% % metric and metric determinant
% tic
% chartName = 'xy';
% chart = wing.apicalSOI.atlas.getChart(chartName);
% gpatch = wing.apicalSOI.g.getPatch(chart.domain.name);
% metric = gpatch.getTransform(chartName).apply;
% 
% detg = gpatch.determinant.apply{1};
% dA = sqrt(detg)*chart.image.stepSize(1)*chart.image.stepSize(2);
% 
% assert(gpatch ~=0, 'calculate metric first, SOI.calcInducedMetric');
% 
% Acorrected = [];
% for ci = 1:length(wing.apicalSegmentation.cells{t})
%     l = wing.apicalSegmentation.cells{t}(ci).label;
%     cidx = wing.apicalSegmentation.L == l;
%     Acorrected(ci) = sum(dA(cidx)); %A = sum(cidx(:))
% end
% toc
% 
% wing.properAreas = Acorrected;

%-----------------
% non-autonomy
%-----------------

t = 1;
bintypes = {'boundary'};%,'center'};

for clIdx = 1:max(cloneApicalL(:))
    
    for i = 1:numel(bintypes)

        type = bintypes{i};

        scrsz = get(0,'ScreenSize');
        h = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
        results{bi}{fi}{clIdx}{i} = measureCloneAutonomy(type, t, clIdx, wing);
        set(gcf,'Color','white');
%         if ~isempty(results)
%             I = getframe(gcf);
%             fname = [fullfile(resultsDir, 'results', bareFname), '_' type '_nonautonomy_' num2str(clIdx) '.tif'];
%             imwrite(I.cdata, fname);
%         end
        close
    end
end

end

end

%%
%----------------------------------
% scatter plot of many clones
%----------------------------------

bintypes = {'boundary distance binning'};%, 'center distance binning'};
orientationtypes = {'normal','radial'};

for binTypeIndex = 1:numel(bintypes)
for typeIndex = 2%1:numel(orientationtypes)

clCenter = [];
clArea = [];
orientation =[];
N = 0;
    
for bi = biRange%1:9%; %bi = 1;%7
    
    for fi = 1:numel(resolutions{bi})
        
        for i = 1:numel(results{bi}{fi})
                
            if ~isempty(results{bi}{fi}{i}{1})
                
                N = N+1;

                clCenter(N,:) = results{bi}{fi}{i}{binTypeIndex}.clCenter;
                clArea(N) = results{bi}{fi}{i}{binTypeIndex}.clArea;
                if typeIndex == 1
                    orientation(N,:) = results{bi}{fi}{i}{binTypeIndex}.orientationRadial(:);
                elseif typeIndex == 2
                    orientation(N,:) = results{bi}{fi}{i}{binTypeIndex}.orientationNormal(:);
                end
            end
        end
    end
end

centerDist = sqrt(sum((clCenter - 512).^2,2));
posMedian = median(centerDist);

figure,
fsize = 16;
idx = centerDist < posMedian;
colors = sqrt(clArea(idx));
prop = orientation;
scatter(prop(idx,1), prop(idx,3), 30, colors, 'filled');
hold on
idx = centerDist >= posMedian;
colors = sqrt(clArea(idx));
scatter(prop(idx,1), prop(idx,3), 40, colors, 'square');
hold off
xlabel('orientation on clone edge','FontSize',fsize);
ylabel('orientation 15 micron away','FontSize',fsize);
%graphtitle = [['bantam, weighted ' orientationtypes{typeIndex} ' orientation '] {[bintypes{binTypeIndex} ' 15']}];
graphtitle = {'anisotropy orientation near bantam clones'};
title(graphtitle,'FontSize',fsize);
axis equal
axis([-1 1 -1 1]);
fname = fullfile(mainDataDir, [graphtitle{:} '.png']);
saveas(gcf, fname);
set(gca,'LineWidth',2);
set(gca,'Fontsize',12);

end
end

%% special channel with clone outline

lw = 2;
field = apicalSOI.getField('data');
pbs = field(t).getPatch('xy_index').apply;
for i = 1:3, pbs{i} = mat2gray(pbs{i}); end;
apicalPB = cat(3,pbs{[3 1 2]});
apicalPB = mat2gray(pbs{1});

h = figure;
imshow(apicalPB, [], 'InitialMagnification', 50);
hold on
for i = 1:numel(cloneBdryApical)
    plot(cloneBdryApical{i}(:,1),cloneBdryApical{i}(:,2), '-r','LineWidth',lw)
end
if hasbasal
    for i = 1:numel(cloneBdryBasal)
        plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2), '-b','LineWidth',lw)
    end
end
hold off

fname = [fullfile(resultsDir, 'results', bareFname), '_apical_' setlabels{specialC} '.tif'];
saveas(h, fname);

%% height map with clone outline

if hasbasal
    basalZ = basalSOI.embedding.patches{1}.apply{3};
    apicalZ = apicalSOI.embedding.patches{1}.apply{3};
    height = basalZ - apicalZ;
else
    height = apicalSOI.embedding.patches{1}.apply{3};
end
height = double(height)*resolutions{bi}{fi}; % convert to microns
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
end
if hasbasal
    for i = 1:numel(cloneBdryBasal)
        plot(cloneBdryBasal{i}(:,1),cloneBdryBasal{i}(:,2), '-k', 'LineWidth',1)
    end
end
hold off

fname = [fullfile(resultsDir, 'results', bareFname), '_heightMap.tif'];
saveas(h, fname);


%% apical and basal density in clone
t = 1;

stats = {};

statsA = regionprops(cloneApicalL,'Area','Centroid');
cellArea = cellLayer.getCellGeometry(t,'area');
areaBins = 0:15:450;

if isempty(apicalSOI.g.patches)
    apicalSOI.NCalcInducedMetric;
    detgA = apicalSOI.g.patches{1}.determinant.apply{1};
end

if hasbasal
    statsB = regionprops(cloneBasalL,'Area','Centroid');

    if isempty(basalSOI.g.patches)
        basalSOI.NCalcInducedMetric;
        detgB = basalSOI.g.patches{1}.determinant.apply{1};
    end
end

nMatchedClones = length(statsA);
cloneLabels = [];

% for special distribution
field = apicalSOI.getField('data_SIP');
specialIm = field.patches{1}.apply{specialC};

for i = 1:nMatchedClones
    
    stats{i}.apicalArea = statsA(i).Area;
    stats{i}.apicalAreaProper = sum(detgA(cloneApicalL == i));
    
    % backwards way of getting the label but whatever
    I = round(statsA(i).Centroid);
    CM = cellLayer.getCellGeometry(t,'centroid')';
    CM(cat(1,cellLayer.cells{t}.state) == 0,:) = 0;
    
    [~,clonecentercellidx] = min(sqrt((CM(:,1)-I(1)).^2 + (CM(:,2)-I(2)).^2));
    cloneLabel = cellLayer.cells{t}(clonecentercellidx).state;
    cloneCellIdx = cat(1,cellLayer.cells{t}.state) == cloneLabel;
    stats{i}.cloneLabel = cloneLabel;
    
    stats{i}.apicalCellCount = sum(cloneCellIdx);
    stats{i}.apicalAreaSegmented = sum(cellArea(cloneCellIdx));
    stats{i}.apicalAreaDist = makeDist(cellArea(cloneCellIdx), areaBins);
    stats{i}.apicalDensity = stats{i}.apicalCellCount./(stats{i}.apicalArea*resolutions{bi}{fi}^2);
    
    stats{i}.special = setlabels{specialC};
    if i == 1
        stats{i}.specialDist = makeDist(specialIm(cloneApicalL == i), 30);
    else
        stats{i}.specialDist = makeDist(specialIm(cloneApicalL == i), stats{1}.specialDist.bins);
    end
    
    if hasbasal
        stats{i}.basalArea = statsB(i).Area;
        stats{i}.basalAreaProper = sum(detgB(cloneBasalL == i));
        stats{i}.basalDensity = stats{i}.apicalCellCount./(stats{i}.basalArea*resolutions{bi}{fi}^2);    

        stats{i}.areaRatio = statsA(i).Area/statsB(i).Area;
        stats{i}.areaRatioProper = stats{i}.apicalAreaProper/stats{i}.basalAreaProper;
    end
end

inclone = cat(1,cellLayer.cells{t}.state)*0;
for i = 1:numel(stats)
    inclone = inclone | cat(1,cellLayer.cells{t}.state) == stats{i}.cloneLabel;
end

%% visualize

figure, 
imshow(apicalSOI.data.patches{1}.apply{specialC},[]);

areas = cellLayer.getCellGeometry(t,'area');
bdryBondIdx = cellLayer.getStateBoundary(t,'clone',1);
cloneBdryBonds = cellLayer.bonds{t}(bdryBondIdx(:,1));

hold on 
options = struct('colorTable', inclone, 'transparent', 0);
cellLayer.visualize(t, options);

for i = 1:numel(cloneBdryBonds)
    V = cellLayer.vertices{t}(cloneBdryBonds(i).vertInd,1:2);
    plot(V(:,1),V(:,2),'-k','LineWidth',3);
end
hold off
axis equal

% fname = [fullfile(tldir, bareFname), '_areaMap.tif'];
% saveas(gcf, fname);

%% make fake clones

figure,
imshow(cat(3, mat2gray(apicalSOI.data.patches{1}.apply{specialC}),cloneMaskApical,cloneMaskApical));

fakeclones = false(size(cloneMaskApical));
for i = 1:2
   BW = roipoly; 
   fakeclones = fakeclones | BW;
end

%%
cellCM = cellLayer.getCellGeometry(t, 'centroid');
CMsub = round(cellCM)';
fakecloneL = bwlabel(fakeclones);
linind = sub2ind(size(cloneMaskApical), CMsub(:,2), CMsub(:,1));
%inclone = cloneMaskApical(linind);
infakeclone = uint16(fakecloneL(linind));

statsA = regionprops(fakecloneL,'Area','Centroid');
nFakeClones = length(statsA);

fakestats = {};

for i = 1:nFakeClones
    
    fakestats{i}.mask = fakecloneL==i;

    fakestats{i}.apicalArea = statsA(i).Area;
    fakestats{i}.apicalAreaProper = sum(detgA(fakestats{i}.mask));
    
    cloneCellIdx = infakeclone==i;    
    fakestats{i}.apicalCellCount = sum(cloneCellIdx);
    fakestats{i}.apicalAreaSegmented = sum(cellArea(cloneCellIdx));
    fakestats{i}.apicalAreaDist = makeDist(cellArea(cloneCellIdx), areaBins);
    
    fakestats{i}.apicalDensity = fakestats{i}.apicalCellCount./(fakestats{i}.apicalArea*resolutions{bi}{fi}^2);
    fakestats{i}.specialDist = makeDist(specialIm(fakestats{i}.mask), stats{1}.specialDist.bins);
    
    if hasbasal
        fakestats{i}.basalAreaProper = sum(detgB(fakestats{i}.mask));
        fakestats{i}.basalDensity = fakestats{i}.apicalCellCount./(fakestats{i}.apicalArea*resolutions{bi}{fi}^2);
        fakestats{i}.areaRatioProper = fakestats{i}.apicalAreaProper/fakestats{i}.basalAreaProper;
    end
end

%% numbers in txtfile

fname = [fullfile(resultsDir, 'results', bareFname), '_densitystats.mat'];
save(fname, 'stats','fakestats');

txtfile = [fullfile(resultsDir, 'results', bareFname), '_stats.txt'];
if exist(txtfile)
    delete(txtfile);
end
diary(txtfile)

disp('------------------ clone stats ------------------');

for i = 1:nMatchedClones
    stats{i}
end

disp('------------------ fake clone stats ------------------');

for i = 1:nFakeClones
    fakestats{i}
end

diary off
    
%% fake clone outlines in image

fakeBdry = bwboundaries(fakeclones');
matchedApicalBdry = bwboundaries((cloneApicalL > 0)');

h = figure;
imshow(apicalSOI.data.patches{1}.apply{1}, [], 'InitialMagnification', 50);
hold on
for i = 1:numel(matchedApicalBdry)
    plot(matchedApicalBdry{i}(:,1),matchedApicalBdry{i}(:,2), '-r','LineWidth',lw)
end
for i = 1:numel(fakeBdry)
    plot(fakeBdry{i}(:,1),fakeBdry{i}(:,2), '-c','LineWidth',lw)
end
if hasbasal
    matchedBasalBdry = bwboundaries((cloneBasalL > 0)');
    for i = 1:numel(matchedBasalBdry)
        plot(matchedBasalBdry{i}(:,1), matchedBasalBdry{i}(:,2), '-b','LineWidth',lw)
    end
end
hold off

fname = [fullfile(resultsDir, 'results', bareFname), '_fakeclones.tif'];
saveas(h, fname);


%% plot distributions 

figure,
cloneMask = cloneApicalL > 0;

% area stats
statsCombined = {};
cellArea = cellLayer.getCellGeometry(t, 'area');
statsCombined{1}.area = makeDist(cellArea(inclone > 0), areaBins);
statsCombined{2}.area = makeDist(cellArea(infakeclone > 0), areaBins);

row = 1;
nrows = 2;
ncols = 2;
plotInfo = struct('plotTitle', 'cell area');
fieldname = 'area';
subplotDist(nrows, ncols, row, plotInfo, statsCombined, fieldname, {'clone','WT'})
                
% special stats
row = 2;
statsCombined{1}.specialDist = makeDist(specialIm(cloneMask), stats{1}.specialDist.bins);
statsCombined{2}.specialDist = makeDist(specialIm(fakeclones), stats{1}.specialDist.bins);

plotInfo = struct('plotTitle', [setlabels{specialC} ' intensity']);
subplotDist(nrows, ncols, row, plotInfo, statsCombined, 'specialDist', {'clone','WT'})
  
fname = [fullfile(resultsDir, 'results', bareFname), '_' setlabels{specialC} '_combinedStats.png'];
saveas(gcf, fname);

%% plot individual distributions

allstats = [stats, fakestats];

figure,

row = 1;
nrows = 2;
ncols = 2;
plotInfo = struct('plotTitle', 'cell area');
fieldname = 'apicalAreaDist';

% plotInfo = struct('plotTitle', plotTitles{row});
subplotDist(nrows, ncols, row, plotInfo, allstats, fieldname, {},...
    {1:numel(stats), (numel(stats)+1) : (numel(stats)+numel(fakestats))},...
                {'clone','WT'});

row = 2;
fieldName = 'specialDist';
plotInfo = struct('plotTitle', [setlabels{specialC} ' intensity']);
subplotDist(nrows, ncols, row, plotInfo, allstats, fieldName, {},...
    {1:numel(stats), (numel(stats)+1) : (numel(stats)+numel(fakestats))},...
                {'clone','WT'});

fname = [fullfile(resultsDir, 'results', bareFname), '_' setlabels{specialC} '_individualStats.png'];
saveas(gcf, fname);
