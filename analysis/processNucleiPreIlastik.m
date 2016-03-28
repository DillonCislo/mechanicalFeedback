clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));

% import data locations file
dataLocations; 

bi = 10;
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
fname = [fullfile(dataDir, bareFnames{i}), '_scaleSS.tif'];

disp(bareFnames{i});

info = imfinfo(fname);
nSlices = numel(info)/3;
im = zeros([info(1).Height, info(1).Width nSlices], 'uint16');

for zi = 1:nSlices
    im(:,:,zi) = imread(fname, (zi-1)*3 + 2);
end

%% try enhancing nuclei 
% DOESNT REALLY SEEM TO WORK

xsect = squeeze(im(400,:,:))';
imshow(imadjust(xsect),[])

%%
xsect = squeeze(im(400,:,:))';
xsect = mat2gray(xsect(1:300,1:400));

sigxy = 10;
sigz = 1.5*sigxy;

ker = -GaussD(sigxy, 2, 1);

logstack = LOG(xsect, [sigz sigxy]);
imshow(logstack,[]);

%%
[X,Y] = meshgrid(-20:20, -30:30);
R = 12;
bla = X.^2 + Y.^2/5 < R.^2;

xsect = squeeze(im(400,:,:))';
xsect = mat2gray(xsect(1:300,1:400));
hatx = xsect - imopen(xsect, bla);
imshow(cat(3, imadjust(xsect), imadjust(hatx), imadjust(xsect-2*hatx)),[]);








