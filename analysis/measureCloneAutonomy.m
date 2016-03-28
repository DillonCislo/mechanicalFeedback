function results = measureCloneAutonomy(type, t, clIdx, wing)

    clMask = wing.cloneApicalL == clIdx;
    clMask = imfill(clMask,'holes');
    
    % better than just apicalPB > 0 to mask out disc folds / edges
    discMask = imfill(imopen(wing.apicalSegmentation.L > 0,strel('disk',1)),'holes');
    
    cellLayer = wing.apicalSegmentation;
    res = wing.resolution;
    specialC = wing.specialChannel;
    setlabels = wing.channelLabels;

    s = regionprops(clMask, 'centroid',...
                                'EquivDiameter','BoundingBox');
    cloneBdryApical = bwboundaries(clMask');
    cloneBdryApical = cloneBdryApical{1};
    
    if numel(s) ~= 1
        warning('disconnected clone: returning');
        results = [];
        return;
    end
    
    ysize = size(wing.cloneApicalL,1);
    xsize = size(wing.cloneApicalL,2);
    
    hasbasal = ~isempty(wing.basalSOI);
    
    % distance from clone boundary
    insideDist = bwdist(~clMask);
    bdryDist = bwdist(clMask) - (insideDist - 1);

    % imshow(bdryDist,[])
    % hold on
    % plot(cloneBdryApical{clIdx}(:,1), cloneBdryApical{clIdx}(:,2), '-r','LineWidth',lw)
    % hold off
    
    % distance from clone center
    clCenter = s.Centroid;
    equivRadius = s.EquivDiameter*res/2;
    [X,Y] = meshgrid(1:xsize, 1:ysize);
    [~,centerDist] = cart2pol(X - clCenter(1), Y - clCenter(2));
    equivCircle = bwboundaries((centerDist < equivRadius/res)');
    equivCircle = equivCircle{1};

    % imshow(centerDist,[])
    % hold on
    % plot(cloneBdryApical{clIdx}(:,1), cloneBdryApical{clIdx}(:,2), '-r','LineWidth',lw)
    % hold off
    
    if strcmp(type,'boundary')
        
        distmap  = bdryDist;
        % bin edges
        binE = linspace(-150,200,31); % min(distmap(:))
        % bin center
        binC = conv(binE, [1 1]/2, 'valid');
        edgePos = 0;
        distlabel = 'distance from clone edge (micron)';
        
    elseif strcmp(type,'center');
        
        %bins = 300*sqrt((bins/300));
        distmap  = centerDist;
        binE = linspace(0,350,31); %min(distmap(:))
        % bin center
        binC = conv(binE, [1 1]/2, 'valid');
        edgePos = equivRadius;
        distlabel = 'distance from clone center (micron)';
    else
        error('type unknown');
    end
    
    % apical PB for myo / jub / etc levels
    field = wing.apicalSOI.getField('data');
    pbs = field(t).getPatch('xy_index').apply;
    apicalPB = mat2gray(pbs{1});
    
    % height
    apicalZ = wing.apicalSOI.embedding.patches{1}.apply{3};
    if hasbasal
        basalZ = wing.basalSOI.embedding.patches{1}.apply{3};
        height = basalZ - apicalZ;
    else
        height = wing.apicalSOI.embedding.patches{1}.apply{3};
    end
    height = double(height)*res;
    
    % cell centers for density and masking
    CM = cellLayer.getCellGeometry(t, 'centroid')';
    CMsub = round(CM);
    CMind = sub2ind([ysize, xsize], CMsub(:,2), CMsub(:,1));

    % cell area for alternative density 
    areas = round(cellLayer.getCellGeometry(t, 'area')');
    areas = areas*res^2;

    % anisotropy
    lmaj = round(cellLayer.getCellGeometry(t, 'majorAxisLength')');
    lmin = round(cellLayer.getCellGeometry(t, 'minorAxisLength')');
    a = ((lmaj - lmin)./(lmaj + lmin))';

    % cell orientation
    %------------------
    
    % cell major axis
    phi = cellLayer.getCellGeometry(t, 'orientation');
    cellAxis = [cosd(phi), sind(phi)];

    % normal field to the clone boundary
    sigma = 5;
    [dx,dy] = GaussGradient(bdryDist, sigma);
    N = sqrt(dx.^2 + dy.^2);
    clNormal = [dx(CMind)./N(CMind), -dy(CMind)./N(CMind)];
    
    [dx,dy] = GaussGradient(centerDist, sigma);
    N = sqrt(dx.^2 + dy.^2);
    clRadial = [dx(CMind)./N(CMind), -dy(CMind)./N(CMind)];
%     phiCl = atan2d(clNormal(:,2),clNormal(:,1));
%     phiCl(phiCl < 0) = 180 + phiCl(phiCl < 0);

    % cell orientation relative to clone boundary
    %cosphirel = sum(clNormal.*cellAxis, 2);
    cosphiNormal = sum(clNormal(:,1).*cellAxis(:,2) - clNormal(:,2).*cellAxis(:,1), 2);
    nemOrderNormal = 2*cosphiNormal.^2-1;

    W = a.^2;%/mean(a.^2);
    cosphiNormalWeight = sum(clNormal(:,1).*cellAxis(:,2) - clNormal(:,2).*cellAxis(:,1), 2);
    nemOrderNormalWeight = W.*(2*cosphiNormalWeight.^2-1);
    
    cosphiRadial = sum(clRadial(:,1).*cellAxis(:,2) - clRadial(:,2).*cellAxis(:,1), 2);
    nemOrderRadial = 2*cosphiRadial.^2-1;
    
    cosphiRadialWeight = sum(clRadial(:,1).*cellAxis(:,2) - clRadial(:,2).*cellAxis(:,1), 2);
    nemOrderRadialWeight = W.*(2*cosphiRadialWeight.^2-1);
    
% %     %imshow(bdryDist,[]);
% figure
% 	imshow(apicalSOI.data.patches{1}.apply{1},[]);
% %     imshow(cellLayer.cellL(t) == 1398,[])
% %     options = struct('cellIndex', false, 'transparent', 'all');%'colorTable',phirel);
% %     %cellLayer.visualize(t, options);
%     hold on
%     quiver(X(CMind),Y(CMind),clNormal(:,1), -clNormal(:,2),1,'g');
%     scale = 10;
%     vx = scale*cellAxis(:,1);
%     vy = -scale*cellAxis(:,2);
%     quiver(CM(:,1)-vx/2,CM(:,2)-vy/2,vx,vy,0,'.b');
%     %scatter(CM(:,1),CM(:,2),'.b');
%     lw = 2;
%     plot(cloneBdryApical(:,1), cloneBdryApical(:,2), '-k','LineWidth',lw);
%     hold off

    % for orientation scatter plot 
    
    edgeMarg = round(1/res);
    d = round(10/res);
    
    clEdgeMask = distmap > edgePos/res - edgeMarg  & distmap < edgePos/res + edgeMarg;
    clEdgeMask = clEdgeMask & discMask & apicalPB > 0;
    
    clAwayMask = distmap > edgePos/res + d - edgeMarg  & distmap < edgePos/res + d + edgeMarg;
    clAwayMask = clAwayMask & discMask & apicalPB > 0;
    
    clAway15Mask = distmap > edgePos/res + round(15/res) - edgeMarg  & distmap < edgePos/res + round(15/res) + edgeMarg;
    clAway15Mask = clAway15Mask & discMask & apicalPB > 0;
    
    edgeNemNormW = mean(nemOrderNormalWeight(clEdgeMask(CMind)))/mean(W(clEdgeMask(CMind)));
    awayNemNormW = mean(nemOrderNormalWeight(clAwayMask(CMind)))/mean(W(clAwayMask(CMind)));
    edgeNemNormWstd = std(nemOrderNormalWeight(clEdgeMask(CMind)))/mean(W(clEdgeMask(CMind)));
    awayNemNormWstd = std(nemOrderNormalWeight(clAwayMask(CMind)))/mean(W(clAwayMask(CMind)));
    
    edgeNemRadialW = mean(nemOrderRadialWeight(clEdgeMask(CMind)))/mean(W(clEdgeMask(CMind)));
    awayNemRadialW = mean(nemOrderRadialWeight(clAwayMask(CMind)))/mean(W(clAwayMask(CMind)));
    edgeNemRadialWstd = std(nemOrderRadialWeight(clEdgeMask(CMind)))/mean(W(clEdgeMask(CMind)));
    awayNemRadialWstd = std(nemOrderRadialWeight(clAwayMask(CMind)))/mean(W(clAwayMask(CMind)));
    
    away15NemNormW = mean(nemOrderNormalWeight(clAway15Mask(CMind)))/mean(W(clAway15Mask(CMind)));
    away15NemRadialW = mean(nemOrderRadialWeight(clAway15Mask(CMind)))/mean(W(clAway15Mask(CMind)));
    
    % bounding box
    ymargin  = 100;
    ybmin = max(1,round(s.BoundingBox(2)) - ymargin);
    ybmax = min(ysize,round(s.BoundingBox(2) + s.BoundingBox(4)) + ymargin);
    aspect = 1.5;
    xbmin = round(s.BoundingBox(1));
    xbmax = round(s.BoundingBox(1) + s.BoundingBox(3));
    xmargin = ((ybmax-ybmin)*aspect - (xbmax-xbmin))/2;
    xbmin = max(1, xbmin - xmargin);
    xbmax = min(xsize, xbmax + xmargin);
    bbox = struct('xbmin',xbmin,'xbmax',xbmax,'ybmin',ybmin,'ybmax',ybmax);

    % histogram
    meanI = [];                 stdI = [];
    density = [];
    density2 = [];
    meanAnisotropy = [];        stdAnisotropy = [];
    meanH = [];                 stdH = [];
    meanApicalZ = [];
    meanAreas = [];
    meanAreasProper = [];
    meanNemNormal = [];         stdNemNormal = [];     
    meanNemRadial = [];         stdNemRadial = [];  
    meanNemNormalW = [];        stdNemNormalW = [];     
    meanNemRadialW = [];        stdNemRadialW = [];  
    
    for i = 1:numel(binC);

        mask = distmap > binE(i) & distmap <= binE(i+1);
        %mask(470:end,:) = false;

        smask = mask & discMask & apicalPB > 0;
        meanI(i) = mean(apicalPB(smask));
        stdI(i) = std(apicalPB(smask));
        
        % cell properties
        CMindMask = smask(CMind);
        
        density(i) = sum(CMindMask)./(sum(smask(:))*res^2);
        density(density == 0) = NaN;

        meanAreas(i) = mean(areas(CMindMask));
        density2(i) = mean(1./areas(CMindMask));

        if isfield(wing, 'properAreas')
            Aprop = wing.properAreas*res^2;
            meanAreasProper(i) = mean(Aprop(CMindMask)); 
        end
        
        meanAnisotropy(i) = mean(a(CMindMask));
        stdAnisotropy(i) = std(a(CMindMask));

        % orientation
        meanNemNormal(i) = mean(nemOrderNormal(CMindMask));
        stdNemNormal(i) = std(nemOrderNormal(CMindMask));
        meanNemRadial(i) = mean(nemOrderRadial(CMindMask));
        stdNemRadial(i) = std(nemOrderRadial(CMindMask));
        
        % anisotropy weighted orientation
        meanNemNormalW(i) = mean(nemOrderNormalWeight(CMindMask))/mean(W(CMindMask));
        stdNemNormalW(i) = std(nemOrderNormalWeight(CMindMask))/mean(W(CMindMask));
        meanNemRadialW(i) = mean(nemOrderRadialWeight(CMindMask))/mean(W(CMindMask));
        stdNemRadialW(i) = std(nemOrderRadialWeight(CMindMask))/mean(W(CMindMask));
        
        % height 
        smask = mask & discMask & apicalZ > 0;
        meanApicalZ(i) = mean(apicalZ(smask));

        smask = mask & (height > 0);
        meanH(i) = mean(height(smask));
        stdH(i) = std(height(smask));

    end
    
    % structure results
    
    results = struct('meanI', meanI, 'density', density,...
                     'areas', meanAreas, 'areasProper', meanAreasProper,...
                     'anisotropy', meanAnisotropy,...
                     'meanApicalZ', meanApicalZ, 'meanH', meanH,...
                     'boundingBox',bbox, 'orientationNormal',...
                     [  edgeNemNormW, edgeNemNormWstd;...
                        awayNemNormW, awayNemNormWstd;...
                        away15NemNormW, 0],...
                     'orientationRadial',...
                     [  edgeNemRadialW, edgeNemRadialWstd;....
                        awayNemRadialW, awayNemRadialWstd;...
                        away15NemRadialW, 0],...
                     'clCenter', s.Centroid, 'clArea', sum(clMask(:)));
    
	%------------------------------------------------------------
	% visualize
    %------------------------------------------------------------
    
    nrows = 3; ncols = 3; ploti = 1;
    titlefsize = 13;
    
    % myosin image
    subplot(nrows,ncols,ploti)
    ploti = ploti + 1;
    
    otherClones = wing.cloneApicalL > 0;
    otherClones(clMask) = 0;
    otherClones = otherClones(ybmin:ybmax, xbmin:xbmax);
    apicalPBcl = apicalPB(ybmin:ybmax, xbmin:xbmax);
    
    imshow(cat(3,apicalPBcl,apicalPBcl,apicalPBcl + 0.5*otherClones));
    
    hold on
    lw = 1;
    plot(cloneBdryApical(:,1) - xbmin,...
        cloneBdryApical(:,2) - ybmin, '-r','LineWidth',lw)
    if strcmp(type, 'center')
        plot(equivCircle(:,1)- xbmin, equivCircle(:,2)- ybmin, '-m','LineWidth',lw)
    end
    hold off
    title(setlabels{specialC}, 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel(wing.rawDataName, 'Interpreter','none');

    % height image
    subplot(nrows,ncols,ploti)
    ploti = ploti + 1;
    zmin = min(height(:));
    zmax = max(height(:));
    imshow(height, [zmin zmax]);
    shading interp
    hold on
    lw = 2;
    plot(cloneBdryApical(:,1),...
        cloneBdryApical(:,2), '-k','LineWidth',lw)
    if strcmp(type, 'center')
        plot(equivCircle(:,1), equivCircle(:,2), '-m','LineWidth',lw)
    end
    hold off
    title('height', 'FontSize', titlefsize, 'FontWeight', 'bold');
    axis([xbmin xbmax ybmin ybmax]);
    
    % anisotropy image
    subplot(nrows,ncols,ploti)
    ploti = ploti + 1;
    options = struct('colorTable', a);
    imshow(apicalPB);
    hold on
    lw = 2;
    cellLayer.visualize(t, options);
    plot(cloneBdryApical(:,1),...
        cloneBdryApical(:,2), '-r','LineWidth',lw)
    if strcmp(type, 'center')
        plot(equivCircle(:,1), equivCircle(:,2), '-m','LineWidth',lw)
    end
    hold off
    title('segmentation', 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel('color: anisotropy', 'Interpreter','none');
    axis([xbmin xbmax ybmin ybmax]);
    
    % myosin plot
    subplot(nrows,ncols,ploti)
    prop = meanI;
    plot(binC*res, meanI, '-')
    hold on
    plot(edgePos*ones([1 2]),[min(meanI) max(meanI)],'--k')
    hold off
    title(setlabels{specialC}, 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel(distlabel);
    ylabel([setlabels{specialC} ' intensity']);
    set(gca,'LineWidth',2);
    axis([binE(1)*res binE(end)*res min(prop) max(prop)])
    ploti = ploti + 1;

    % height plot
    subplot(nrows,ncols,ploti)
    meanApicalZmu = -meanApicalZ*res;
    prop = meanApicalZmu-max(meanApicalZmu);
    propError = stdH;
    plot(binC*res, prop, 'r')
    legendentries = {'apical','apical - basal'};
    hold on
    if hasbasal
        plot(binC*res, meanH - max(meanH));
        errorbar(binC*res, meanH - max(meanH), propError);
    else
        legendentries = legendentries(1);
    end
    plot(edgePos*ones([1 2]),[min(meanH) max(meanH)]-max(meanH),'--k')
    hold off
    title('clone height', 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel(distlabel);
    ylabel('height');
    legend(legendentries, 'Location', 'SouthEast');
    set(gca,'LineWidth',2);
    axis([binE(1)*res binE(end)*res min(prop)-max(propError) max(prop)+max(propError)]);
    ploti = ploti + 1;
    
    % anisotropy plot
    subplot(nrows,ncols,ploti)
    prop = meanAnisotropy;
    propError = stdAnisotropy;
    plot(binC*res, prop)
    hold on
    plot(edgePos*ones([1 2]),[min(prop) max(prop)],'--k')
    errorbar(binC*res, prop, propError)
    hold off
    title('apical anisotropy', 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel(distlabel);
    ylabel('anisotropy');
    set(gca,'LineWidth',2);
    axis([binE(1)*res binE(end)*res min(prop)-max(propError) max(prop)+max(propError)]);
    ploti = ploti + 1;
    
    % density plot
    density3 = 1./meanAreas;
    subplot(nrows,ncols,ploti)
    plot(binC*res, density)
    hold on
    plot(binC*res, density3, 'r')
    legendentries = {'n/Amask', '1/<Acell>', '1/<Aprop>'};
    %plot(bins(1:end-1)*res, density2, 'g') '<1/Acell>', 
    if isfield(wing, 'properAreas')
        density4 = 1./meanAreasProper;
        plot(binC*res, density4, 'm')
    else
        legendentries = legendentries(1:2);
    end
    plot(edgePos*ones([1 2]),[min(density) max(density)],'--k')
    hold off
    title('apical density', 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel(distlabel);
    ylabel('density (\mu m^{-2})');
    legend(legendentries, 'Location', 'NorthEast');
    set(gca,'LineWidth',2);
    axis([binE(1)*res binE(end)*res min(density) max(density3)]);
    ploti = ploti + 1;
    
     % orientation plot
    subplot(nrows,ncols,ploti)
    prop = meanNemNormal;
    propError = stdNemNormal;
    plot(binC*res, prop, 'b')
    hold on
    plot(binC*res, meanNemRadial, 'r');
    plot(edgePos*ones([1 2]),[min(prop) max(prop)],'--k')
    %errorbar(bins(1:end-1)*res, prop, propError)
    hold off
    legend({'radial','normal'});
    title('orientation', 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel(distlabel);
    ylabel('2 cos^2\phi - 1');
    set(gca,'LineWidth',2);
    axis([binE(1)*res binE(end)*res [min(prop) max(prop)]]);
    ploti = ploti + 1;
    
    % weighted orientation plot
    subplot(nrows,ncols,ploti)
    %titlefsize = 24;
    prop = meanNemNormalW;
    propError = stdNemNormalW;
    plot(binC*res, prop, 'b','LineWidth',2)
    hold on
    %plot(binC*res, meanNemRadialW, 'r');
    plot(edgePos*ones([1 2]),[min(prop) max(prop)],'--k','LineWidth',2)
    %errorbar(bins(1:end-1)*res, prop, propError)
    %scatter(edgePos, edgeNemNormW, 100,'xr', 'LineWidth',2);
    %scatter(edgePos + 15, away15NemNormW, 100,'xr', 'LineWidth',2);
    %scatter(edgePos, edgeNemRadialW, 100,'xr', 'LineWidth',2);
    %scatter(edgePos + 15, away15NemRadialW, 100,'xr', 'LineWidth',2);
    hold off
    %legend({'radial','normal'});
    title('anisotropy orientation', 'FontSize', titlefsize, 'FontWeight', 'bold');
    xlabel(distlabel, 'FontSize', titlefsize);
    ylabel('tangential alignment S', 'FontSize', titlefsize);
    set(gca,'LineWidth',2);
    axis([binE(1)*res binE(end)*res [min(prop) max(prop)]]);
    
    colormap jet
    
end