% Copyright (c) 2015, Mateusz Wójcik (mateuszjanwojcik@gmail.com)
% This is free software and it is distributed under the BSD license - see LICENSE file.
% If you use this code for your research, please cite the paper mentioned here:
%  https://github.com/mjwojcik/Art2MonitoringHybridSystem#reference

function runSimulation( DATA, dims, art2params, hybridparams)
% Assuming the variable DATA expresses parameteres of machine condition,
% this procedure shows how MonitoringHybridSystem can monitor machine operational states.
% To use this procedure, provide DATA matrix as a vector of multidimensional data points
% and specify three dimensions for which plotting functions will be called (e.g. dims = [1 2 3]).
% The parameters art2params and hybridparams are optional and can be specified also inside the procedure.

    tmpfolder = [pwd() filesep 'logs' filesep 'tmp-' datestr(now, 'yyyy-mm-dd_HH:MM:SS')];
    mkdir(tmpfolder);
    save([tmpfolder filesep 'before_simulation.mat']);
  
    logFile = [tmpfolder filesep 'simulation.log'];
    logger = log4m.getLogger(logFile);
    logger.setCommandWindowLevel(logger.INFO);
    logger.setLogLevel(logger.DEBUG);
            
    if (nargin < 4)
        hybridparams = MonitoringHybridSystem.getDefaultParams(size(DATA,2));
        % hybridparams.DeltaTstart = 1000?;
        % hybridparams.DeltaTcheckMinMax = 100;
        % hybridparams.DeltaTcheckStability = 500;

        art2params = Art2.getDefaultParams(-1);
        art2params.learningRate = 0.01;
        art2params.learningLength = 10;
        art2params.theta = 0.001; 
    end    
    
    hsystem = MonitoringHybridSystem(hybridparams, art2params);
        
    nbOfPoints = size(DATA,1);
    resultAreas = zeros(nbOfPoints,1);
    resultClusters = zeros(nbOfPoints,1);
    areaCount = 0;
    for i=1:nbOfPoints
        X = DATA(i,:);
        [newarea, newcluster, resultAreas(i), resultClusters(i)] = hsystem.processPoint(X);
        if (newarea)
            areaCount = areaCount + 1;
            if (areaCount > 1)
                plotAll(DATA, resultAreas, resultClusters, areaCount, i, tmpfolder, dims);
            end
        elseif (newcluster)
            plotAll(DATA, resultAreas, resultClusters, areaCount+1, i, tmpfolder, dims);
        end
    end
    plotAll(DATA, resultAreas, resultClusters, areaCount+1, nbOfPoints, tmpfolder, dims);
    
    save([tmpfolder filesep 'after_simulation.mat']);
    
    logger.close();
end

function plotAll(DATA, resultAreas, resultClusters, lastAreaCount, lastPointNumber, path, dims)
    h = figure('Position', [0 0 800 400], 'Name', sprintf('The classification after %d data points', lastPointNumber));
    hold all;
        
    for areaidx=1:lastAreaCount-1
        filter = resultAreas == areaidx;
        plot3(DATA(filter, dims(1)), DATA(filter, dims(2)), DATA(filter, dims(3)), getSymbol(areaidx), 'color', getColor(areaidx));
    end         
    
    mainfilter = (resultAreas == lastAreaCount);
    maxClusterIdx = max(resultClusters(mainfilter));
    for clusteridx=1:maxClusterIdx
        filter = mainfilter & (resultClusters == clusteridx);
        plot3(DATA(filter, dims(1)), DATA(filter, dims(2)), DATA(filter, dims(3)), '.', 'color', getColor(clusteridx+20));
    end
    
    hold off;
    view(45, 45);
    set(h, 'PaperPositionMode', 'auto');
    destination_filename = sprintf('%s/%d', path, lastPointNumber);
    print(h, '-r0', [destination_filename '.png'], '-dpng');
    saveas(h,[destination_filename '.fig']);     
end

function m = getSymbol(idx)
    markers = {'+','o','*','x','s','d','^','v','>','<','p','h'};
    m = markers{mod(idx,numel(markers))+1};
end

function m = getColor(idx)
    K = [lines;hsv(128)];
    m = K(mod(idx,numel(K))+1,:);
end
