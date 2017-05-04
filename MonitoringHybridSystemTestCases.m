% Copyright (c) 2015, Mateusz Wójcik (mateuszjanwojcik@gmail.com)
% This is free software and it is distributed under the BSD license - see LICENSE file.
% If you use this code for your research, please cite the paper mentioned here:
%  https://github.com/mjwojcik/Art2MonitoringHybridSystem#reference

classdef MonitoringHybridSystemTestCases < BaseTestCases
    methods (Test)        
        
        function art2StereographicProjectionSystemSimulation(testCase)
            D = BaseTestCases.generate2DdataSet1();
            art2params = Art2.getDefaultParams(size(D,2)+1);
            art2params.theta = 0.01;
            art2params.f2f1WeightRatio = 0.95;            
            art2params.learningLength = 5;
            art2params.vigilance = 0.997;
            art2 = Art2(art2params);

            scallingBounds = MonitoringHybridSystem.calculateBounds(D);    
            scaled_D = MonitoringHybridSystem.scaleMatrix(D, scallingBounds);
            
            nbOfPoints = size(scaled_D, 1);
            eventCounter = 0;
            results = zeros(1, nbOfPoints);
            classMapping = zeros(1,4);
            for i=1:nbOfPoints
                scaled_X = scaled_D(i,:);
                scaled_X = MonitoringHybridSystem.transformUsingStereographicProjection(scaled_X);
                [event, results(i)] = art2.processPoint(scaled_X);
                if (event == true)
                    eventCounter = eventCounter + 1;
                end
                X1 = D(i,1);
                if (X1 < 0.18)
                    classMapping(results(i)) = 1;
                elseif (X1 < 0.4)
                    classMapping(results(i)) = 2;
                elseif (X1 < 0.65)
                    classMapping(results(i)) = 3;
                else
                    classMapping(results(i)) = 4;
                end
                    
            end
            testCase.assertEqual(eventCounter, 4);
            Sets = cell(1, 4);
            Sets{classMapping(1)}.x = D(results==1,1);
            Sets{classMapping(1)}.y = D(results==1,2);
            Sets{classMapping(2)}.x = D(results==2,1);
            Sets{classMapping(2)}.y = D(results==2,2);
            Sets{classMapping(3)}.x = D(results==3,1);
            Sets{classMapping(3)}.y = D(results==3,2);
            Sets{classMapping(4)}.x = D(results==4,1);
            Sets{classMapping(4)}.y = D(results==4,2);
            
            testCase.verifyTrue(max(Sets{1}.x) < 0.21);
            testCase.verifyTrue(max(Sets{1}.y) < 0.21);
            testCase.verifyTrue(min(Sets{1}.x) > 0.1);
            testCase.verifyTrue(min(Sets{1}.y) > 0.1);
            
            testCase.verifyTrue(max(Sets{2}.x) < 0.4);
            testCase.verifyTrue(max(Sets{2}.y) < 0.4);
            testCase.verifyTrue(min(Sets{2}.x) > 0.19);
            testCase.verifyTrue(min(Sets{2}.y) > 0.19);

            testCase.verifyTrue(max(Sets{3}.x) < 0.8);
            testCase.verifyTrue(max(Sets{3}.y) < 0.8);
            testCase.verifyTrue(min(Sets{3}.x) > 0.45);
            testCase.verifyTrue(min(Sets{3}.y) > 0.45);

            testCase.verifyTrue(max(Sets{4}.x) < 1);
            testCase.verifyTrue(max(Sets{4}.y) < 1);
            testCase.verifyTrue(min(Sets{4}.x) > 0.65);
            testCase.verifyTrue(min(Sets{4}.y) > 0.65);

            % plotData(D, results); % debug function
        end
        
        
        function hybridSystemDetermineSpecialPoints(testCase)
            D = BaseTestCases.generate2DdataSet3();
            D = D(1:1100,:);
            D_scaled = MonitoringHybridSystem.scaleMatrix(D, MonitoringHybridSystem.calculateBounds(D));
            specialPoints = MonitoringHybridSystem.determineSpecialPoints(D_scaled, 30);
            %plotData(specialPoints, 5*ones(1,size(specialPoints,2)));
            testCase.assertEqual(size(specialPoints), [2 2]);
        end      
       
        function hybridSystemDetermineSpecialPoints2(testCase)
            D1 = BaseTestCases.generate2DdataSet4_1();
            
            specialPoints = MonitoringHybridSystem.determineSpecialPoints(D1, 30);
            testCase.assertEqual(size(specialPoints), [4 2]);
            
            sortedpoints = sortrows(specialPoints);
            testCase.verifyLessThanOrEqual(sortedpoints(1,1), 1 - 0.59);
            testCase.verifyGreaterThanOrEqual(sortedpoints(1,1), 1 - 0.61);
            testCase.verifyLessThanOrEqual(sortedpoints(1,2), 0.61);
            testCase.verifyGreaterThanOrEqual(sortedpoints(1,2), 0.59);

            testCase.verifyLessThanOrEqual(sortedpoints(2,1), 1 - 0.54);
            testCase.verifyGreaterThanOrEqual(sortedpoints(2,1), 1 - 0.56);
            testCase.verifyLessThanOrEqual(sortedpoints(2,2), 0.56);
            testCase.verifyGreaterThanOrEqual(sortedpoints(2,2), 0.54);

            testCase.verifyLessThanOrEqual(sortedpoints(3,1), 1 - 0.49);
            testCase.verifyGreaterThanOrEqual(sortedpoints(3,1), 1 - 0.51);
            testCase.verifyLessThanOrEqual(sortedpoints(3,2), 0.51);
            testCase.verifyGreaterThanOrEqual(sortedpoints(3,2), 0.49);

            testCase.verifyLessThanOrEqual(sortedpoints(4,1), 1 - 0.44);
            testCase.verifyGreaterThanOrEqual(sortedpoints(4,1), 1 - 0.46);
            testCase.verifyLessThanOrEqual(sortedpoints(4,2), 0.46);
            testCase.verifyGreaterThanOrEqual(sortedpoints(4,2), 0.44);

            
            art2_1 = Art2(Art2.getDefaultParams(2));   
            art2_1.theta = 0; 
            art2_1.setupWeights(specialPoints);
            testCase.verifyEqual(art2_1.nbOfClasses, 4);
            testCase.verifyGreaterThan(art2_1.vigilance, 0.99);
            testCase.verifyLessThanOrEqual(art2_1.vigilance, 1);
            
            % plotArt2LTM(D1, art2_1); % debug function
            
            D2 = BaseTestCases.generate2DdataSet4_2();
            
            specialPoints = MonitoringHybridSystem.determineSpecialPoints(D2, 30);
            testCase.assertEqual(size(specialPoints), [2 2]);
            
            sortedpoints = sortrows(specialPoints);
            testCase.verifyLessThanOrEqual(sortedpoints(1,1), 0.2);
            testCase.verifyGreaterThanOrEqual(sortedpoints(1,1), 0);
            testCase.verifyLessThanOrEqual(sortedpoints(2,1), 1);
            testCase.verifyGreaterThanOrEqual(sortedpoints(2,1), 0.6);
            testCase.verifyLessThanOrEqual(sortedpoints(1,2), 1);
            testCase.verifyGreaterThanOrEqual(sortedpoints(1,2), 0.8);
            testCase.verifyLessThanOrEqual(sortedpoints(2,2), 0.4);
            testCase.verifyGreaterThanOrEqual(sortedpoints(2,2), 0);
                       
            art2_2 = Art2(Art2.getDefaultParams(2));   
            art2_2.theta = 0; 
            art2_2.setupWeights(specialPoints);
            testCase.verifyEqual(art2_2.nbOfClasses, 2);
            testCase.verifyLessThan(art2_2.vigilance, 0.99);
            testCase.verifyLessThan(art2_2.vigilance, art2_1.vigilance);           
            
            % plotArt2LTM(D2, art2_2); % debug function            
        end    
  
        function hybridSystemCreateOCC(testCase)
            D = BaseTestCases.generate2DdataSet3;
            D_scaled = MonitoringHybridSystem.scaleMatrix(D, MonitoringHybridSystem.calculateBounds(D));
                        
            hybridparams = MonitoringHybridSystem.getDefaultParams(size(D_scaled,2));   
            hybridparams.DeltaTcheckStability = 250;
            hybridparams.MaxClusterCountRatio = 1;
            art2params = Art2.getDefaultParams(-1);            
            hsystem = MonitoringHybridSystem(hybridparams, art2params);
            
            testCase.verifyEqual(hsystem.history, zeros(500,2));
            testCase.verifyEqual(hsystem.historyCount, 0);
            
            nbOfPoints = size(D, 1);
            for i=1:nbOfPoints
                X = D_scaled(i,:);
                hsystem.addToHistory(X);
            end
            
            testCase.verifyEqual(size(hsystem.history), [3750,2]);
            testCase.verifyEqual(hsystem.historyCount, nbOfPoints);
            
            hsystem.ART2.nbOfClasses = 5;
            hsystem.addNewOCC();
                                               
            testCase.verifyEqual(size(hsystem.history), [3750,2]);
            testCase.verifyEqual(hsystem.historyCount, 0);
            testCase.verifyEqual(hsystem.occCount, 1);
            
            % plotOCCresult(hsystem.vOCC(1).treshold, hsystem.vOCC(1).em_model, D_scaled);
            
            testCase.verifyEqual(hsystem.isAlreadyKnown([0.1 0.1]), 1);
            testCase.verifyEqual(hsystem.isAlreadyKnown([0.2 0.2]), 1);
            testCase.verifyEqual(hsystem.isAlreadyKnown([0.3 0.3]), 2);
            testCase.verifyEqual(hsystem.isAlreadyKnown([0.45 0.45]), 1);
            testCase.verifyEqual(hsystem.isAlreadyKnown([0.65 0.65]), 1);
        end
        
        function hybridSystemSimulation2withoutStereographicProjection(testCase)
            D = BaseTestCases.generate2DdataSet2;
                        
            hybridparams = MonitoringHybridSystem.getDefaultParams(size(D,2));
            hybridparams.EnableStereographicProjection = false;            
            hybridparams.MaxClusterCountRatio = 1.5;
            
            art2params = Art2.getDefaultParams(-1);
            art2params.theta = 0.001;

            % D_scaled = MonitoringHybridSystem.scaleMatrix(D, MonitoringHybridSystem.calculateBounds(D(1:hybridparams.DeltaTstart,:))); % debug matrix         
            
            hsystem = MonitoringHybridSystem(hybridparams, art2params);
                        
            nbOfPoints = size(D, 1);
            eventLog = [];
            eventLog2 = [];
            vigilance = -1;
            resultAreas = zeros(nbOfPoints,1);
            resultClusters = zeros(nbOfPoints,1);
            for i=1:nbOfPoints
                X = D(i,:);
                [newarea, newcluster, resultAreas(i), resultClusters(i)] = hsystem.processPoint(X);
                if (newarea == true)
                    eventLog = [eventLog i]; %#ok<AGROW>
                end
                if (newcluster == true)
                    eventLog2 = [eventLog2 i]; %#ok<AGROW>
                    % plotArt2LTM(D_scaled(1:i,:), hsystem.ART2); % debug function
                end
                if (i == hsystem.DeltaTstart)
                    vigilance = hsystem.ART2.vigilance;
                    testCase.verifyEqual(hsystem.ART2.nbOfClasses, 2);
                end
            end
            testCase.assertEqual(hsystem.occCount, 2);
            testCase.assertEqual(length(eventLog), 2);
            testCase.verifyLessThan(eventLog(1), 1900);
            testCase.verifyGreaterThan(eventLog(2), 1900);
            
            testCase.verifyEqual(resultAreas(1:500), -1 * ones(500,1));
            testCase.verifyEqual(resultAreas(501:1900), ones(1400,1));
            testCase.verifyEqual(resultAreas(1901:3600), 2 * ones(1700,1));

            testCase.verifyEqual(resultClusters(eventLog(1)+1:1900), -1 * ones(1900-eventLog(1),1));
            testCase.verifyEqual(resultClusters(eventLog(2)+1:3600), -1 * ones(3600-eventLog(2),1));

            eventLog2 = eventLog2(eventLog2 > 1900);
            testCase.assertGreaterThanOrEqual(length(eventLog2), 2);
            testCase.verifyEqual(eventLog2(1), 1901);
            
            testCase.verifyGreaterThan(hsystem.ART2.vigilance, vigilance);
            testCase.verifyLessThan(hsystem.ART2.vigilance, 1);
        end        

        function hybridSystemSimulation2withStereographicProjection(testCase)
            D = BaseTestCases.generate2DdataSet2;
                        
            hybridparams = MonitoringHybridSystem.getDefaultParams(size(D,2));
            hybridparams.EnableStereographicProjection = true;            
            hybridparams.MaxClusterCountRatio = 1.5;
            
            art2params = Art2.getDefaultParams(-1);
            art2params.theta = 0.001;

            % D_scaled = MonitoringHybridSystem.scaleMatrix(D, MonitoringHybridSystem.calculateBounds(D(1:hybridparams.DeltaTstart,:))); % debug matrix         
            
            hsystem = MonitoringHybridSystem(hybridparams, art2params);
                        
            nbOfPoints = size(D, 1);
            eventLog = [];
            eventLog2 = [];
            vigilance = -1;
            resultAreas = zeros(nbOfPoints,1);
            resultClusters = zeros(nbOfPoints,1);
            for i=1:nbOfPoints
                X = D(i,:);
                [newarea, newcluster, resultAreas(i), resultClusters(i)] = hsystem.processPoint(X);
                if (newarea == true)
                    eventLog = [eventLog i]; %#ok<AGROW>
                end
                if (newcluster == true)
                    eventLog2 = [eventLog2 i]; %#ok<AGROW>
                    % plotArt2LTM(D_scaled(1:i,:), hsystem.ART2); % debug function
                end
                if (i == hsystem.DeltaTstart)
                    vigilance = hsystem.ART2.vigilance;
                    testCase.verifyEqual(hsystem.ART2.nbOfClasses, 2);
                end
            end
            testCase.assertEqual(hsystem.occCount, 2);
            testCase.assertEqual(length(eventLog), 2);
            testCase.verifyLessThan(eventLog(1), 1900);
            testCase.verifyGreaterThan(eventLog(2), 1900);
            
            testCase.verifyEqual(resultAreas(1:500), -1 * ones(500,1));
            testCase.verifyEqual(resultAreas(501:1900), ones(1400,1));
            testCase.verifyEqual(resultAreas(1901:3600), 2 * ones(1700,1));

            testCase.verifyEqual(resultClusters(eventLog(1)+1:1900), -1 * ones(1900-eventLog(1),1));
            testCase.verifyEqual(resultClusters(eventLog(2)+1:3600), -1 * ones(3600-eventLog(2),1));

            eventLog2 = eventLog2(eventLog2 > 1900);
            testCase.assertGreaterThanOrEqual(length(eventLog2), 2);
            testCase.verifyEqual(eventLog2(1), 1901);
            
            testCase.verifyGreaterThan(hsystem.ART2.vigilance, vigilance);
            testCase.verifyLessThan(hsystem.ART2.vigilance, 1);
        end
        
        function hybridSystemSimulation3(testCase)
            D = BaseTestCases.generate2DdataSet3;
                        
            hybridparams = MonitoringHybridSystem.getDefaultParams(size(D,2));
            hybridparams.DeltaTstart = 1100;
            hybridparams.EnableStereographicProjection = true;            
            hybridparams.MaxClusterCountRatio = 1.5;
            
            art2params = Art2.getDefaultParams(-1);
            art2params.theta = 0.001;

            % D_scaled = MonitoringHybridSystem.scaleMatrix(D, MonitoringHybridSystem.calculateBounds(D(1:hybridparams.DeltaTstart,:))); % debug matrix         
            
            hsystem = MonitoringHybridSystem(hybridparams, art2params);
                        
            nbOfPoints = size(D, 1);
            eventLog = [];
            eventLog2 = [];
            vigilance = -1;
            resultAreas = zeros(nbOfPoints,1);
            resultClusters = zeros(nbOfPoints,1);
            for i=1:nbOfPoints
                X = D(i,:);
                [newarea, newcluster, resultAreas(i), resultClusters(i)] = hsystem.processPoint(X);
                if (newarea == true)
                    eventLog = [eventLog i]; %#ok<AGROW>
                end
                if (newcluster == true)
                    eventLog2 = [eventLog2 i]; %#ok<AGROW>
                    % plotArt2LTM(D_scaled(1:i,:), hsystem.ART2); % debug function
                end
                if (i == hsystem.DeltaTstart)
                    vigilance = hsystem.ART2.vigilance;
                    testCase.verifyEqual(hsystem.ART2.nbOfClasses, 2);
                end
            end
            testCase.assertEqual(hsystem.occCount, 2);
            testCase.assertEqual(length(eventLog), 2);
            testCase.verifyLessThan(eventLog(1), 1900);
            testCase.verifyGreaterThan(eventLog(2), 1900);
            
            testCase.verifyEqual(resultAreas(1:1100), -1 * ones(1100,1));
            testCase.verifyEqual(resultAreas(1101:1900), ones(800,1));
            testCase.verifyEqual(resultAreas(1901:3600), 2 * ones(1700,1));

            testCase.verifyEqual(resultClusters(eventLog(1)+1:1900), -1 * ones(1900-eventLog(1),1));
            testCase.verifyEqual(resultClusters(eventLog(2)+1:3600), -1 * ones(3600-eventLog(2),1));

            eventLog2 = eventLog2(eventLog2 > 1900);
            testCase.assertGreaterThanOrEqual(length(eventLog2), 2);
            testCase.verifyEqual(eventLog2(1), 1901);
            
            testCase.verifyGreaterThan(hsystem.ART2.vigilance, vigilance);
            testCase.verifyLessThan(hsystem.ART2.vigilance, 1);
        end
        
        function hybridSystemSimulation4withoutStereographicProjection(testCase)
            D = [BaseTestCases.generate2DdataSet4_1() ; BaseTestCases.generate2DdataSet4_2()];
                        
            hybridparams = MonitoringHybridSystem.getDefaultParams(size(D,2));
            hybridparams.EnableStereographicProjection = false;            
            hybridparams.DeltaTstart = 14000;
            hybridparams.DeltaTcheckMinMax = 100;
            hybridparams.DeltaTcheckStability = 1000;
            hybridparams.MaxClusterCountRatio = 1;
            hybridparams.ScalingBoundaryPercentage = 1;
            
            art2params = Art2.getDefaultParams(-1);
            art2params.theta = 0.000001;

            % D_scaled = MonitoringHybridSystem.scaleMatrix(D, MonitoringHybridSystem.calculateBounds(D(1:hybridparams.DeltaTstart,:))); % debug matrix         
            
            hsystem = MonitoringHybridSystem(hybridparams, art2params);

            nbOfPoints = size(D, 1);
            vigilance = -1;

            for i=1:nbOfPoints
                X = D(i,:);
                hsystem.processPoint(X);
                if (i == hsystem.DeltaTstart)
                    vigilance = hsystem.ART2.vigilance;
                    testCase.verifyEqual(hsystem.ART2.nbOfClasses, 4);
                    % plotArt2LTM(D_scaled(1:i,:), hsystem.ART2); % debug function
                end
            end
            testCase.verifyLessThan(hsystem.ART2.vigilance, vigilance);
            testCase.verifyGreaterThanOrEqual(hsystem.occCount, 2);
        end  
      
        function hybridSystemSimulation4withStereographicProjection(testCase)
            D = [BaseTestCases.generate2DdataSet4_1() ; BaseTestCases.generate2DdataSet4_2()];
                        
            hybridparams = MonitoringHybridSystem.getDefaultParams(size(D,2));
            hybridparams.EnableStereographicProjection = true;            
            hybridparams.DeltaTstart = 14000;
            hybridparams.DeltaTcheckMinMax = 100;
            hybridparams.DeltaTcheckStability = 1000;
            hybridparams.MaxClusterCountRatio = 1;
            hybridparams.ScalingBoundaryPercentage = 1;
            
            art2params = Art2.getDefaultParams(-1);
            art2params.theta = 0.000001;

            % D_scaled = MonitoringHybridSystem.scaleMatrix(D, MonitoringHybridSystem.calculateBounds(D(1:hybridparams.DeltaTstart,:))); % debug matrix         
            
            hsystem = MonitoringHybridSystem(hybridparams, art2params);

            nbOfPoints = size(D, 1);
            vigilance = -1;

            for i=1:nbOfPoints
                X = D(i,:);
                hsystem.processPoint(X);
                if (i == hsystem.DeltaTstart)
                    vigilance = hsystem.ART2.vigilance;
                    testCase.verifyEqual(hsystem.ART2.nbOfClasses, 4);
                    % plotArt2LTM(D_scaled(1:i,:), hsystem.ART2); % debug function
                end
            end
            testCase.verifyLessThan(hsystem.ART2.vigilance, vigilance);
            testCase.verifyGreaterThanOrEqual(hsystem.occCount, 2);
        end 
    end   
end 

%% helpers

function plotArt2LTM(DATA, art2)  %#ok<DEFNU>
    figure('Position', [0 0 600 500], 'Name', 'ART-2 LTM memory visualization');
    hold all;

    V = applyNorm(DATA);
    plot(V(:,1),V(:,2), '.k');

    for k=1:art2.nbOfClasses
        vf1f2 = art2.f1f2(:,k);
        vf1f2 = vf1f2 ./ norm(vf1f2);
        quiver(0, 0, vf1f2(1), vf1f2(2), '-b');                        
    end
    for k=1:art2.nbOfClasses
        vf2f1 = art2.f2f1(k,:);
        vf2f1 = vf2f1 ./ norm(vf2f1);
        quiver(0, 0, vf2f1(1), vf2f1(2), '--k');
    end

    axis([0 1 0 1]);
    hold off;
end     

function [ N ] = applyNorm( W )
    N = zeros(size(W));
    for i=1:size(W,1)
        N(i,:)=W(i,:)./norm(W(i,:));
    end
end

function plotData(DATA, results) %#ok<DEFNU>
    hold all;
    for k=1:size(DATA,1);  
        P = DATA(k,:);
        plot(P(1), P(2), '.', 'color', getColor(results(k)));                    
    end
    hold off;
end   

function m = getColor(idx)
    K = [lines;hsv(128)];
    m = K(mod(idx,size(K,1))+1,:);
end

function plotOCCresult(tr, em_model, data) %#ok<DEFNU>
    [X, Y, Z] = calcdatagrid(em_model, 100);
    
    FigHandle = figure;
    plot(data(:,1),data(:,2),'.','MarkerSize', 15, 'MarkerEdgeColor',[.5 .5 .5]);
    hold on;
    contour(X,Y,Z, [tr tr], 'k', 'LineWidth',3);
    set(FigHandle, 'Position', [100, 100, 800, 800]);
end

function [X, Y, Z] = calcdatagrid( em_model, bitmap_quality)
    vmin = -0.1;
    vmax =  0.8;
    
    [X, Y] = meshgrid(linspace(vmin, vmax, bitmap_quality), linspace(vmin, vmax, bitmap_quality));
    Z = MonitoringHybridSystem.em_f([X(:) Y(:)], em_model);
    Z = reshape(Z, [bitmap_quality,bitmap_quality]);
end
