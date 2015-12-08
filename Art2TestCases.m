% Copyright (c) 2015, Mateusz Wójcik (mateuszjanwojcik@gmail.com)
% This is free software and it is distributed under the BSD license - see LICENSE file.
% If you use this code for your research, please cite the paper mentioned here:
%  https://github.com/mjwojcik/Art2MonitoringHybridSystem#reference

classdef Art2TestCases < BaseTestCases
    methods (Test)        
        
        function fastLearningNoiseSuppressionStepByStepTest(testCase)
            art2params = Art2.getDefaultParams(2);
            art2params.e = 0;
            art2params.theta = 0.7;
            art2params.f2f1WeightRatio = 1.0;
            
            art2 = Art2(art2params);            
            
            testCase.verifyEqual(art2.f2f1(1,:), [0 0]);
            testCase.verifyEqual(art2.f2f1(10,:), [0 0]);
            testCase.verifyEqual(art2.f2f1(100,:), [0 0]);
            testCase.verifyEqual(art2.f1f2(:,1), [7.07 ; 7.07], 'AbsTol', 0.01);
            testCase.verifyEqual(art2.f1f2(:,10), [7.07 ; 7.07], 'AbsTol', 0.01);
            testCase.verifyEqual(art2.f1f2(:,100), [7.07 ; 7.07], 'AbsTol', 0.01);
            
            S = [0.8 0.6];
            
            art2.f1_init(S);            
            testCase.verifyEqual(art2.f1_U, [0 0]);
            testCase.verifyEqual(art2.f1_W, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_P, [0 0]);
            testCase.verifyEqual(art2.f1_X, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_Q, [0 0]);
            testCase.verifyEqual(art2.f1_V, [0.8 0]);   
           
            art2.f1_update(S, -1);
            testCase.verifyEqual(art2.f1_U, [1 0]);
            testCase.verifyEqual(art2.f1_W, [10.8 0.6]);
            testCase.verifyEqual(art2.f1_P, [1 0]);
            testCase.verifyEqual(art2.f1_X, [0.998 0.055], 'AbsTol', 0.001);
            testCase.verifyEqual(art2.f1_Q, [1 0]);
            testCase.verifyEqual(art2.f1_V, [10.998 0], 'AbsTol', 0.001);   
            
            Y = art2.simulateF1F2();
            testCase.verifyEqual(art2.nbOfClasses, 0);
            testCase.verifyEqual(size(Y), [1 1]);
            
            testCase.verifyEqual(art2.calculateVigilance(1), 1);
            testCase.verifyTrue(art2.isActiveNeuronOK(1));
            
            art2.learningLength = 0;
            art2.learning(S, 1);
            testCase.verifyEqual(art2.f1_W, [10.8 0.6]);
            testCase.verifyEqual(art2.f1_X, [0.998 0.055], 'AbsTol', 0.001);
            testCase.verifyEqual(art2.f1_Q, [1 0]);
            testCase.verifyEqual(art2.f1_V, [10.998 0], 'AbsTol', 0.001);   

            art2.adjustWeights(1);
            testCase.verifyEqual(art2.f2f1(1,:), [.6*.9*1+(1-.6*.9*.1)*0      .6*.9*0+(1-.6*.9*.1)*0]);
            testCase.verifyEqual(art2.f1f2(:,1), [.6*.9*1+(1-.6*.9*.1)*7.07 ; .6*.9*0+(1-.6*.9*.1)*7.07], 'AbsTol', 0.01);
            art2.f1_update(S, 1);
            testCase.verifyEqual(art2.f1_U, [1 0]);
            testCase.verifyEqual(art2.f1_W, [10.8 0.6]);
            testCase.verifyEqual(art2.f1_P, [1+0.9*0.54 0]);
            testCase.verifyEqual(art2.f1_X, [0.998 0.055], 'AbsTol', 0.001);
            testCase.verifyEqual(art2.f1_Q, [1 0]);
            testCase.verifyEqual(art2.f1_V, [0.998+10*1 0], 'AbsTol', 0.001);   

            art2.adjustWeights(1);
            testCase.verifyEqual(art2.f2f1(1,:), [.6*.9*1+(1-.6*.9*.1)*0.54   .6*.9*0+(1-.6*.9*.1)*0]);
            testCase.verifyEqual(art2.f1f2(:,1), [.6*.9*1+(1-.6*.9*.1)*7.229 ; .6*.9*0+(1-.6*.9*.1)*6.689], 'AbsTol', 0.01);
            art2.f1_update(S, 1);
            testCase.verifyEqual(art2.f1_U, [1 0]);
            testCase.verifyEqual(art2.f1_W, [10.8 0.6]);
            testCase.verifyEqual(art2.f1_P, [1+0.9*1.05 0], 'AbsTol', 0.001);  
            testCase.verifyEqual(art2.f1_X, [0.998 0.055], 'AbsTol', 0.001);
            testCase.verifyEqual(art2.f1_Q, [1 0]);
            testCase.verifyEqual(art2.f1_V, [0.998+10*1 0], 'AbsTol', 0.001);           

            art2.learningLength = 1000;
            art2.learning(S, 1);
            testCase.verifyEqual(art2.f2f1(1,:), 1./(1-art2.d)*[1 0], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2.f1f2(:,1), 1./(1-art2.d)*[1;0], 'AbsTol', 0.00001);          
            testCase.verifyEqual(art2.f1_U, [1 0]);
            testCase.verifyEqual(art2.f1_Q, [1 0]);
        end

        function fastLearningNoiseSuppressionDetermineClassTest(testCase)
            art2params = Art2.getDefaultParams(2);
            art2params.e = 0;
            art2params.theta = 0.7;
            art2params.f2f1WeightRatio = 1.0;
            
            art2 = Art2(art2params);                        
            
            S = [0.8 0.6];
            
            testCase.verifyEqual(art2.determineClass(S, false), 1);
            testCase.verifyEqual(art2.f1_U, [1 0]);
            testCase.verifyEqual(art2.f1_W, [10.8 0.6]);
            testCase.verifyEqual(art2.f1_P, [1 0]);
            testCase.verifyEqual(art2.f1_X, [0.998 0.055], 'AbsTol', 0.001);
            testCase.verifyEqual(art2.f1_Q, [1 0]);
            testCase.verifyEqual(art2.f1_V, [10.998 0], 'AbsTol', 0.001);   
        end

        function fastLearningNoiseSuppressionTest(testCase)
            art2params = Art2.getDefaultParams(2);
            art2params.e = 0;
            art2params.theta = 0.7;
            art2params.f2f1WeightRatio = 1.0;
            art2params.learningLength = 1000;
            
            art2 = Art2(art2params);                        
            
            S = [0.8 0.6];
            
            [newcluster, neuronIdx] = art2.process(S, false);
            testCase.verifyEqual(art2.nbOfClasses, 1);
            testCase.verifyEqual(neuronIdx, 1);
            testCase.verifyTrue(newcluster);
            testCase.verifyEqual(art2.f2f1(1,:), 1./(1-art2.d)*[1 0], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2.f1f2(:,1), 1./(1-art2.d)*[1;0], 'AbsTol', 0.00001);          
            testCase.verifyEqual(art2.f1_U, [1 0]);
            testCase.verifyEqual(art2.f1_Q, [1 0]);
        end
        
        function fastLearningNoNoiseSuppressionPart1StepByStepTest(testCase)
            art2params = Art2.getDefaultParams(2);
            art2params.e = 0;
            art2params.theta = 0.1;
            art2params.f2f1WeightRatio = 1.0;
            
            art2 = Art2(art2params);            
            
            testCase.verifyEqual(art2.f2f1(1,:), [0 0]);
            testCase.verifyEqual(art2.f2f1(10,:), [0 0]);
            testCase.verifyEqual(art2.f2f1(100,:), [0 0]);
            testCase.verifyEqual(art2.f1f2(:,1), [7.07 ; 7.07], 'AbsTol', 0.01);
            testCase.verifyEqual(art2.f1f2(:,10), [7.07 ; 7.07], 'AbsTol', 0.01);
            testCase.verifyEqual(art2.f1f2(:,100), [7.07 ; 7.07], 'AbsTol', 0.01);
            
            S1 = [0.8 0.6];
            
            art2.f1_init(S1);            
            testCase.verifyEqual(art2.f1_U, [0 0]);
            testCase.verifyEqual(art2.f1_W, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_P, [0 0]);
            testCase.verifyEqual(art2.f1_X, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_Q, [0 0]);
            testCase.verifyEqual(art2.f1_V, [0.8 0.6]);   
           
            art2.f1_update(S1, -1);
            testCase.verifyEqual(art2.f1_U, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_W, [8.8 6.6]);
            testCase.verifyEqual(art2.f1_P, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_X, [0.8 0.6], 'AbsTol', 0.001);
            testCase.verifyEqual(art2.f1_Q, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_V, [0.8 0.6] + 10.*[0.8 0.6], 'AbsTol', 0.001);
            
            Y = art2.simulateF1F2();
            testCase.verifyEqual(art2.nbOfClasses, 0);
            testCase.verifyEqual(size(Y), [1 1]);            

            art2.learningLength = 1;
            art2.learning(S1, 1);
            testCase.verifyEqual(art2.f2f1(1,:), 0.6*0.9.*[0.8 0.6] + (1-.6*.9*.1).*[0 0], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2.f1_U, [0.8 0.6]);
            testCase.verifyEqual(art2.f1_W, [0.8 0.6] + 10.*[0.8 0.6]);
            testCase.verifyEqual(art2.f1_P, [0.8 0.6] + 0.9.*[0.432 0.324]);
            testCase.verifyEqual(art2.f1_X, [0.8 0.6], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2.f1_Q, [0.8 0.6], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2.f1_V, [0.8 0.6] + 10.*[0.8 0.6], 'AbsTol', 0.001);   

            art2.adjustWeights(1);
            testCase.verifyEqual(art2.f2f1(1,:), 0.6*0.9.*[0.8 0.6] + (1-.6*.9*.1)*0.54.*[0.8 0.6], 'AbsTol', 0.00001);

            art2.learningLength = 1000;
            art2.learning(S1, 1);
            testCase.verifyEqual(art2.f2f1(1,:), 1./(1-art2.d)*[0.8 0.6], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2.f1f2(:,1), 1./(1-art2.d)*[0.8;0.6], 'AbsTol', 0.00001); 

            testCase.verifyEqual(art2.nbOfClasses, 0);
            art2.nbOfClasses = 1;
            
            S2 = [0.6 0.8];
            
            art2.f1_init(S2);            
            testCase.verifyEqual(art2.f1_U, [0 0]);
            testCase.verifyEqual(art2.f1_W, [0.6 0.8]);
            testCase.verifyEqual(art2.f1_P, [0 0]);
            testCase.verifyEqual(art2.f1_X, [0.6 0.8]);
            testCase.verifyEqual(art2.f1_Q, [0 0]);
            testCase.verifyEqual(art2.f1_V, [0.6 0.8]);   
           
            art2.f1_update(S2, -1)
            testCase.verifyEqual(art2.f1_U, [0.6 0.8]);
            testCase.verifyEqual(art2.f1_W, [6.6 8.8]);
            testCase.verifyEqual(art2.f1_P, [0.6 0.8]);
            testCase.verifyEqual(art2.f1_X, [0.6 0.8], 'AbsTol', 0.001);
            testCase.verifyEqual(art2.f1_Q, [0.6 0.8]);
            testCase.verifyEqual(art2.f1_V, [0.6 0.8] + 10.*[0.6 0.8], 'AbsTol', 0.001);
            
            Y = art2.simulateF1F2();
            testCase.verifyEqual(size(Y), [2 1]);    
            testCase.verifyEqual(Y, [9.6 ; 9.9], 'AbsTol', 0.01);
        end
                
        function fastLearningNoNoiseSuppressionPart2Test(testCase)
            art2params = Art2.getDefaultParams(2);
            art2params.e = 0;
            art2params.theta = 0.1;
            art2params.f2f1WeightRatio = 0.71;            
            art2params.learningLength = 1000;
            
            art2_a = Art2(art2params);   
            
            testCase.verifyEqual(art2_a.f2f1(1,:), [0 0]);
            testCase.verifyEqual(art2_a.f2f1(10,:), [0 0]);
            testCase.verifyEqual(art2_a.f2f1(100,:), [0 0]);
            testCase.verifyEqual(art2_a.f1f2(:,1), [5.02 ; 5.02], 'AbsTol', 0.01);
            testCase.verifyEqual(art2_a.f1f2(:,10), [5.02 ; 5.02], 'AbsTol', 0.01);
            testCase.verifyEqual(art2_a.f1f2(:,100), [5.02 ; 5.02], 'AbsTol', 0.01);
            
            testCase.verifyEqual(art2_a.nbOfClasses, 0);
            [newcluster, neuronIdx] = art2_a.process([0.8 0.6], false);
            testCase.verifyEqual(art2_a.nbOfClasses, 1);
            testCase.verifyEqual(neuronIdx, 1);
            testCase.verifyTrue(newcluster);
            testCase.verifyEqual(art2_a.f2f1(1,:), 1./(1-art2_a.d)*[0.8 0.6], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2_a.f1f2(:,1), 1./(1-art2_a.d)*[0.8;0.6], 'AbsTol', 0.00001); 

            art2_a.learningLength = 1;
            
            [newcluster, neuronIdx] = art2_a.process([0.6 0.8], false);
            testCase.verifyEqual(art2_a.nbOfClasses, 1);
            testCase.verifyEqual(neuronIdx, 1);
            testCase.verifyEqual(newcluster, false);
            testCase.verifyEqual(art2_a.f2f1(1,:), [7.9   6.1], 'AbsTol', 0.1);
            testCase.verifyEqual(art2_a.f1f2(:,1), [7.9 ; 6.1], 'AbsTol', 0.1);         
            
            art2params.vigilance = 0.995;
            art2_b = Art2(art2params);   
                       
            [newcluster, neuronIdx] = art2_b.process([0.8 0.6], false);
            testCase.verifyEqual(art2_b.nbOfClasses, 1);
            testCase.verifyEqual(neuronIdx, 1);
            testCase.verifyTrue(newcluster);
            testCase.verifyEqual(art2_b.f2f1(1,:), 1./(1-art2_b.d)*[0.8 0.6], 'AbsTol', 0.00001);
            testCase.verifyEqual(art2_b.f1f2(:,1), 1./(1-art2_b.d)*[0.8;0.6], 'AbsTol', 0.00001); 
            
            [newcluster, neuronIdx] = art2_b.process([0.6 0.8], false);
            testCase.verifyEqual(art2_b.nbOfClasses, 2);
            testCase.verifyEqual(neuronIdx, 2);
            testCase.verifyTrue(newcluster);
            testCase.verifyEqual(art2_b.f2f1(2,:), 1./(1-art2_b.d)*[0.6 0.8], 'AbsTol', 0.00001); 
            testCase.verifyEqual(art2_b.f1f2(:,2), 1./(1-art2_b.d)*[0.6;0.8], 'AbsTol', 0.00001); 
        end
        
        function art2SystemSimulation1(testCase)
            D = BaseTestCases.generate2DdataSet1();
            art2params = Art2.getDefaultParams(size(D,2));
            art2params.theta = 0.01;
            art2params.f2f1WeightRatio = 0.95;            
            art2params.learningLength = 15;
            art2params.vigilance = 0.9985;
            art2 = Art2(art2params);

            
            nbOfPoints = size(D, 1);
            eventCounter = 0;
            results = zeros(1, nbOfPoints);
            for i=1:nbOfPoints
                X = D(i,:);
                [event, results(i)] = art2.processPoint(X);
                if (event == true)
                    eventCounter = eventCounter + 1;
                end
            end
            
            %testCase.assertEqual(eventCounter, 3);
            Set1x = D(results==1,1);
            Set1y = D(results==1,2);
            Set2x = D(results==2,1);
            Set2y = D(results==2,2);
            Set3x = D(results==3,1);
            Set3y = D(results==3,2);
            
            testCase.verifyTrue(min(Set1x) < 0.2);
            testCase.verifyTrue(min(Set1y) < 0.2);
            testCase.verifyTrue(max(Set1x) > 0.8);
            testCase.verifyTrue(max(Set1y) > 0.8);
            
            testCase.verifyTrue(min(Set2x) < 0.2);
            testCase.verifyTrue(min(Set2y) < 0.2);
            testCase.verifyTrue(max(Set2x) > 0.8);
            testCase.verifyTrue(max(Set2y) > 0.8);

            testCase.verifyTrue(min(Set3x) < 0.2);
            testCase.verifyTrue(min(Set3y) < 0.2);
            testCase.verifyTrue(max(Set3x) > 0.8);
            testCase.verifyTrue(max(Set3y) > 0.8);
            
            % plotArt2LTM(D, art2); % debug function
            % plotData(D, results); % debug function
        end

        function art2SystemSimulation2(testCase)
            D = BaseTestCases.generate2DdataSet2();
            art2params = Art2.getDefaultParams(size(D,2));
            art2params.theta = 0.01;
            art2params.f2f1WeightRatio = 0.95;            
            art2params.learningLength = 5;
            art2params.vigilance = 0.975;
            art2 = Art2(art2params);

            nbOfPoints = size(D, 1);
            eventLog = [];
            for i=1:nbOfPoints
                X = D(i,:);
                [event, ~] = art2.processPoint(X);
                if (event == true)
                    eventLog = [eventLog i]; %#ok<AGROW>
                end
            end
            testCase.assertEqual(length(eventLog), 3);
            testCase.verifyEqual(eventLog(1), 1);
            testCase.verifyEqual(eventLog(2), 26);
            testCase.verifyEqual(eventLog(3), 1901);
            
            % plotArt2LTM(D, art2); % debug function
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
