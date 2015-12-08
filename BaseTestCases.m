% Copyright (c) 2015, Mateusz Wójcik (mateuszjanwojcik@gmail.com)
% This is free software and it is distributed under the BSD license - see LICENSE file.
% If you use this code for your research, please cite the paper mentioned here:
%  https://github.com/mjwojcik/Art2MonitoringHybridSystem#reference

classdef BaseTestCases < matlab.unittest.TestCase
    methods (TestClassSetup)
        function disableLogger(testCase)  %#ok<MANU>
            logger = log4m.getLogger();
            logger.setCommandWindowLevel(logger.OFF);
            logger.setLogLevel(logger.OFF);
            logger.setActive(false);
        end
    end
    
    methods(Static)
        function D = generate2DdataSet1()
            rng('default');
            [C1, C2, C3, C4, D1, D2, D3, D4] = PrepareData();
            % PlotData(D1, D2, D3, D4);
            
            P1 = D1(:,1:300);
            P2 = [D1(:,301:end) D2 D3 D4(:,1:end-300)];
            P3 = D4(:,end-299:end);

            R1 = randperm(300);
            R2 = randperm(C1+C2+C3+C4-600);
            R3 = randperm(300);

            D = [P1(:,R1) P2(:,R2) P3(:,R3)]'; 
        end
        
        function D = generate2DdataSet2()
            rng('default');
            [C1, C2, C3, C4, D1, D2, D3, D4] = PrepareData();
            D1=[D1(1,:);1-D1(2,:)];
            D2=[D2(1,:);1-D2(2,:)];
            D3=[D3(1,:);1-D3(2,:)];
            D4=[D4(1,:);1-D4(2,:)];

            % PlotData(D1, D2, D3, D4);

            P1 = [D1 D2];
            P2 = D3(:,1:100);
            P3 = [D3(:,100:end) D4];

            R1 = randperm(C1+C2);
            R2 = randperm(100);
            R3 = randperm(C3-100+C4);

            D = [P1(:,R1) P2(:,R2) P3(:,R3)]'; 
        end

        function D = generate2DdataSet3()
            rng('default');
            [C1, C2, C3, C4, D1, D2, D3, D4] = PrepareData();

            % PlotData(D1, D2, D3, D4);

            P1 = [D1 D2];
            P2 = D3(:,1:100);
            P3 = [D3(:,100:end) D4];

            R1 = randperm(C1+C2);
            R2 = randperm(100);
            R3 = randperm(C3-100+C4);

            D = [P1(:,R1) P2(:,R2) P3(:,R3)]'; 
        end

        function D = generate2DdataSet4_1()
            rng('default');
            C1 = 5000;
            D1 = randomPoints( C1, 0.0005, 1 - 0.60, 0.60);
            D2 = randomPoints( C1, 0.0005, 1 - 0.55, 0.55);
            D3 = randomPoints( C1, 0.0005, 1 - 0.50, 0.50);
            D4 = randomPoints( C1, 0.0005, 1 - 0.45, 0.45);

            % PlotData(D1, D2, D3, D4); 

            P1 = [D1 D2 D3 D4];

            R1 = randperm(4*C1);

            D = P1(:,R1)';
        end
        
        function D = generate2DdataSet4_2()
            rng('default');
            C1 = 2000;
            D1 = randomPoints( C1, 0.10, 1 - 0.9, 0.9);
            D2 = randomPoints( C1, 0.20, 1 - 0.2, 0.2);

            % PlotData(D1, D2); 

            P1 = [D1 D2];

            R1 = randperm(C1+C1);

            D = P1(:,R1)'; 
        end
   end
end

function [C1, C2, C3, C4, D1, D2, D3, D4] = PrepareData()
    C1 = 900;
    C2 = 1000;
    C3 = 900;
    C4 = 800;
    D1 = randomPoints( C1, 0.15, 0.84, 0.84);
    D2 = randomPoints( C2, 0.10, 0.59, 0.59);
    D3 = randomPoints( C3, 0.05, 0.26, 0.26);
    D4 = randomPoints( C4, 0.04, 0.15, 0.15);
end

function PlotData(D1, D2, D3, D4) %#ok<DEFNU>
    figure;
    hold all;
    plot (D1(1,:),D1(2,:),'.');
    plot (D2(1,:),D2(2,:),'.');
    if (nargin == 4)
        plot (D3(1,:),D3(2,:),'.');
        plot (D4(1,:),D4(2,:),'.');
    end
    hold off;
    axis([0 1 0 1]);
end

function [ D ] = randomPoints( n, radius, xc, yc)
    theta = rand(1,n)*(2*pi);
    r = sqrt(rand(1,n))*radius;
    x = xc + r.*cos(theta);
    y = yc + r.*sin(theta);
    D = [x; y];
end
