% Copyright (c) 2015, Mateusz Wójcik (mateuszjanwojcik@gmail.com)
% This is free software and it is distributed under the BSD license - see LICENSE file.
% If you use this code for your research, please cite the paper mentioned here:
%  https://github.com/mjwojcik/Art2MonitoringHybridSystem#reference

classdef Art2 < handle
 
    properties (Access=public) 
        % LTM
        f1f2
        f2f1

        % STM
        f1_U
        f1_W
        f1_P
        f1_X
        f1_Q
        f1_V
        
        % main parameters
        a
        b
        c
        d
        e
        theta
        vigilance
        
        % others
        capacity
        dataDim
        iterationCounter
        learningLength     
        learningRate
        nbOfClasses
        f2f1WeightRatio        
    end
    
    methods(Static)
        function [ art2params ] = getDefaultParams(dataDim)
        art2params.capacity = 100;
        art2params.a = 10;
        art2params.b = 10;
        art2params.c = 0.1;
        art2params.d = 0.9;
        art2params.e = eps;
        art2params.learningLength = 5;
        art2params.learningRate = 0.6;
        art2params.dataDim = dataDim;
        art2params.theta = 1/sqrt(dataDim);
        art2params.vigilance = 0.9;
        art2params.f2f1WeightRatio = 0.95;
        end
    end
    
    methods (Access=private)     
        function f1_update_U(net)
            net.f1_U = net.f1_V./(net.e+norm(net.f1_V));
        end
        
        function f1_update_W(net, S)
            net.f1_W = S + net.f1_U .* net.a;
        end 
        
        function f1_update_P(net, activeNeuronIdx)
            if (activeNeuronIdx == -1)
                net.f1_P = net.f1_U;
            else                
                net.f1_P = net.f1_U + net.f2f1(activeNeuronIdx,:).*net.d;
            end
        end
        
        function f1_update_X(net)
            net.f1_X = net.f1_W ./ (net.e + norm(net.f1_W));
        end
        
        function f1_update_Q(net)
            net.f1_Q = net.f1_P ./ (net.e + norm(net.f1_P));
        end
        
        function f1_update_V(net)
            net.f1_V = net.noiseSuppressiom(net.f1_X) + net.noiseSuppressiom(net.f1_Q) .* net.b;                
        end        
                
        function [R] = noiseSuppressiom(net, X)
            R = zeros(size(X));
            for i = 1:length(X)
                if(X(i) >= net.theta)
                    R(i)=X(i);
                end
            end
        end
    end
    
    methods (Access=public)      
        function f1_init(net, S)
            net.f1_U = zeros(size(S));
            net.f1_P = zeros(size(S));
            net.f1_Q = zeros(size(S));
            net.f1_update_W(S);
            net.f1_update_X();
            net.f1_update_V();            
        end
        
        function f1_update(net, S, activeNeuronIdx)
            net.f1_update_U();
            net.f1_update_W(S);
            net.f1_update_P(activeNeuronIdx);
            net.f1_update_X();
            net.f1_update_Q();
            net.f1_update_V();            
        end
        
        function adjustWeights(net, activeNeuronIdx)
            net.f1f2(:,activeNeuronIdx) = net.f1_U' .* net.learningRate * net.d + net.f1f2(:,activeNeuronIdx) .* (1 + net.learningRate * net.d * (net.d - 1));
            net.f2f1(activeNeuronIdx,:) = net.f1_U  .* net.learningRate * net.d + net.f2f1(activeNeuronIdx,:) .* (1 + net.learningRate * net.d * (net.d - 1));
        end
        
        function learning(net, S, activeNeuronIdx)
            net.f1_update_W(S);
            net.f1_update_X();
            net.f1_update_Q();
            net.f1_update_V();    
            
            for i = 1:net.learningLength
                net.adjustWeights(activeNeuronIdx);
                net.f1_update(S, activeNeuronIdx);
            end
        end
    
        function normR = calculateVigilance(net, activeNeuronIdx)
            net.f1_update_U();
            net.f1_update_P(activeNeuronIdx);
                      
            R = (net.f1_U + net.f1_P .* net.c) ./ (net.e + norm(net.f1_U) + net.c * norm(net.f1_P));
            normR = norm(R);
        end
        
        function result = isActiveNeuronOK(net, activeNeuronIdx)
            result = net.calculateVigilance(activeNeuronIdx) + net.e >= net.vigilance;
            if(result)
                log4m.getLogger.trace('ART2::isActiveNeuronOK', sprintf('Class #%d OK', activeNeuronIdx));            
            else
                log4m.getLogger.trace('ART2::isActiveNeuronOK', sprintf('Class #%d FAILED', activeNeuronIdx));            
            end
        end
        
        function Y = simulateF1F2(net)
            Y = net.f1f2(:,1:net.nbOfClasses+1)'*net.f1_P';
            log4m.getLogger.trace('ART2::simulateF1F2', ['Y = ' sprintf('%f ', Y) ']']);
        end
        
        function activeNeuronIdx = determineClass(net, S, force)
            activeNeuronIdx = -1;
            net.f1_init(S);
            net.f1_update(S, activeNeuronIdx);
            
            Y = net.simulateF1F2();
            
            if (force)     
                maxCompValue = 0;
                for i=1:net.nbOfClasses
                    compValue = net.calculateVigilance(i) + net.e;
                    if (compValue > maxCompValue)
                        maxCompValue = compValue;
                    end
                end
                
                if (net.vigilance < maxCompValue)
                    net.vigilance = maxCompValue;
                end                 
                
                activeNeuronIdx = net.nbOfClasses+1;
                return;
            end

            while (true)
                [max_value,activeNeuronIdx] = max(Y);
                if(isnan(max_value))
                    error('!!! no capacity left');
                end

                if(net.isActiveNeuronOK(activeNeuronIdx))
                    return;
                else
                    Y(activeNeuronIdx)=NaN;
                end
            end
        end
                      
        function logItarationStatus(net, S)
            log4m.getLogger.debug('ART2::run', '***************************************************');
            log4m.getLogger.debug('ART2::run',['**************** Example ' sprintf('%.4d', net.iterationCounter) ' *********************']);
            log4m.getLogger.debug('ART2::run', '***************************************************');
            log4m.getLogger.debug('ART2::run', ['Source example = [ ' sprintf('%f ', S) ']']);
        end
        
        function [newcluster, activeNeuronIdx] = process(net, S, force)
            newcluster = false;
            activeNeuronIdx = net.determineClass(S, force);

            if (activeNeuronIdx == net.nbOfClasses+1)
                net.nbOfClasses = net.nbOfClasses+1;
                activeNeuronIdx = net.nbOfClasses;
                newcluster = true;
                log4m.getLogger.info(sprintf('Nb %d ART2::process', net.iterationCounter), sprintf('New class J = %d was added', activeNeuronIdx));                        
            else
                log4m.getLogger.debug(sprintf('Nb %d ART2::process', net.iterationCounter), sprintf('Class J = %d was determined', activeNeuronIdx));                        
            end

            net.learning(S, activeNeuronIdx);       
        end        
        
        function reset(net)
            net.f1f2 = net.f2f1WeightRatio .* ones(net.dataDim,net.capacity)./((1-net.d)*sqrt(net.dataDim));  
            net.f2f1 = zeros(net.capacity,net.dataDim);
            net.nbOfClasses = 0;
            net.iterationCounter = 0;
        end
        
        function net = Art2(art2params)
            net.capacity = art2params.capacity;
            net.a = art2params.a;
            net.b = art2params.b;
            net.c = art2params.c;
            net.d = art2params.d;
            net.e = art2params.e;
            net.learningLength = art2params.learningLength;
            net.learningRate = art2params.learningRate;
            net.dataDim = art2params.dataDim;
            net.theta = art2params.theta;
            net.vigilance = art2params.vigilance;
            net.f2f1WeightRatio = art2params.f2f1WeightRatio;
            net.reset();
        end
        
        function setupWeights(net, M)
            net.reset();
            nbOfPoints = size(M,1);
            learningRateBackup = net.learningRate;
            learningLengthBackup = net.learningLength;
            net.learningRate = 1;
            net.learningLength = 1;
            for pointNb = 1:nbOfPoints
                net.process(M(pointNb,:), true);
            end
            net.learningRate = learningRateBackup;
            net.learningLength = learningLengthBackup;
        end

        function [newcluster, clusterid] = processPoint(net, S)                   
            net.iterationCounter = net.iterationCounter + 1;
            net.logItarationStatus(S);
            
            [newcluster, clusterid] = process(net, S, false);

            if (log4m.getLogger.isDebugOn())
                net.printDebugInfo();            
            end
        end
        
        function printDebugInfo(net)
            logmessage = sprintf('f1f2 weights:');
            for i=1:net.nbOfClasses
                logmessage=strcat(logmessage, sprintf('\n#%d:',i), sprintf(' %f ',net.f1f2(:,i)));
            end
            log4m.getLogger.debug('ART2::DebugInfo', logmessage);
            logmessage = sprintf('f2f1 weights:');
            for i=1:net.nbOfClasses
                logmessage=strcat(logmessage, sprintf('\n#%d:',i), sprintf(' %f ',net.f2f1(i,:)));
            end
            log4m.getLogger.debug('ART2::DebugInfo', logmessage);
        end      
        
        function logDescription(net)
            desc = sprintf( ...
               ['Properties:\n' ...
                'a = %f, b = %f, c = %f, d = %f, e = %f\n' ...
                'learningRate = %f\n'...
                'learningLength = %d\n' ...
                'vigilance = %f\n' ...
                'theta = %f\n'], ...
                                  net.a, net.b, net.c, net.d, net.e, ...
                                  net.learningRate, net.learningLength, ...
                                  net.vigilance, net.theta);
            log4m.getLogger.info('Network description', desc);            
        end
    end
end
