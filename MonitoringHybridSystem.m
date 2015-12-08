% Copyright (c) 2015, Mateusz Wójcik (mateuszjanwojcik@gmail.com)
% This is free software and it is distributed under the BSD license - see LICENSE file.
% If you use this code for your research, please cite the paper mentioned here:
%  https://github.com/mjwojcik/Art2MonitoringHybridSystem#reference

classdef MonitoringHybridSystem < handle 

    properties(Access=public)
        
        %% user input 
        DeltaTstart
        DeltaTcheckMinMax
        DeltaTcheckStability
        MinClusterCount
        MaxClusterCountRatio
        MinimalAreaSize
        OCC_MAX_N
        OCC_tresholdModifier
        VigilanceModifier
        DataDim
        EnableStereographicProjection
       
        %% internal members
        MaxClusterCountLastCheck
        MaxClusterCount
        ART2
        ART2prevNbOfClasses
        
        vOCC
        occCount
        
        T
        isinit
        
        scallingBounds
        
        history
        historyCount
        
        NextCheckMinMaxART2
        NextCheckStabilityART2
    end
    
    methods(Static)
    
        function [ R ] = measureMatrix( M )    
            s = size(M,2);
            R.vMinArr = zeros(1,s);
            R.vMaxArr = zeros(1,s);
            for i = 1:s
                V=M(:,i);
                R.vMinArr(i) = min(V);
                if (R.vMinArr(i) > 0)
                    R.vMinArr(i) = 0;
                end
                R.vMaxArr(i) = max(V);
            end
        end
        
        function [ R ] = scaleMatrix( M, scallingBounds )      
            R = M;
            s = size(M,2);
            for i = 1:s
                V = M(:,i);
                vMin = scallingBounds.vMinArr(i);
                vMax = scallingBounds.vMaxArr(i);
                R(:,i) = scaleVector(V, vMin, vMax);
            end
        end
        
        function [ R ] = transformUsingStereographicProjection( M )
            k = size(M,1);
            l = size(M,2);
            R = zeros(k,l+1);
            for i = 1:k
                V = M(i,:);
                R(i,:) = transformVector(V);
            end
        end
        
        function specialPoints = determineSpecialPoints(data, OCC_MAX_N)
            [~, vb_model, ~] = GaussianMixtureLib.vbgm(data', OCC_MAX_N);
            specialPoints = vb_model.m(:,vb_model.alpha > 1.01)';
        end
        
        function [ model ] = em_invcov( model )
            model.invcovs = model.Sigma;
            k = length(model.weight);
            for i=1:k
                model.invcovs(:,:,i) = inv(model.Sigma(:,:,i));
            end 
        end
        
        function R = em_f(X, model)
            [N,d] = size(X);
            k = length(model.weight);
            p = zeros(N,k);

            Z = (2*pi).^(-d/2);
            for i = 1:k
                dif = X - repmat(model.mu(:,i)',N,1);
                c = squeeze(model.invcovs(:,:,i));
                p(:,i) = sqrt(det(c))*Z*exp(-0.5*sum((dif*c).*dif,2));
            end

            R = sum(p.*repmat(model.weight,N,1), 2);
        end
        
        function [ hybridparams ] = getDefaultParams(dataDim)
            hybridparams.DeltaTstart = 500;
            hybridparams.DeltaTcheckMinMax = 50;
            hybridparams.DeltaTcheckStability = 200;
            hybridparams.MinClusterCount = 2;
            hybridparams.MaxClusterCountRatio = 1;
            hybridparams.MinimalAreaSize = 50;
            hybridparams.OCC_MAX_N = 30;
            hybridparams.OCC_tresholdModifier = 0.01;
            hybridparams.VigilanceModifier = 0.4;
            hybridparams.DataDim = dataDim;
            hybridparams.EnableStereographicProjection = true;
        end
    end   
    
    methods
                
        function log(obj, message)
            log4m.getLogger.info(sprintf('R[%.3f] A[%d] C[%d] NM[%d] NS[%d]', obj.ART2.vigilance, obj.occCount, obj.ART2.nbOfClasses, obj.NextCheckMinMaxART2, obj.NextCheckStabilityART2), message);
        end
        
        function ART2maintenanceMinMax(obj)
            if (obj.ART2.nbOfClasses < obj.MinClusterCount || (obj.ART2.nbOfClasses > obj.MaxClusterCount && obj.ART2.nbOfClasses ~= obj.MaxClusterCountLastCheck))
                if (obj.ART2.nbOfClasses < obj.MinClusterCount)
                    obj.ART2.vigilance = obj.ART2.vigilance + (1-obj.ART2.vigilance) * obj.VigilanceModifier;
                    obj.log('ART2maintenanceMinMax - vigilance increased');
                else
                    obj.ART2.vigilance = obj.ART2.vigilance - (1-obj.ART2.vigilance) * obj.VigilanceModifier;  
                    obj.MaxClusterCountLastCheck = obj.ART2.nbOfClasses;
                    obj.log('ART2maintenanceMinMax - vigilance decreased');
                end
                obj.NextCheckStabilityART2 = obj.NextCheckStabilityART2 + obj.DeltaTcheckMinMax;                
            end
            obj.NextCheckMinMaxART2 = obj.ART2.iterationCounter + obj.DeltaTcheckMinMax;           
        end
        
        function newarea = ART2maintenanceStability(obj)   
            if (obj.ART2.nbOfClasses == obj.ART2prevNbOfClasses)
                if (obj.historyCount >= obj.MinimalAreaSize)
                    obj.log(sprintf('ART2maintenanceStability - stable - creating new OCC from %d points', obj.historyCount));            
                    obj.addNewOCC();
                    newarea = true;                    
                else                    
                    obj.log(sprintf('ART2maintenanceStability - stable - but not enough points (%d)', obj.historyCount));            
                    newarea = false;
                end
            else
                obj.log('ART2maintenanceStability - unstable');            
                obj.ART2prevNbOfClasses = obj.ART2.nbOfClasses;
                newarea = false;
            end
            obj.NextCheckStabilityART2 = obj.T + obj.DeltaTcheckStability;
        end
        
        function addNewOCC(obj)
            obj.occCount = obj.occCount + 1;
            data = obj.history(1:obj.historyCount,:);
            occComponents = min(obj.ART2.nbOfClasses * 10, 100);
            [~, em_model,~] = GaussianMixtureLib.emgm(data',occComponents);
            em_model = obj.em_invcov(em_model);
            obj.vOCC(obj.occCount).em_model = em_model;
            obj.vOCC(obj.occCount).treshold = min(obj.em_f(data, em_model));
            obj.vOCC(obj.occCount).treshold = obj.OCC_tresholdModifier * obj.vOCC(obj.occCount).treshold;
                        
            obj.historyCount = 0;
            obj.ART2.reset();
            obj.NextCheckMinMaxART2 = obj.DeltaTcheckMinMax;
            obj.ART2prevNbOfClasses = 0;
        end
                
        function areaidx = isAlreadyKnown(obj, X)
            for i=1:obj.occCount
                if (obj.em_f(X, obj.vOCC(i).em_model) >= obj.vOCC(i).treshold)
                    areaidx = i;
                    return;
                end
            end
            areaidx = obj.occCount + 1;
        end
        
        
        function [newarea, newcluster, areaidx, clusteridx] = processPoint(obj, X)
            % newarea - is new area added when the point X is processed?
            %           (bool value)
            % newcluster - is new cluster added when the point X is
            %              processed?
            %              (bool value) 
            % areaidx - index of area which the point X belongs to
            % clusterid - index of ART-2 pattern which classified
            %              the point X
            
            
            obj.T = obj.T + 1;
            log4m.getLogger.setTime(obj.T);
            newarea = false;
            newcluster = false;
            areaidx = -1;
            clusteridx = -1;
            if (obj.isinit)      
                X = scaleVector(X, obj.scallingBounds.vMinArr, obj.scallingBounds.vMaxArr);
                areaidx = obj.isAlreadyKnown(X);
                if (areaidx <= obj.occCount)
                    log4m.getLogger.debug(sprintf('R[%.3f] A[%d] C[%d]', obj.ART2.vigilance, obj.occCount, obj.ART2.nbOfClasses), sprintf('processPoint - the point [%f] is known as area idx %d', X, areaidx));
                else
                    if (obj.EnableStereographicProjection)
                        S = obj.transformUsingStereographicProjection(X);
                    else
                        S = X;
                    end
                    [newcluster, clusteridx] = obj.ART2.processPoint(S);
                    obj.addToHistory(X);
                    
                    if (obj.ART2.iterationCounter >= obj.NextCheckMinMaxART2)
                        obj.ART2maintenanceMinMax();
                    end
                    if (obj.T >= obj.NextCheckStabilityART2)
                        newarea = obj.ART2maintenanceStability();
                    end
                end
            else
                obj.addToHistory(X);
                if (obj.T == obj.DeltaTstart)
                    obj.scallingBounds = obj.measureMatrix( obj.history );    
                    obj.history = obj.scaleMatrix(obj.history, obj.scallingBounds);
                    specialPoints = obj.determineSpecialPoints(obj.history, obj.OCC_MAX_N);
                    if (obj.EnableStereographicProjection)
                        specialPoints = obj.transformUsingStereographicProjection(specialPoints);
                    end
                                        
                    obj.ART2.setupWeights(specialPoints);
                    specialPointsCount = size(specialPoints,1);
                    obj.MaxClusterCount = floor(max(specialPointsCount*obj.MaxClusterCountRatio, obj.MinClusterCount));
                    obj.log(sprintf('determineSpecialPoints [%d, %d]', specialPointsCount, obj.MaxClusterCount));
                    obj.NextCheckMinMaxART2 = obj.DeltaTcheckMinMax;
                    obj.NextCheckStabilityART2 = obj.T + obj.DeltaTcheckStability;
                    obj.ART2prevNbOfClasses = 0;
                    obj.isinit = true;
                end
            end
        end
        
        function addToHistory(obj, X)
            if (size(obj.history,1)==obj.historyCount)
                obj.history=[obj.history ; zeros(obj.DeltaTcheckStability, obj.DataDim)];
            end
            obj.historyCount = obj.historyCount + 1;
            obj.history(obj.historyCount,:) = X; 
        end
        
        function obj = MonitoringHybridSystem(hybridparams, art2params)
            obj.DeltaTstart = hybridparams.DeltaTstart;
            obj.DeltaTcheckMinMax = hybridparams.DeltaTcheckMinMax;
            obj.DeltaTcheckStability = hybridparams.DeltaTcheckStability;
            obj.MinClusterCount = hybridparams.MinClusterCount;
            obj.MaxClusterCountRatio = hybridparams.MaxClusterCountRatio;
            obj.MinimalAreaSize = hybridparams.MinimalAreaSize;
            obj.OCC_MAX_N = hybridparams.OCC_MAX_N;
            obj.OCC_tresholdModifier = hybridparams.OCC_tresholdModifier;
            obj.VigilanceModifier = hybridparams.VigilanceModifier;
            obj.DataDim = hybridparams.DataDim;
            obj.EnableStereographicProjection = hybridparams.EnableStereographicProjection;

            obj.T = 0;
            obj.historyCount = 0;
            obj.occCount = 0;
            obj.MaxClusterCountLastCheck = -1;

            log4m.getLogger.setTime(obj.T);
            if (obj.EnableStereographicProjection)
                art2params.dataDim = obj.DataDim + 1;
            else
                art2params.dataDim = obj.DataDim;
            end
            obj.ART2 = Art2(art2params);            
            obj.history = zeros(obj.DeltaTstart, obj.DataDim);
        end    
    end
end

function [V] = scaleVector( V, vMin, vMax )    
    from = 0;
    to = 1/sqrt(2);

    V = V-vMin;
    if (vMax > 0)
        V = V./vMax;
    end
    V = V .* (to-from) + from;
end

function [ R ] = transformVector( V )
    l = length(V);
    R = zeros(l+1,1);
    kwadraty = 0;
    for i=1:l
        kwadraty = kwadraty + V(i)^2;
    end
    for i=1:l
        R(i) = 2*V(i)/(kwadraty+1);
    end
    R(l+1) = (1-kwadraty)/(1+kwadraty);
end
