%% RCL Class 
% It implements the behavior of Restricted Candidate List of GRASP
% algorithm
classdef RCL  < handle
    %% Properties
    properties
        ID;             % the candidate list (it's an array of Harvest Unit ID's) 
        woodVolume;     % wood volume of each candidate
        farmID;         % the farm ID of each candidate
        cost;           % the cost of each candidate
        prob;           % the probability of chosing each candidate
        numCandidates;  % number of candidates
    end %properties
    %% Methods
    methods       
        %% Select Candidate Function
        %This function selects a candidate from the list using the
        % probabilites of RCL::prob 
        function [hu, f] = selectCandidate(obj)
             [i,j] = size(obj.prob);
             dim = 2;
             if (i > j) 
                 dim = 1;
             end             
             index = sum( bsxfun(@ge, rand(), cumsum(obj.prob./sum(obj.prob))),dim) + 1;  
             hu = obj.ID(index);
             f = obj.farmID(index);
        end %function        
        %% Filter RCL for Moving To Another Farm
        % Filter RCL using weight basis
        % if Last Farm
        function filterForNextFarm (obj, NVmax,criticalDemand,Vt,vcut, meanVelocity, lastFarm)
            index = zeros(length(obj.ID),1);
            NV = zeros(length(meanVelocity),1);
            delta_v = NVmax(lastFarm) - vcut;
            for f = 1:length(NVmax)
                NVmax(f) = NVmax(f) + delta_v*meanVelocity(f)/meanVelocity(lastFarm);
                if (NVmax(f) > criticalDemand(f) + Vt(f))
                    NV(f) = criticalDemand(f) + Vt(f);
                else
                    NV(f) = NVmax(f);
                end %if
            end %for
            
            for f = 1:length(NVmax)
                i1 = obj.woodVolume <= NV(f);
                i2 = obj.farmID == f;
                i = i1 & i2;
                if (max(i)== 0 && criticalDemand(f) > 0)
                    i = obj.woodVolume == min(obj.woodVolume(i2));
                end %if
                index = index | i;
            end % for
            obj.ID = obj.ID(index);
            obj.woodVolume = obj.woodVolume(index);
            obj.farmID = obj.farmID(index);
            obj.cost = obj.cost(index);
            obj.prob = obj.prob(index);
            obj.numCandidates = length(obj.ID);
        end %function
        %% Filter RCL for same far
        function filterForSameFar(obj, NV, criticalDemand , vcut)
            index = obj.woodVolume <= NV - vcut;
            if (max(index) == 0 && (criticalDemand - vcut) > 0)
                [~,index] = min(obj.woodVolume);
            end %if
                
            obj.ID = obj.ID(index);
            obj.woodVolume = obj.woodVolume(index);
            obj.farmID = obj.farmID(index);
            obj.cost = obj.cost(index);
            obj.prob = obj.prob(index);
            obj.numCandidates = length(obj.ID);
        end %function
        %% No Bias
         function b = noBias(obj)
            b = ones(length(obj.ID),1);
         end %function
         %% Linear Bias 
        function b = linearBias(obj)
            b = 1./obj.ID;
        end %function
        %% Create RCL candidates probabilities 
        function p = calcProb (obj,b)
            p = b./sum(b);
        end %function
    end %methods    
end %classdef