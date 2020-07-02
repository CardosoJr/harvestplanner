%% Farm Class
% classe que descreve a fazenda

% TODO:
% atualizar metodo showData
classdef Farm
    %%
    properties
        
        numHarvestUnit;   % number of harvest units
        listHarvestUnit;  % list of harvest units. size [numHarvestUnit]
        name;             % farm name
        id;               % farm id 
        xCarbCoord;       % x Coordinates from Carbonization Center 
        yCarbCoord;       % y Coordinates from Carbonization Center
        lastRCL;          % last RCL calculate
        meanCutVelocity;  % Mean Cut Velocity for all harvest units
        lastHarvestUnit;
        
    end %properties
    %% Public Methods
    methods
        %% Constructor
        % Construct a object Farm
        function obj = Farm(name,id, harvestUnits, CarbCoord)
            
            if nargin > 0 
                obj.meanCutVelocity = 0;
                obj.name = name;
                obj.id = id;
                obj.numHarvestUnit = length(harvestUnits);
                obj.listHarvestUnit = HarvestUnit.empty(obj.numHarvestUnit,0);
                for i = 1:obj.numHarvestUnit
                    obj.listHarvestUnit(i) = harvestUnits(i);
                    obj.meanCutVelocity = obj.meanCutVelocity + harvestUnits(i).velocity;
                end %for 
                obj.meanCutVelocity = obj.meanCutVelocity/obj.numHarvestUnit;
                obj.xCarbCoord = CarbCoord(1);
                obj.yCarbCoord = CarbCoord(2);
            else 
                obj.listHarvestUnit = HarvestUnit();
                obj.name = '';
                obj.id = 0;
            end % if                       
        end%function
        %% Show Data Method 
        % Show all data from this farm in the console
        function showData(obj)
           fprintf('Printing Data from Farm %s; Id: %d\n', obj.name, obj.id);
           fprintf('\t Carbonization Center Coordinates: [%f %f]\n', obj.xCarbCoord, obj.yCarbCoord);
           fprintf('\n');
           fprintf('Printing data from its Harvest Units\n');
           for i =1:obj.numHarvestUnit
                obj.listHarvestUnit(i).showData();
           end %for
        end%function
        %% Calculate Cost Method
        % Calculates the total cost to harvest every harvest unit in this
        % farm following the order of list Harvest Unit
        % Inputs:
        %   obj   -> "This" pointer to use this specific object
        %   month -> Month id of this particular farm schedule
        % Outputs:
        %   cost -> The total cost to harvest this farm
        function cost = calcCost(obj, month) 
            cost = 0;
            for i = 1:obj.numHarvestUnit-1
                cost = cost + obj.listHarvestUnit(i).calcCost();
                cost = cost + obj.listHarvestUnit(i).calcTransitionCost(obj.listHarvestUnit(i+1), [obj.xCarbCoord obj.yCarbCoord], month);
            end%for
            i = i + 1;
            cost = cost + obj.listHarvestUnit(i).calcCost();
        end%function
        
        %% Calculate Cost Parallel Method 
        % Calculates the total cost, using a parallel for function, to harvest every harvest unit in this
        % farm following the order of list Harvest Unit
         function cost = parCalcCost(obj, month) 
            cost = 0;
            list1 = obj.listHarvestUnit;
            list2 = obj.listHarvestUnit;
            xCarb = obj.xCarbCoord;
            yCarb = obj.yCarbCoord;
            parfor i = 1:obj.numHarvestUnit-1
                cost = cost + list1(i).calcCost();
                cost = cost + list1(i).calcTransitionCost(list2(i+1), [xCarb yCarb], month);
            end %parfor
            cost = cost + obj.listHarvestUnit(end).calcCost();
        end%function
        
        %% Calculate Transition Cost to Another Farm Method
        % Calculates the total cost of moving the machinary from this farm
        % to another
        % Inputs:
        %   obj    -> "This" pointer to use this specific object
        %   obj_fa -> Farm object of the next farm on the month schedule
        %   month  -> Month id of the next farm schedule
        % Outputs:
        %   cost -> The transition cost
        function cost = calcTransitionCost(obj, obj_fa, month)
             cost = ((obj.listHarvestUnit(end).xCoord - obj_fa.listHarvestUnit(1).xCoord)^2 + ...
             (obj.listHarvestUnit(end).yCoord - obj_fa.listHarvestUnit(1).yCoord)^2)^0.5;
             cost = cost + obj_fa.listHarvestUnit(1).cutDiff(month)*((obj_fa.xCarbCoord - obj_fa.listHarvestUnit(1).xCoord)^2 + ...
             (obj_fa.yCarbCoord - obj_fa.listHarvestUnit(1).yCoord)^2)^0.5;
                        
        end %function   
        %% Calculate the total cut wood volume in all harvest units
        function volume = calcWoodVolume(obj)
            volume = 0;
            for i = 1:obj.numHarvestUnit
                volume = volume + (obj.listHarvestUnit(i).woodVolume);
            end % for
        end % function 
        
        %% Calculate Farm Shedule Viability on Demand Perspective
        % This method calculates the viability to schedule this farm using
        % all harvest units in listHarvestUnit
        % The equation to calculate viability is 
        % e = (V_cut + Storage)/mean(D) in the period of wood use
        % Inputs:
        %   obj   -> "This" pointer to use this specific object
        %   month -> month id of this particular farm schedule 
        % Outputs:
        %   e -> The demand viability parameter
%         function e = calcViabilityDemand(obj, month) 
%             
%             V_cut = obj.calcWoodVolume();
%             meanD = sum(obj.demand(month+3:month+5))/3;
%             meanE = sum(obj.storage(month+3:month+5))/3;
%             if (meanD ~= 0)
%                 e = (V_cut + meanE)/meanD;
%             else
%                 e = 100;
%             end % if
%         end % function     
        %% Calculate Time To Harvest
        % Inputs: 
        %   obj    -> "This" pointer to use this specific object
        % Outputs:
        %   t -> The total time to harvest this farm
        function t = calcHarvestTime(obj)
            t = 0;            
            for i  = 1:obj.numHarvestUnit
                t = t + obj.listHarvestUnit(i).calcHarvestTime();
            end %for
        end % function       
        %% Plot Farm Paths
        function plotPath(obj, color, marker)
            for hu = 1:obj.numHarvestUnit - 1
                obj.listHarvestUnit(hu).plotCoordinates(color,marker,5);
                obj.listHarvestUnit(hu).plotPath(obj.listHarvestUnit(hu+1),color);
            end %for
            obj.listHarvestUnit(end).plotCoordinates(color,marker,5);
        end %function       
        %% Write Data To File    
        function writeDataToFile(obj,fileID)
           fprintf(fileID,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
           fprintf(fileID,'Printing Data from Farm %s; Id: %d\n', obj.name, obj.id);
           fprintf(fileID,'\t Carbonization Center Coordinates: [%f %f]\n', obj.xCarbCoord, obj.yCarbCoord);
           fprintf(fileID,'\n');
           fprintf(fileID,'Printing data from its Harvest Units\n');
           for i =1:obj.numHarvestUnit
                obj.listHarvestUnit(i).writeDataToFile(fileID);
           end %for
        end %function       
        %% Build RCL by size aspect
        % This method builds a RCL list of a specific Harvest Unit
        % condidering the next harvest unit will be in the same Farm.
        % This RCL is limited by size 'p'
        function rcl = buildRCLBySize (obj, p, huID, month)           
            if huID > 0 
                C = zeros(obj.numHarvestUnit-1,2);
                j = 1;
                for i = 1:huID-1
                    C(j,1) = obj.listHarvestUnit(i).calcCost() + obj.listHarvestUnit(huID).calcTransitionCost(obj.listHarvestUnit(i),...
                        [obj.xCarbCoord obj.yCarbCoord],month);
                    C(j,2) = obj.listHarvestUnit(i).woodVolume;
                    j = j + 1;
                end %for
                for i = huID+1:obj.numHarvestUnit
                    C(j,1) = obj.listHarvestUnit(i).calcCost() + obj.listHarvestUnit(huID).calcTransitionCost(obj.listHarvestUnit(i),...
                        [obj.xCarbCoord obj.yCarbCoord],month);
                    C(j,2) = obj.listHarvestUnit(i).woodVolume;
                    j = j + 1;
                end %for
                [C_hat, I] = sort(C(:,1));
                % build RCL object
                rcl = RCL();
                rcl.cost = C_hat(1:p);
                rcl.ID = I(1:p);
                rcl.farmID = obj.id*ones(p,1);
                rcl.woodVolume = C(I(1:p),2);
                rcl.numCandidates = p;
                % Calculate each candidate probability using LINEAR BIAS
                rcl.prob = rcl.calcProb(rcl.linearBias());               
            else
                rcl = obj.buildInitialRCL();
            end %if
        end %function        
        %% Build RCL by quality aspect
        % This method builds a RCL list of a specific Harvest Unit
        % condidering the next harvest unit will be in the same Farm.
        % This RCL is limited by greedy function values between
        % [c_min, c_min + alfa * (c_max - c_min)]
        function rcl = buildRCLByQuality (obj,alfa, huID, month)
            if huID > 0 
                C = zeros(obj.numHarvestUnit-1,2);
                j = 1;
                for i = 1:huID-1
                    C(j,1) = obj.listHarvestUnit(i).calcCost() + obj.listHarvestUnit(huID).calcTransitionCost(obj.listHarvestUnit(i),...
                        [obj.xCarbCoord obj.yCarbCoord],month);
                    C(j,2) = obj.listHarvestUnit(i).woodVolume;
                    j = j + 1;
                end %for
                for i = huID+1:obj.numHarvestUnit
                    C(j,1) = obj.listHarvestUnit(i).calcCost() + obj.listHarvestUnit(huID).calcTransitionCost(obj.listHarvestUnit(i),...
                        [obj.xCarbCoord obj.yCarbCoord],month);
                    C(j,2) = obj.listHarvestUnit(i).woodVolume;
                    j = j + 1;
                end %for
                [C_hat, I] = sort(C(:,1));
                c_min = min(C(:,1));
                c_max = max(C(:,1));
                i = C_hat <= c_min + alfa * (c_max - c_min);
                % build RCL object
                rcl = RCL();
                rcl.cost = C_hat(i,1);
                rcl.ID = I(i);
                rcl.farmID = obj.id*ones(length(rcl.cost),1);
                rcl.woodVolume = C(I(i),2);
                rcl.numCandidates = length(rcl.cost);
                % Calculate each candidate probability using LINEAR BIAS
                rcl.prob = rcl.calcProb(rcl.linearBias());
                
            else
                rcl = obj.buildInitialRCL();
            end % if
        end %function      
        %% Build Initial RCL 
        % This method builds a initial random RCL
        % Include alfa ?
        function rcl = buildInitialRCL(obj)
            C = zeros(obj.numHarvestUnit,2);
            for hu = 1:obj.numHarvestUnit
                C(hu,1) = rand();
                C(hu,2) = obj.listHarvestUnit(hu).woodVolum;
            end %for
            [C_hat, I] = sort(C(:,1)); 
            % build RCL object
            rcl = RCL();
            rcl.cost = C_hat;
            rcl.ID = I;
            rcl.farmID = obj.id*ones(obj.numHarvestUnit,1);
            rcl.woodVolume = C(I,2);
            rcl.numCandidates = obj.numHarvestUnit;
            % Calculate each candidate probability using LINEAR BIAS
            rcl.prob = rcl.calcProb(rcl.linearBias());
        end %function        
    end %methods
end %classdef