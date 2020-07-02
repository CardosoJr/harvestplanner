%% Month Schedule Class
% classe que descreve o sequenciamento de um mes

% TODO: 
% - construtor() 
% - showData()
classdef MonthSchedule
    %%
    properties 
        
        numFarm;        % Number of farms of this month schedule
        listFarm;       % List of farms in scheduled order. size [numFarm]
        monthID;        % Month
        
    end %properties
    %%
    methods 
        %% Constructor
        function obj = MonthSchedule()
            obj.numFarm = 0;
            obj.listFarm = Farm();
            obj.monthID = 0;
        end%function
        %% Calculate Cost Method
        % Calculate the total cost from this month schedule
        % Inputs: 
        %   obj   -> "This" pointer to use this specific object
        % Outputs:
        %   cost  -> The total cost of this month schedule
        function cost = calcCost(obj)
            cost = 0;
            if (obj.numFarm ~= 0)
                for i = 1:obj.numFarm - 1
                    cost = cost + obj.listFarm(i).calcCost(obj.monthID);
                    cost = cost + obj.listFarm(i).calcTransitionCost(obj.listFarm(i+1), obj.monthID);
                end %for
                cost = cost + obj.listFarm(end).calcCost(obj.monthID);
            end %if
        end %function
        
        %% Calculate Cost Parallel Method
        % Calculates the total cost from this month schedule using a
        % parallel for function
        function cost = parCalcCost(obj)
            cost = 0;
            if (obj.numFarm ~= 0)
                list = obj.listFarm;
                list2 = obj.listFarm;
                month = obj.monthID;
                parfor i = 1:obj.numFarm - 1
                    cost = cost + list(i).parCalcCost(month);
                    cost = cost + list(i).calcTransitionCost(list2(i+1), month);
                end%for
                cost = cost + obj.listFarm(end).calcCost(obj.monthID);
            end %if
        end %function
        %% Calculate Month Transition Cost 
        % Calculate the transition cost from the last Harvest Unit from
        % this month schedule to the first Harvest Unit of the next Month
        % Schedule
        % Inputs: 
        %   obj    -> "This" pointer to use this specific object
        %   obj_ms -> The MonthSchedule object of the next schedule in the
        %   solution
        % Outputs:
        %   cost   -> The total cost of this month schedule
        function cost = calcTransitionCost(obj, obj_ms) 
           cost = 0;
           if (obj.numFarm ~= 0 && obj_ms.numFarm ~= 0)
               cost = ((obj.listFarm(end).listHarvestUnit(end).xCoord - obj_ms.listFarm(1).listHarvestUnit(1).xCoord)^2 + ...
                   (obj.listFarm(end).listHarvestUnit(end).yCoord - obj_ms.listFarm(1).listHarvestUnit(1).yCoord)^2  )^0.5;

               cost = cost + obj_ms.listFarm(1).listHarvestUnit(1).cutDiff(obj_ms.monthID)*(( ...
                   obj_ms.listFarm(1).listHarvestUnit(1).xCoord - obj_ms.listFarm(1).xCarbCoord)^2 + ...
                   (obj_ms.listFarm(1).listHarvestUnit(1).yCoord - obj_ms.listFarm(1).yCarbCoord)^2)^0.5;
           end %if
        end   
        %% Calculate Month Demand Viability
        % calculate the viability of this month schedule
        % returns - 1 whether the month is not viable, or the avarege e of all
        % farms whether it is viable
        % Inputs: 
        %   obj -> "This" pointer to use this specific object
        % Outputs:
        %   e   -> The demand viability parameter 
%         function e =  calcViabilityDemand(obj)
%             e = 0;
%             for i = 1:obj.numFarm
%                 farm_e = obj.listFarm(i).calcViabilityDemand(obj.monthID);
%                 if ( e < 1 || e > 3) 
%                     e = -1; 
%                     break;
%                 else
%                     e = e + farm_e/obj.numFarm;
%                 end %if
%             end %for
%         end%function
        
        %% Calculate Month Capacity Viability
        % calculate the viability of this month schedule
        % Inputs: 
        %   obj    -> "This" pointer to use this specific object
        %   obj_sd -> SchedulerData object
         % Outputs:
        %   c   -> The capacity viability parameter 
%         function c =  calcViabilityCapacity(obj, obj_sd)
%             t = 0;
%             tmax = obj_sd.machTotalTime * obj_sd.dm(obj.monthID)*obj_sd.eo(obj.monthID);
%             for i = 1:obj.numFarm
%                 t = t + obj.listFarm(i).calcHarvestTime();
%             end %for
%             % The time to harvest all farms in this month exceeds the
%             % total available time
%             if (t > tmax) 
%                 c = t - tmax;
%             else
%                 c = 0;
%             end %if 
%         end %function
        %% Show Data       
        function showData(obj)
            fprintf('Printing Data from Month Schedule %d\n', obj.monthID);
            for i =1:obj.numFarm
                obj.listFarm(i).showData();
            end %for
        end %function       
        %% Plot MonthShedule Path
        function plotPath(obj, color)
            marker = {'o','s','x','d','v','h','.','+'};
            obj.listFarm(1).listHarvestUnit(1).plotCoordinates([0 0 0],marker(5),15);
            obj.listFarm(end).listHarvestUnit(end).plotCoordinates([0 0 0],marker(6),15);
            for f = 1:obj.numFarm - 1
                obj.listFarm(f).plotPath(color,marker(obj.listFarm(f).id));
                obj.listFarm(f).listHarvestUnit(end).plotPath(obj.listFarm(f+1).listHarvestUnit(1),color);
            end %for
            obj.listFarm(end).plotPath(color,marker(obj.listFarm(end).id))
        end %function     
        %% Write Data To File
        function writeDataToFile(obj, fileID)
            fprintf(fileID,'\n==================================================================\n');
            fprintf(fileID,'Printing Data from Month Schedule %d\n', obj.monthID);
            for i =1:obj.numFarm
                obj.listFarm(i).writeDataToFile(fileID);
            end %for
        end %function
       %% Build RCL by size aspect
        % This method builds a RCL list of a specific Harvest Unit
        % condidering the next harvest unit will be in another Farm.
        % This RCL is limited by size 'p'
        function rcl = buildRCLBySize (obj, p, huID, farm, farmVec)
            C = zeros(sum([obj.listFarm(farmVec).numHarvestUnit]),4);
            j = 1;
            for index = farmVec
               f = obj.listFarm(index);
               for hu = 1:f.numHarvestUnit
                    C(j,1) = f.listHarvestUnit(hu).calcCost() + obj.listFarm(farm).listHarvestUnit(huID).calcTransitionCost(f.listHarvestUnit(hu),...
                        [f.xCarbCoord f.yCarbCoord],obj.monthID);
                    C(j,2) = hu;
                    C(j,3) = f.id;
                    C(j,4) = f.listHarvestUnit(hu).woodVolume;
                    j = j + 1;
               end %for                
            end %for
            [C_hat, I] = sort(C(:,1));
            % build RCL object
            rcl = RCL();
            rcl.cost = C_hat(1:p,1);
            rcl.ID = C(I(1:p),2);
            rcl.farmID =C(I(1:p),3);
            rcl.woodVolume = C(I(1:p),4);
            rcl.numCandidates = p;
            % Calculate each candidate probability using LINEAR BIAS
            rcl.prob = rcl.calcProb(rcl.linearBias());
        end %function        
        %% Build RCL by quality aspect
        % This method builds a RCL list of a specific Harvest Unit
        % condidering the next harvest unit will be in another Farm.
        % This RCL is limited by greedy function values between
        % [c_min, c_min + alfa * (c_max - c_min)]
        function rcl = buildRCLByQuality (obj,alfa, farm, farmVec)
            C = zeros(sum([obj.listFarm(farmVec).numHarvestUnit]),4);
            j = 1;
            for index = farmVec
               f = obj.listFarm(index);
               for hu = 1:f.numHarvestUnit
                    C(j,1) = f.listHarvestUnit(hu).calcCost() + obj.listFarm(farm).lastHarvestUnit.calcTransitionCost(f.listHarvestUnit(hu),...
                        [f.xCarbCoord f.yCarbCoord],obj.monthID);
                    C(j,2) = hu;
                    C(j,3) = f.id;
                    C(j,4) = f.listHarvestUnit(hu).woodVolume;
                    j = j + 1;
               end %for    
            end %for
            [C_hat, I] = sort(C(:,1));
            c_min = min(C(:,1));
            c_max = max(C(:,1));
            i = C_hat <= c_min + alfa * (c_max - c_min);
            % build RCL object
            rcl = RCL();
            rcl.cost = C_hat(i,1);
            rcl.ID = C(I(i),2);
            rcl.farmID = C(I(i),3);
            rcl.woodVolume = C(I(i),4);
            rcl.numCandidates = length(rcl.cost);
            % Calculate each candidate probability using LINEAR BIAS
            rcl.prob = rcl.calcProb(rcl.linearBias());
        end %function   
    end %methods 
end %classdef