%% Neighborhood Class 
% Defines all operations to manipulate solution objects

classdef Neighborhood < handle
    %%
    properties
        numOperation = 10;
        operation; % operation array: contains handles to all operations functions. size[numOperation]
    
    end % properties
    %% Static Methods
    methods(Static)
        %% Swap Harvet Units Method
        % Swap two harvest unit from the random farm at a random month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = swapHarvestUnits(s)
            numMonth = s.numMonthsToSchedule;
            %Choose randomly a month
            month = randi([1 numMonth]);
            numFarm = s.listMonthSchedule(month).numFarm;
            %Choose randomly a farm
            farm = randi([1 numFarm]);
            numHarvestUnit = s.listMonthSchedule(month).listFarm(farm).numHarvestUnit;
            if (numHarvestUnit == 1)
                w = -1;
                return;
            end %if
            %Choose randomly a harvest unit
            harvestUnit1 = randi([1 numHarvestUnit]);
            harvestUnit2 = harvestUnit1;
            %Choose randomly another harvest unit different from the first
            %one
            while(harvestUnit2 == harvestUnit1)
                harvestUnit2 = randi([1 numHarvestUnit]);
            end %while
            %Save the first harvest unit
            hu = s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(harvestUnit1);
            %Swap the harvest units 
            s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(harvestUnit1) = ...
            s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(harvestUnit2);
            s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(harvestUnit2) = hu;
            
            w = 1;
            %Calculate the cost
            
        end %function
        %% Swap N Harvet Units Method
        % Swap two harvest unit from the random farm at a random month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = swapNHarvestUnits(s)
            numMonth = s.numMonthsToSchedule;
            %Choose randomly a month
            month = randi([1 numMonth]);
            numFarm = s.listMonthSchedule(month).numFarm;
            %Choose randomly a farm
            farm = randi([1 numFarm]);
            numHarvestUnit = s.listMonthSchedule(month).listFarm(farm).numHarvestUnit;
            if(numHarvestUnit  < 4)
                [s,w] = Neighborhood.swapHarvestUnits(s);   
                return;
            end %if 
            % Get randomly number of harvest units to swap
            num = randi([1 floor(numHarvestUnit/2)]);
            % Get first and second block of harvest units 
            start1 = randi([1 numHarvestUnit]);
            start2 = randi([0 numHarvestUnit - 2*num]) + start1 + num;
            for h = 1:num
                hu = s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(h+start1-1,numHarvestUnit)+1);
                s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(h+start1-1,numHarvestUnit)+1) = ...
                     s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(h+start2-1,numHarvestUnit)+1);
                s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(h+start2-1,numHarvestUnit)+1) = hu;
            end %for
            % Calculate the cost
            
            w = 1;
        end %function
        %% Shift Harvest Unit Method
        % Shift a random number of harvest units from the random farm at a
        % random month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = shiftHarvestUnits(s)
            numMonth = s.numMonthsToSchedule;
            %Choose randomly a month
            month = randi([1 numMonth]);
            numFarm = s.listMonthSchedule(month).numFarm;
            %Choose randomly a farm
            farm = randi([1 numFarm]);
            numHarvestUnit = s.listMonthSchedule(month).listFarm(farm).numHarvestUnit;
            if (numHarvestUnit == 1)
                w = -1;
                return
            end %if
            % Choose randomly a start point, the displacement and the
            % number of harvest units to shift
            displacement = randi([1 numHarvestUnit-1]);
            start = randi([1 numHarvestUnit]);
            num = randi([1 (numHarvestUnit)]);
            
            for d = 1:displacement
                for h = 1:num 
                    hu = s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(num+1-h+start+d-1,numHarvestUnit)+1);
                    s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(num+1-h+start+d-1,numHarvestUnit)+1) = ...
                        s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(num+1-h+1+start+d-1,numHarvestUnit)+1);
                     s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(num+1-h+1+start+d-1,numHarvestUnit)+1) = hu;
                end %for
            end %for
            
            % Calculate Cost
            w = 1;
        end% function
        %% Reverse Harvest Unit Method
        % Reverse the order of a number of harvest units in a 
        % farm on a month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = reverseHarvestUnits(s)
            numMonth = s.numMonthsToSchedule;
            %Choose randomly a month
            month = randi([1 numMonth]);
            numFarm = s.listMonthSchedule(month).numFarm;
            %Choose randomly a farm
            farm = randi([1 numFarm]);
            numHarvestUnit = s.listMonthSchedule(month).listFarm(farm).numHarvestUnit;
            if (numHarvestUnit == 1)
                w = -1;
                return;
            end %if
            %Choose randomly number of harvest units to shift and a start
            %point
            start = randi([1 numHarvestUnit]);
            numShift = randi([2 numHarvestUnit]);
            for h = 1:floor(numShift/2)
                hu = s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(h+start-1,numHarvestUnit)+1);
                s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(h+start-1,numHarvestUnit)+1) = ...
                    s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(numShift+start-h+1-1,numHarvestUnit)+1);
                s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(mod(numShift+start-h+1-1,numHarvestUnit)+1) = hu;
            end %for
            
            % Calculate Cost
            w = 1;
        end% function
         %% Swap Farm Method
        % Swap two farms on a month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = swapFarms(s)
            numMonth = s.numMonthsToSchedule;
            %Choose randomly a month
            month = randi([1 numMonth]);
            numFarm = s.listMonthSchedule(month).numFarm;
            if( numFarm == 1)
                w = -1;
                return
            end %if
            %Choose randomly a farms
            farm1 = randi([1 numFarm]);
            farm2 = randi([1 numFarm]);
            while(farm2 == farm1)
                farm2 = randi([1 numFarm]);
            end %while
            %Save the first farm
            f = s.listMonthSchedule(month).listFarm(farm1);
            %Swap farms
            s.listMonthSchedule(month).listFarm(farm1) = s.listMonthSchedule(month).listFarm(farm2);
            s.listMonthSchedule(month).listFarm(farm2)= f;
            
            w =1;
        end% function
         %% Reverse Farm Method
        % Reverse the order of a number of farms in a month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = reverseFarms(s)
            numMonth = s.numMonthsToSchedule;
            %Choose randomly a month
            month = randi([1 numMonth]);
            numFarm = s.listMonthSchedule(month).numFarm;
            if( numFarm == 1)
                w = -1;
                return
            end %if
            %Choose randomly number of harvest units to shift and a start
            %point
            start = randi([1 numFarm]);
            numShift = randi([2 numFarm]);
            for h = 1:floor(numShift/2)
                f = s.listMonthSchedule(month).listFarm(mod(h+start-1,numFarm)+1);
                s.listMonthSchedule(month).listFarm(mod(h+start-1,numFarm)+1) = ...
                    s.listMonthSchedule(month).listFarm(mod(numShift+start-h+1-1,numFarm)+1);
                s.listMonthSchedule(month).listFarm(mod(numShift+start-h+1-1,numFarm)+1) = f;
            end %for
            % Calculate the cost
            w = 1;
        end% function
         %% Shift Farm Method
        % Reverse the order of a number of farms in a month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = shiftFarms(s)
            numMonth = s.numMonthsToSchedule;
            %Choose randomly a month
            month = randi([1 numMonth]);
            numFarm = s.listMonthSchedule(month).numFarm;
            if( numFarm == 1)
                w = -1;
                return
            end %if
            % Choose randomly a start point, the displacement and the
            % number of harvest units to shift
            displacement = randi([1 numFarm-1]);
            start = randi([1 numFarm]);
            num = randi([1 (numFarm)]);
            
            for d = 1:displacement
                for h = 1:num 
                    f = s.listMonthSchedule(month).listFarm(mod(num+1-h+start+d-1,numFarm)+1);
                    s.listMonthSchedule(month).listFarm(mod(num+1-h+start+d-1,numFarm)+1) = ...
                        s.listMonthSchedule(month).listFarm(mod(num+1-h+1+start+d-1,numFarm)+1);
                     s.listMonthSchedule(month).listFarm(mod(num+1-h+1+start+d-1,numFarm)+1) = f;
                end %for
            end %for
            
            %Calculate the cost
            w = 1;
        end% function
        %% Swap Harvest Units Between Months Method
        % Shift a random number of harvest units from the random farm at a
        % random month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = swapHarvestUnitsBetweenMonths(s)
            w = -1;
            disp('shiftHarvestUnit')
        end% function
        %% Swap N Harvest Units Between Months Method
        % Shift a random number of harvest units from the random farm at a
        % random month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = swapNHarvestUnitsBetweenMonths(s)
            w = -1;
            disp('shiftHarvestUnit')
        end% function
        %% Swap Farms Between Months Method
        % Shift a random number of harvest units from the random farm at a
        % random month
        % Inputs:
        %   s -> solution object to operate
        % Outputs:
        %   s -> resulting solution object
        function [s,w] = swapFarmsBetweenMonths(s)
            w = -1;
            disp('shiftHarvestUnit')
        end% function
    end % static methods
    %% Methods
    methods
        %% Create Operation Vector Method
        function createOperationVector(obj)
           % Create handles from all static methods
           h1 = @obj.swapHarvestUnits;
           h2 = @Neighborhood.swapNHarvestUnits;
           h3 = @Neighborhood.shiftHarvestUnits;
           h4 = @Neighborhood.reverseHarvestUnits;
           h5 = @Neighborhood.swapFarms;         
           h6 = @Neighborhood.reverseFarms;
           h7 = @Neighborhood.shiftFarms;
           h8 = @Neighborhood.swapHarvestUnitsBetweenMonths;
           h9 = @Neighborhood.swapNHarvestUnitsBetweenMonths;
           h10 = @Neighborhood.swapFarmsBetweenMonths;
           % create a cell array with all handle function
           obj.operation = {h1 h2 h3 h4 h5 h6 h7 h8 h9 h10};
           
        end 
        %% Get Handle Method
        
        function handle = getHandle(obj,index)
            handle = obj.operation{index};
        end
    end % methods
    
end