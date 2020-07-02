%% Solution Class
% classe que descreve uma solucao inteira -> composicao de varios
% MonthSchedule

% TODO:
% - showData()
% - construtor()
classdef Solution < handle
    %%
    properties  
        ID;                     % Solution ID
        numMonthsToSchedule;    % number of months to schedule   
        numFarm;                % number of farms
        listMonthSchedule;      % list of all month schedules in the planning horizon. size[numMonthsToSchedule]
        totalCost;              % the total cost of this solution  
        demandViability;        %
        capacityViability;      %     
        demand;                 % Demand matrix of all farms in the planning horizon. size[numFarms numTotalMonths]
        storage;                % Storage matarix of all farms inthe planning horizont. size[numFarms numTotalMonths]
 
        
    end %properties
    %%
    methods 
        %% Constructor Method
        % Construct a Solution object
        function obj = Solution(initialStorage, demand)
            if nargin > 0
                obj.demand = demand;
                obj.storage = zeros(size(obj.demand));
                obj.numFarm = length(initialStorage(:,1));
                for farm = 1:obj.numFarm
                    initial = initialStorage(farm);
                    for i = 1:3
                        if (initial > obj.demand(farm,i))
                            obj.storage(farm,i) = obj.demand(farm,i);
                            initial = initial - obj.demand(farm,i);
                        else
                            obj.storage(i) = initial;
                            %initial = 0;
                            break;
                        end % if              
                    end % for 
                end % for
            else
               obj.ID = 0;
               obj.numMonthsToSchedule = 0;
               obj.listMonthSchedule = MonthSchedule();
               obj.totalCost = 0; 
            end %if    
        end%function
        %% Show Data Method
        % Show all data in consoel
        function showData(obj) 
            fprintf('\nPrinting Solution %d Data\n', obj.ID);
            fprintf('\tTotal Cost: %f\n', obj.totalCost);
            fprintf('\tNumber of Months to Schedule: %d\n', obj.numMonthsToSchedule);
            for i = 1:obj.numMonthsToSchedule
               obj.listMonthSchedule(i).showData(); 
            end
        end%function
        
        %% Calculate Cost Method
        % calculate the total cost of this solution including all month
        % schedules
        function cost = calcCost(obj)
           obj.totalCost = 0;
           for i = 1: obj.numMonthsToSchedule - 1
              % sum cost of  month schedule
              obj.totalCost = obj.totalCost + (obj.listMonthSchedule(i).calcCost()); 
              % sum costs of transition to another month schedule      
              obj.totalCost = obj.totalCost + (obj.listMonthSchedule(i).calcTransitionCost(obj.listMonthSchedule(i+1))); 
           end % for
           % Sum cost of the last month schedule
           obj.totalCost = obj.totalCost + obj.listMonthSchedule(end).calcCost();
           cost = obj.totalCost;
        end%function
        
        %% Calculate Cost Parallel Method
        % calculate the total cost of this solution including all month
        % schedules using a parallel for function
        function cost = parCalcCost(obj)
            obj.totalCost = 0;
            cost = 0;
            listMonth = obj.listMonthSchedule;
            listMonth2 = obj.listMonthSchedule;
           parfor i = 1: obj.numMonthsToSchedule - 1
              % sum cost of  month schedule
              cost = cost + (listMonth(i).parCalcCost()); 
              % sum costs of transition to another month schedule      
              cost = cost + (listMonth(i).calcTransitionCost(listMonth2(i+1))); 
           end % parfor
           % Sum cost of the last month schedule
           obj.totalCost = cost + obj.listMonthSchedule(end).calcCost();
           cost = obj.totalCost;
        end %function
        
        %% Calculate Demand Viability Method
        % calcuate the total viability of this solution including all month
        % schedules
%         function e = calcViabilityDemand(obj) 
%             e = 0;
%             for i = 1:obj.numMonthsToSchedule
%                 month_e = obj.listMonthSchedule(i).calcViabilityDemand();
%                 if (month_e < 1 || month_e > 3) 
%                     e = -1;
%                     break
%                 else
%                     e = e + month_e/obj.numMonthsToSchedule;
%                 end %if       
%             end %for
%         end%function
        
        %% Calculate Month Demand Viability
        % calculate the viability of this month schedule
        % Inputs: 
        %   obj -> "This" pointer to use this specific object
        % Outputs:
        %   c   -> The capacity viability parameter 
%         function c =  calcViabilityCapacity(obj, obj_sd)
%             c = 0;
%             for i = 1:obj.numMonth
%                 c = c + obj.listMonthSchedule(i).calcViabilityCapacity(obj_sd); 
%             end %for
%         end%function
        
        %% Plot Solution Path 
        % This method plots in a graph the solution's path it's used by a
        % harvest machine
        % Inputs:
        %   obj     -> "This" pointer to use specific object 
        %   figure  -> Figure Object in which the path will be plotted        
        function plotPath(obj)
            %titleStr = strcat('Solution ',int2str(obj.ID), ' Path');
            %title(titleStr);          
            color = hsv(obj.numMonthsToSchedule);
            for m = 1:obj.numMonthsToSchedule-1
                figure;
                hold on;
                name = strcat('Month ', int2str(obj.listMonthSchedule(m).monthID));
                title(name);
                obj.listMonthSchedule(m).plotPath(color(m,:));
               % obj.listMonthSchedule(m).listFarm(end).listHarvestUnit(end).plotPath(obj.listMonthSchedule(m+1).listFarm(1).listHarvestUnit(1),color(m,:));
            end % for
            figure;
            hold on;
            name = strcat('Month ', int2str(obj.listMonthSchedule(1+m).monthID));
            title(name);
            obj.listMonthSchedule(end).plotPath(color(m,:));
        end %function
        %% Write Data To File 
        %Inputs:
        %   fileID: id from file to write
        function writeDataToFile(obj, fileID)
            fprintf(fileID,'\n*************************************************************************************\n');
            fprintf(fileID,'Printing Solution %d Data\n', obj.ID);
            fprintf(fileID,'\tTotal Cost: %f\n', obj.totalCost);
            fprintf(fileID,'\tNumber of Months to Schedule: %d\n', obj.numMonthsToSchedule);
            for i = 1:obj.numMonthsToSchedule
               obj.listMonthSchedule(i).writeDataToFile(fileID); 
            end % for
        end % function
        
        %% Update Storage Method
        % Inputs:
        %   month      -> Month ID the Wood Volume was obtained
        function updateStorage(obj,month)
            for farm = 1:obj.numFarm
                totalVolume = obj.listMonthSchedule(month).listFarm(farm).calcWoodVolume();
                farmIndex = obj.listMonthSchedule(month).listFarm(farm).id;
                for i = month+3:month+5
                    totalVolume = totalVolume + obj.storage(farmIndex,i);
                    if (totalVolume > obj.demand(farmIndex,i))
                        obj.storage(farmIndex,i) = obj.demand(farmIndex,i);
                        totalVolume = totalVolume - obj.demand(farmIndex,i);
                    else
                        obj.storage(farmIndex,i) = totalVolume;
                        break;
                    end %if                   
                end %for
            end %for
        end %function
    end %methods 
end %classdef