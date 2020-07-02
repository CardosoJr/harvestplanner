%% SolutionConstructer 
% this objects constructs initial solutions to the MTILS algorithm

% TODO: 
%   - escolher talhoes com dificuldade de corte menor com mais probabilidade
%   - deixar mais eficiente atribuicoes
%   - calcular apenas delta e
classdef SolutionConstructor 
    %% Static Methods
    methods(Static)
        %% IntialSolutionConstruct
        % Method to construct a viable initial solution
        % Inputs:
        %   SchedData -> SchedulerData object
        %   MS.listFarm     -> Array of Farm objects
        %   e_min     -> The viability demand parameter minimum the
        %   algorithm aims
        %   e_max     -> The viability demand parameter maximum the
        %   algorithm aims
        % Output:
        %   s         -> Viable initial solution ( Solution object)
        %   garbage   -> Array of Farm objects which includes all harvest
        %   units that were not selected

%         function [s, garbage] = InitialSolutionConstruct(SchedData, MS.listFarm, e_max, e_min)
%             
%             numFarm = SchedData.numFarm;
%             numMonth = SchedData.numMonthsToSchedule;
%             
%             % Construct solution object
%             s = Solution();
%             % Allocate array of MonthSchedule
%             s.listMonthSchedule = MonthSchedule.empty(numMonth,0);
%             % For all months to schedule
%             for month = 1:numMonth;
%                 fprintf('\nMES %d\n', month);
%                 accumulatedTime = 0;
%                 capacity = SchedData.machTotalTime * SchedData.dm(month)*SchedData.eo(month);
%                 % Allocate array of Farm object and intialize MonthSchedule
%                 % object
%                 s.listMonthSchedule(month).listFarm = Farm.empty(numFarm,0);
%                 s.listMonthSchedule(month).monthID = month;
%                 s.listMonthSchedule(month).numFarm = 0;
%                 
%                 farmVec = 1:numFarm; % Array of Farm indexes
%                 farmRelation = zeros(numFarm,1); % Array to relate new order of MS.listFarm to the original set
%                 f = 1;  % MS.listFarm index in listMonthSchedule 
%                 
%                 % guarantee an  1 <= e >= 3 to all MS.listFarm in this month
%                 while(f <= numFarm)
%                     
%                     % Randomly choose a farm
%                     [farm,index] = datasample(farmVec,1);
%                     farmVec(index) = []; % clear this farm from the original set to no choose it again
%                     
%                     s.listMonthSchedule(month).listFarm(f) = MS.listFarm(farm);
%                     s.listMonthSchedule(month).listFarm(f).listHarvestUnit = HarvestUnit.empty(MS.listFarm(farm).numHarvestUnit,0);
%                     s.listMonthSchedule(month).numFarm = numFarm;
%                     s.listMonthSchedule(month).listFarm(f).numHarvestUnit = 0;
%                     
%                     e = s.listMonthSchedule(month).listFarm(f).calcViabilityDemand(month);
%                     
%                     if (s.listMonthSchedule(month).listFarm(f).demand(month+3) == 0 || ...
%                           s.listMonthSchedule(month).listFarm(f).demand(month+4) == 0 || ...
%                           s.listMonthSchedule(month).listFarm(f).demand(month+5) == 0 )
%                       
%                         min = 1.5;
%                     else
%                         min = e_min;                      
%                     end %if
%                     
%                     % Here we have to allocate more harvest units to this farm
%                     % since it is not able to suply all 4 next months
%                     h = 1; % harvest unit index in listHarvestUnit
%                     if ( e > min)
%                         farmRelation(f) = -1;
%                     end %if
%                     while (e < min)            
%                         if (accumulatedTime > capacity)
%                             disp('A viable solution was not reached')
%                             return
%                         end %if
%                                             
%                         % uniformly randomly choose a harvest unit
%                         
%                         hu = randi([1 MS.listFarm(farm).numHarvestUnit]);
%                         de = MS.listFarm(farm).listHarvestUnit(hu).woodVolume/mean(MS.listFarm(farm).demand(month+3:month+5));
%                         aux = 1;
%                         min_de = de;
%                         hu_min = hu;
%                         if(e+de > 0.7*3)
%                             while (e+de > 0.7*3)
%                                 if ( aux >= MS.listFarm(farm).numHarvestUnit) 
%                                     break;
%                                 end % if
%                                 hu = randi([1 MS.listFarm(farm).numHarvestUnit]);
%                                 de = MS.listFarm(farm).listHarvestUnit(hu).woodVolume/mean(MS.listFarm(farm).demand(month+3:month+5));
%                                 aux = aux + 1;
%                                 if (min_de > de) 
%                                     min_de = de;
%                                     hu_min = hu;
%                                 end %if 
%                             end %while
%                             hu = hu_min;
%                         end%if
%                         
%                         fprintf('1 Fazenda %d: %d\n',farm ,MS.listFarm(farm).numHarvestUnit);
%                         s.listMonthSchedule(month).listFarm(f).listHarvestUnit(h) = MS.listFarm(farm).listHarvestUnit(hu);
%                         s.listMonthSchedule(month).listFarm(f).numHarvestUnit = s.listMonthSchedule(month).listFarm(f).numHarvestUnit + 1;
%                         % Recalculate the Viability Demand parameter
%                         e = s.listMonthSchedule(month).listFarm(f).calcViabilityDemand(month); 
%                         % Accumulate the time used to harvest
%                         accumulatedTime = accumulatedTime + MS.listFarm(farm).listHarvestUnit(hu).calcHarvestTime();
%                         % Clear this harvest unit from the original set
%                         MS.listFarm(farm).listHarvestUnit(hu) = []; 
%                         MS.listFarm(farm).numHarvestUnit = MS.listFarm(farm).numHarvestUnit - 1;   
%                         h = h + 1;
%                     end %while
%                     farmRelation(f) = farm;
%                     f = f + 1;
%                 end %while
%                 
%                f = 1;
%                farmVec = 1:numFarm; % Array of Farm indexes
%                % Enhance the harvesting
%                % Enhance all MS.listFarm and not exceed the capacity
%                while (f <= numFarm && accumulatedTime <= capacity )
%                    
%                    % Randomly choose a farm
%                     [farm,index] = datasample(farmVec,1);
%                     farmVec(index) = []; % clear this farm from the original set to no choose it again
%                     farmOriginal = farmRelation(farm);
%                     
%                     if (farmOriginal ~= -1)
%                         e = s.listMonthSchedule(month).listFarm(farm).calcViabilityDemand(month);
%                         % reach the maximum
%                         while (e < e_max && accumulatedTime <= capacity)      
% 
%                             % uniformly randomly choose a harvest unit
%                             hu = randi([1 MS.listFarm(farmOriginal).numHarvestUnit]);
%                             de = MS.listFarm(farmOriginal).listHarvestUnit(hu).woodVolume/mean(MS.listFarm(farmOriginal).demand(month+3:month+5));
%                             aux = 1;
%                             min_de = de;
%                             hu_min = hu;
%                             if(e+de > 0.7*3)
%                                 while (e+de > 0.7*3)
%                                     if ( aux >= MS.listFarm(farmOriginal).numHarvestUnit) 
%                                         break;
%                                     end % if
%                                     hu = randi([1 MS.listFarm(farmOriginal).numHarvestUnit]);
%                                     de = MS.listFarm(farmOriginal).listHarvestUnit(hu).woodVolume/mean(MS.listFarm(farmOriginal).demand(month+3:month+5));
%                                     aux = aux + 1;
%                                     if (min_de > de) 
%                                         min_de = de;
%                                         hu_min = hu;
%                                     end %if 
%                                 end %while
%                                 hu = hu_min;
%                             end%if                       
%                             fprintf('2 Fazenda %d: %d\n',farmOriginal,MS.listFarm(farmOriginal).numHarvestUnit);
%                             % Accumulate the time used to harvest
%                             timeToHarvest = MS.listFarm(farmOriginal).listHarvestUnit(hu).calcHarvestTime();
%                             accumulatedTime = accumulatedTime + timeToHarvest;
%                             if (accumulatedTime  <= capacity)
% 
%                                 s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(end+1) = MS.listFarm(farmOriginal).listHarvestUnit(hu);
%                                 s.listMonthSchedule(month).listFarm(farm).numHarvestUnit = s.listMonthSchedule(month).listFarm(farm).numHarvestUnit + 1;
%                                 % Recalculate the Viability Demand parameter
%                                 e = s.listMonthSchedule(month).listFarm(farm).calcViabilityDemand(month);
%                                 % Clear this harvest unit from the original set
%                                 MS.listFarm(farmOriginal).listHarvestUnit(hu) = []; 
%                                 MS.listFarm(farmOriginal).numHarvestUnit = MS.listFarm(farmOriginal).numHarvestUnit - 1;
%                             end %if
%                         end %while
%                         f = f + 1;
%                     end %if
%                end %while  
%                % Update the storage values
%                 for farm = 1:numFarm
%                   s.listMonthSchedule(month).listFarm(farm).updateStorage(month);
%                   index = farmRelation(farm);
%                   MS.listFarm(index).storage = s.listMonthSchedule(month).listFarm(farm).storage;
%                end %for
%             end %for
%            garbage = MS.listFarm;
%             
%         end %function
        %% Create completely Randomized Solution
        % Method to construct a viable initial solution
        % Inputs:
        %   SchedData -> SchedulerData object
        %   MS.listFarm     -> Array of Farm objects
        %   e_min     -> The viability demand parameter minimum the
        %   algorithm aims
        %   e_max     -> The viability demand parameter maximum the
        %   algorithm aims
        % Output:
        %   s         -> Viable initial solution ( Solution object)
        %   garbage   -> Array of Farm objects which includes all harvest
        %   units that were not selected
        
        function [s, garbage] = randomizedSolution(SchedData, MS)
            
            numFarm = SchedData.numFarm;
            numMonth = SchedData.numMonthsToSchedule;
            gap_max = SchedData.max_construct;
            gap_min = SchedData.min_construct;
            % Construct solution object
            s = Solution(SchedData.initialStorage,SchedData.demand);
            % Construct garbage object
            garbage = Solution();
            
            % Allocate array of MonthSchedule
            s.listMonthSchedule = MonthSchedule.empty(numMonth,0);
            s.numMonthsToSchedule = numMonth;
            % For all months to schedule
            for month = 1:numMonth;                
                accumulatedTime = 0;
                capacity = SchedData.machTotalTime * SchedData.dm(month)*SchedData.eo(month);
                % Allocate array of Farm object and intialize MonthSchedule
                % object
                s.listMonthSchedule(month).listFarm = Farm.empty(numFarm,0);
                s.listMonthSchedule(month).monthID = month;
                s.listMonthSchedule(month).numFarm = numFarm;
                
                farmVec = 1:numFarm; % Array of Farm indexes
                farmRelation = zeros(numFarm,1); % Array to relate new order of MS.listFarm to the original set
                
                % First Loop: This loop objective is to schedule harvest
                % units capable of suply the critical month demand
                for f = 1:numFarm
                    
                    % Randomly choose a farm
                    [farm,index] = datasample(farmVec,1);
                    farmVec(index) = []; % clear this farm from the original set to no choose it again
                    
                    s.listMonthSchedule(month).listFarm(f) = MS.listFarm(farm);
                    s.listMonthSchedule(month).listFarm(f).listHarvestUnit = HarvestUnit.empty(1,0);
                    s.listMonthSchedule(month).listFarm(f).numHarvestUnit = 0;
                    
                    % intially the cut wood volume is zero
                    v = 0;
                    % the demand on the critical month is the month demand
                    % minus the storage volume of that specific month
                    criticalDemand = s.demand(farm,month+3) - s.storage(farm,month+3);  
                                      
                    % Here we have to allocate more harvest units to this farm
                    % since it is not able to suply the critical month
                    h = 1; % harvest unit index in listHarvestUnit
                    if ( criticalDemand == 0)
                        farmRelation(f) = -1;
                    end %if
                    while (v < criticalDemand)            
                        if (accumulatedTime > capacity || MS.listFarm(farm).numHarvestUnit < 1)
                            disp('No viable solution was reached')
                            return
                        end %if
                                            
                        % uniformly randomly choose a harvest unit
                        hu_min = randi([1 MS.listFarm(farm).numHarvestUnit]);
                        dv_min = MS.listFarm(farm).listHarvestUnit(hu_min).woodVolume;

                        aux = 1;
                        while (v+dv_min > gap_min*criticalDemand && aux <= MS.listFarm(farm).numHarvestUnit)
                            hu = randi([1 MS.listFarm(farm).numHarvestUnit]);
                            dv = MS.listFarm(farm).listHarvestUnit(hu).woodVolume;
                            aux = aux + 1;
                            if (dv_min > dv) 
                                dv_min = dv;
                                hu_min = hu;
                            end %if 
                        end %while
                        
                        s.listMonthSchedule(month).listFarm(f).listHarvestUnit(h) = MS.listFarm(farm).listHarvestUnit(hu_min);
                        s.listMonthSchedule(month).listFarm(f).numHarvestUnit = s.listMonthSchedule(month).listFarm(f).numHarvestUnit + 1;
                        % Recalculate the volume tha was cut
                        v = v + dv_min; 
                        % Accumulate the time used to harvest
                        accumulatedTime = accumulatedTime + MS.listFarm(farm).listHarvestUnit(hu_min).calcHarvestTime();
                        % Clear this harvest unit from the original set
                        MS.listFarm(farm).listHarvestUnit(hu_min) = []; 
                        MS.listFarm(farm).numHarvestUnit = MS.listFarm(farm).numHarvestUnit - 1;   
                        h = h + 1;
                    end %while
                    farmRelation(f) = farm;
                end % for
                
               f = 1;
               farmVec = 1:numFarm; % Array of Farm indexes
               % Enhance the harvesting
               % Enhance all MS.listFarm and not exceed the capacity
               while (f <= numFarm && accumulatedTime <= capacity )
                   
                   % Randomly choose a farm
                    [farm,index] = datasample(farmVec,1);
                    farmVec(index) = []; % clear this farm from the original set to no choose it again
                    farmOriginal = farmRelation(farm);
                    
                    if (farmOriginal ~= -1)
                        v = s.listMonthSchedule(month).listFarm(farm).calcWoodVolume();
                        % The accumulated demand is the demand of all 3 usage
                         % months minus the total storage in this period
                        accumulatedDemand = sum(s.demand(farmOriginal,month+3:month+5)) - ...
                            sum(s.storage(farmOriginal,month+3:month+5));        
                        % reach the maximum
                        while (v <= gap_max*accumulatedDemand && accumulatedTime <= capacity && MS.listFarm(farmOriginal).numHarvestUnit > 0)      

                            % uniformly randomly choose a harvest unit
                            hu_min = randi([1 MS.listFarm(farmOriginal).numHarvestUnit]);
                            dv_min = MS.listFarm(farmOriginal).listHarvestUnit(hu_min).woodVolume;
                            aux = 1;             
                            while (v + dv_min > gap_max*accumulatedDemand)
                                if ( aux >= MS.listFarm(farmOriginal).numHarvestUnit) 
                                    break;
                                end % if
                                hu = randi([1 MS.listFarm(farmOriginal).numHarvestUnit]);
                                dv = MS.listFarm(farmOriginal).listHarvestUnit(hu).woodVolume;
                                aux = aux + 1;
                                if (dv_min > dv) 
                                    dv_min = dv;
                                    hu_min = hu;
                                end %if 
                            end %while
                            if(v+dv_min > accumulatedDemand)
                                break;
                            end
                            % Accumulate the time used to harvest
                            timeToHarvest = MS.listFarm(farmOriginal).listHarvestUnit(hu_min).calcHarvestTime();
                            accumulatedTime = accumulatedTime + timeToHarvest;
                            if (accumulatedTime  <= capacity)

                                s.listMonthSchedule(month).listFarm(farm).listHarvestUnit(end+1) = MS.listFarm(farmOriginal).listHarvestUnit(hu_min);
                                s.listMonthSchedule(month).listFarm(farm).numHarvestUnit = s.listMonthSchedule(month).listFarm(farm).numHarvestUnit + 1;
                                % Recalculate the Viability Demand parameter
                                v = v + dv_min;
                                % Clear this harvest unit from the original set
                                MS.listFarm(farmOriginal).listHarvestUnit(hu_min) = []; 
                                MS.listFarm(farmOriginal).numHarvestUnit = MS.listFarm(farmOriginal).numHarvestUnit - 1;
                            end %if
                        end %while
                        f = f + 1;
                        fprintf('Demanda: %f | Desperdicio: %f\n', accumulatedDemand, v - accumulatedDemand)
                    end %if
               end %while  
               % Update the storage values           
               s.updateStorage(month);
              
               % Check and delete MS.listFarm without harvest units
               f = 1;
               for  aux = 1:s.listMonthSchedule(month).numFarm
                    if (s.listMonthSchedule(month).listFarm(f).numHarvestUnit == 0)
                        s.listMonthSchedule(month).listFarm(f) = [];
                        s.listMonthSchedule(month).numFarm = s.listMonthSchedule(month).numFarm - 1;
                        f = f - 1;
                    end%if
                    f = f + 1;
               end%for
            end %for month
           % Check and delete MonthSchedules without MS.listFarm 
           month = 1; 
           for  aux = 1:numMonth
               if (s.listMonthSchedule(month).numFarm == 0)
                   s.listMonthSchedule(month) = [];
                   s.numMonthsToSchedule = s.numMonthsToSchedule - 1;
                   month = month - 1;
               end %if
               month = month + 1;
           end %for
           garbage = MS.listFarm;       
        end       
        %% Greedy and Randomized Solution Constructor 1 
        % This solution constructor implements the greedy algorithm from
        % GRASP metaheuristics
        % It implements the RCL limited by the quality of solutions
         function [s, garbage] = greedyAndRandomizedSolution1(SchedData, MS)
            
            numFarm = SchedData.numFarm;
            numMonth = SchedData.numMonthsToSchedule;
            alfa = SchedData.alfa;
            epsilon = SchedData.epsilon;
            gap_max = SchedData.max_construct;
            % Construct solution object
            s = Solution(SchedData.initialStorage,SchedData.demand);
            % Construct garbage object
            garbage = Solution();
            % Allocate array of MonthSchedule
            s.listMonthSchedule = MonthSchedule.empty(numMonth,0);
            s.numMonthsToSchedule = numMonth; 
            % Construct auxiliary objects
            
            % For all months to schedule
            for month = 1:numMonth   
                %Define capacity and total accumulated time to harvest variables
                accumulatedTime = 0;                
                capacity = SchedData.machTotalTime * SchedData.dm(month)*SchedData.eo(month);   
                % Allocate array of Farm object and intialize MonthSchedule
                % object
                s.listMonthSchedule(month).listFarm = Farm.empty(numFarm,0);
                s.listMonthSchedule(month).monthID = month;
                s.listMonthSchedule(month).numFarm = numFarm;
                MS.monthID = month;
                farmVec = 1:numFarm; % Array of Farm indexes
                % Randomly choose the first farm
                [farm,~] = datasample(farmVec,1);
                % Imaginary first harvest unit
                hu = 0;
                % if hu is zero, then the RCL is build randomly
                RCL = MS.listFarm(farm).buildRCLByQuality(alfa,hu,month);
                % Calculate NVmax parameters
                delta_capacity = capacity;           % the capacity left when all critical demands are suplied
                Vt = zeros(numFarm,1);               % Demand(t+4 to t+5) - Storage(t+4 to t+5) of each farm
                sumTVt = 0;                          % Sum of Vt's multiplied by meanVelocity of each farm
                criticalDemand = zeros (numFarm,1);  % Demand(t+3) - Storage(t+3) of each farm
                meanVelocity = zeros (numFarm,1);    % Mean cut velocity of each farm
                for f = 1:numFarm
                    delta_capacity = delta_capacity - ...
                        (s.demand(MS.listFarm(f).id,month+3) - s.storage(MS.listFarm(f).id,month+3)) / MS.listFarm(f).meanCutVelocity;
                    
                    Vt(f) = (sum(s.demand(f,month+4:month+5)) - sum(s.storage(f,month+4:month+5)));
                    sumTVt = sumTVt + Vt(f)/MS.listFarm(f).meanCutVelocity;
                    criticalDemand(f) = s.demand(f,month+3) - s.storage(f,month+3);   
                    meanVelocity(f) = MS.listFarm(f).meanCutVelocity;
                end %for farm 
                % if delta capacity is negative, there is no capacity to even suply the critical demand
                if (delta_capacity <= 0)
                    disp('No viable solution was reached');
                    return;
                end %if
                % Calculate the Maximum Needed Volume for each farm
                NVmax = zeros(numFarm,1);
                for f = 1:numFarm
                    NVmax(f) = (( Vt(f) /sumTVt) * delta_capacity) + criticalDemand(f);
                end %for    
                % Initialize variables
                v = 0;
                lastFarm = -1;
                %Select the first harvest unit 
                NV = NVmax(farm);
                % Needed Wood Volume cant be greater than the whole period demand
                if (NV > (criticalDemand(farm) + Vt(farm)))
                   NV = (criticalDemand(farm) + Vt(farm)); 
                end % if
                RCL.filterForSameFar(NV,criticalDemand(farm),v);
                [hu, ~] = RCL.selectCandidate(); 
                % Suply needed wood volume for each farm
                for f = 1:numFarm           
                    s.listMonthSchedule(month).listFarm(f) = MS.listFarm(farm);
                    s.listMonthSchedule(month).listFarm(f).listHarvestUnit = HarvestUnit.empty(1,0);
                    s.listMonthSchedule(month).listFarm(f).numHarvestUnit = 0;                    
                    % Use unused NVM from last farm
                    if (f > 1)
                        NVmax(farm) = (NVmax(lastFarm) - v)*meanVelocity(farm)/meanVelocity(lastFarm); % Use the time left of last iteration
                        % Calculate the Needed Volume (It has a safety gap from the Maximum Needed Volume
                        NV = NVmax(farm);
                        % Needed Wood Volume cant be greater than the whole period demand
                        if (NV > (criticalDemand(farm) + Vt(farm)))
                           NV = (criticalDemand(farm) + Vt(farm)); 
                        end % if
                    end % if                   
                    % Save farm ID for next iteration
                    lastFarm = farm;
                    % intially the cut wood volume is zero
                    v = 0;                   
                    % Here we have to allocate more harvest units to this farm since it is not able to suply the critical month
                    h = 1; % harvest unit index in listHarvestUnit         
                    while (v < NV)  
                        if (accumulatedTime > capacity || MS.listFarm(farm).numHarvestUnit < 1)
                            disp('No viable solution was reached')
                            return
                        end %if                        
                        dv = MS.listFarm(farm).listHarvestUnit(hu).woodVolume;                         
                        s.listMonthSchedule(month).listFarm(f).listHarvestUnit(h) = MS.listFarm(farm).listHarvestUnit(hu);
                        s.listMonthSchedule(month).listFarm(f).numHarvestUnit = s.listMonthSchedule(month).listFarm(f).numHarvestUnit + 1;
                        % Recalculate the volume tha was cut
                        v = v + dv; 
                        % Accumulate the time used to harvest
                        accumulatedTime = accumulatedTime + MS.listFarm(farm).listHarvestUnit(hu).calcHarvestTime();
                        % Create RCL to move to another harvest unit
                        RCL = MS.listFarm(farm).buildRCLByQuality(alfa,hu,month);
                        % Clear this harvest unit from the original set
                        MS.listFarm(farm).lastHarvestUnit = MS.listFarm(farm).listHarvestUnit(hu);
                        MS.listFarm(farm).listHarvestUnit(hu) = []; 
                        MS.listFarm(farm).numHarvestUnit = MS.listFarm(farm).numHarvestUnit - 1;
                        % Increase counter
                        h = h + 1;
                        % Select new harvest unit 
                        RCL.filterForSameFarm(NV,criticalDemand(farm),v);
                        hu = RCL.selectCandidate();
                    end %while
                    fprintf('Demanda Critica: %f | Demanda Rest.: %f | Desperdicio: %f\n',criticalDemand(farm), Vt(farm), v - (criticalDemand(farm) + Vt(farm)));             
                    % Clear from farmVec the current farm
                    farmVec(farmVec == farm) = [];
                    % Create RCL to move to another farm
                    RCL = MS.buildRCLByQuality(alfa,farm,farmVec);
                    RCL.filterForNextFarm(NVmax,criticalDemand,Vt,v,meanVelocity,lastFarm);
                    [hu, farm] = RCL.selectCandidate();
                end % for farm                
                % Update the storage values           
                s.updateStorage(month);
                % Check and delete MS.listFarm without harvest units
                f = 1;
                for  aux = 1:s.listMonthSchedule(month).numFarm
                    if (s.listMonthSchedule(month).listFarm(f).numHarvestUnit == 0)
                        s.listMonthSchedule(month).listFarm(f) = [];
                        s.listMonthSchedule(month).numFarm = s.listMonthSchedule(month).numFarm - 1;
                        f = f - 1;
                    end%if
                    f = f + 1;
                end%for                
            end % for month    
            % Check and delete MonthSchedules without MS.listFarm 
            month = 1; 
            for  aux = 1:numMonth
               if (s.listMonthSchedule(month).numFarm == 0)
                   s.listMonthSchedule(month) = [];
                   s.numMonthsToSchedule = s.numMonthsToSchedule - 1;
                   month = month - 1;
               end %if
               month = month + 1;
            end %for
            garbage = MS.listFarm;  
         end %function
        %% Solution Construct 4
        % This solution constructor implements the greedy algorithm from
        % GRASP metaheuristics
        % It implements the RCL limited by the quality of solutions using a
        % bias function to random them
         function [s, garbage] = greedyAndRandomizedSolution2(SchedData, MS, gap_max, gap_min, alfa)
         end %function 
    end %static methods
end %classdef