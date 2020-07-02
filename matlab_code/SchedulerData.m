%% Scheduler Data Class
% Classe que contem importantes para o sequenciamento


% TODO:
% - atualizar showData()
classdef SchedulerData
    %%
    properties
        numMonthsToSchedule;        % number of months to schedule
        numTotalMonth;             % total number of months to
        numFarm;                   % number of farms in this scheduling problem
        numHarvestUnit;            % size[numFarms]
        eo;                         %operating efficiency. size[numMonths]
        dm;                         %machinary disponibility. size[numMonths]
        ds;                         % number of days in each month. size[numMonths]
        demand;                     % demand in each month. size[numMonths, numFarms]
        machTotalTime;              % The total time to operate the machinery considering both 8-hour shifts
        initialStorage;             % The initial total wood volume per farm. size[numFarms]
        
        max_k;
        max_vns_k;
        max_vnd_k;
        max_shake_operations;
        min_shake_operations;
        max_construct;
        min_construct;
        alfa;                       % Greedy/Randomized parameter
        epsilon;
        
    end%properties
    %%
    methods
        %% Constructor 
        function obj = SchedulerData(scData,numHU, machTotalTime, demand, initialStorage)
           
            obj.numFarm = length(numHU);
            
            obj.numHarvestUnit = zeros(obj.numFarm,1);
            
            obj.numHarvestUnit(:) = numHU(:);
            obj.numMonthsToSchedule = length(scData);            
            obj.numTotalMonth = obj.numMonthsToSchedule + 3;
            
            obj.eo = zeros(obj.numMonthsToSchedule,1);
            obj.dm = zeros(obj.numMonthsToSchedule,1);
            obj.ds = zeros(obj.numMonthsToSchedule,1);
            
            obj.eo(:) = scData(:,1);
            obj.dm(:) = scData(:,2);
            obj.ds(:) = scData(:,3);
            obj.machTotalTime = machTotalTime;       
            
            obj.demand = zeros(obj.numTotalMonth, obj.numFarm);
            obj.demand = demand(1:end,1:end);  
            
            obj.initialStorage = initialStorage;
        end % constructor
        %% Show Data Method
        function showData(obj)
            fprintf('Printing configuration data to Scheduler\n');
            fprintf('     Number of Farms: %d. Number of Months: %d\n', obj.numFarm, obj.numMonthsToSchedule);
            for i = 1:obj.numFarm
                fprintf('     Farm %d: %d Harvest Units\n',i, obj.numHarvestUnit(i));
            end %for
            fprintf('     Operating Efficiency: ');
            for i = 1:obj.numMonthsToSchedule
                fprintf('%d | ', obj.eo(i));                
            end %for
            fprintf('\n     Machinary Disponibility: ')
            for i = 1:obj.numMonthsToSchedule
                fprintf('%d | ', obj.dm(i));
            end %for
            fprintf('\n     Number of Days: ');
            for i = 1:obj.numMonthsToSchedule
                fprintf('%d | ', obj.ds(i));
            end %for
            fprintf('\n');
            fprintf('GRASP max iterations: %d\n', obj.max_k);
            fprintf('VNS max iterations: %d\n', obj.max_vns_k);
            fprintf('Shake max operations: %d\n', obj.max_shake_operations);
            fprintf('Shake min operations: %d\n', obj.min_shake_operations);
            fprintf('GAP max construct: %f\n', obj.max_construct);
            fprintf('GAP min construct: %f\n', obj.min_construct);
        end%function
        %% Show Data Method 2
        function showData2(obj)
            fprintf('Printing configuration data\n');
            fprintf('     Number of Farms: %d\n', obj.numFarm);
            fprintf('     Number of Months to Schedule: %d\n', obj.numMonthsToSchedule);
            fprintf('     Total Number of Months: %d\n', obj.numTotalMonth);           
        end%function
        
        %% Write Data To File Method
        
        function writeDataToFile(obj,fileID)
            fprintf(fileID,'Printing configuration data\n');
            fprintf(fileID,'     Number of Farms: %d\n', obj.numFarm);
            fprintf(fileID,'     Number of Months to Schedule: %d\n', obj.numMonthsToSchedule);
            fprintf(fileID,'     Total Number of Months: %d\n', obj.numTotalMonth);
            
            for i = 1:obj.numFarm
                fprintf(fileID,'     Farm %d: %d Harvest Units\n',i, obj.numHarvestUnit(i));
            end %for
            fprintf(fileID,'     Operating Efficiency: ');
            for i = 1:obj.numMonthsToSchedule
                fprintf(fileID,'%d | ', obj.eo(i));                
            end %for
            fprintf(fileID,'\n     Machinary Disponibility: ');
            for i = 1:obj.numMonthsToSchedule
                fprintf(fileID,'%d | ', obj.dm(i));
            end %for
            fprintf(fileID,'\n     Number of Days: ');
            for i = 1:obj.numMonthsToSchedule
                fprintf(fileID,'%d | ', obj.ds(i));
            end %for
            fprintf(fileID,'\n');
            
            fprintf(fileID,'GRASP max iterations: %d\n', obj.max_k);
            fprintf(fileID,'VNS max iterations: %d\n', obj.max_vns_k);
            fprintf(fileID,'Shake max operations: %d\n', obj.max_shake_operations);
            fprintf(fileID,'Shake min operations: %d\n', obj.min_shake_operations);
            fprintf(fileID,'GAP max construct: %f\n', obj.max_construct);
            fprintf(fileID,'GAP min construct: %f\n', obj.min_construct);  
        end %function
    end%methods
end %classdef