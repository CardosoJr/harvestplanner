%% GRASP
% Implements the GRASP metaheuristc and its variants
classdef GRASP < handle
    %% Properties
    properties
        alfaArray;             % Array with all possible alfas to be used in reactiveGrasp
        avgSolution;          % Contain the average from all solutions using each alfa from alfaArray. Size [alfaArray]
        alfaProb;             % Probabilities for each alfa to be used. Size [alfaArray]
        
    end %properties
    %% Public Methods
    methods
        %% Constructor
        function obj = GRASP()
            obj.alfaArray = 0;
            obj.avgSolution = 0;
            obj.alfaProb = 0;
        end %function
        %% Define Alfa Array Method
        function defAlfaArray(obj, alfa)
            obj.alfaArray = alfa;
            obj.avgSolution = zeros(length(alfa));
            obj.alfaProb = ones(length(alfa))/length(alfa);
        end %function
    end %methods
    %% Static Methods
    methods (Static)
        %% Normal GRASP
        function [Best_Solution, fobj] = normalGRASP(MS, schedData)
            vns = VNS();
            fobj = zeros(schedData.max_k,1);
            % Find the best solution
            for k = 1:schedData.max_k  
                 % Create Initial Solution
                 [x_hat, garbage] = SolutionConstructor.randomizedSolution(schedData, MS,schedData.max_construct,schedData.min_construct);
                 x_hat.ID = k;
                 x_hat.writeDataToFile(fileID);

                 fprintf('\n%d Solution Construct\nCalling GeneralVNS method\n', k);
                 x_hat = vns.GeneralVNS(x_hat,schedData.max_vnd_k,schedData.max_vns_k,schedData.max_shake_operations,schedData.min_shake_operations);
                 if (x_hat.calcCost() < best_solution.totalCost)
                    Best_Solution = x_hat;
                 end %if
                 fprintf('~~~~~~~~~~ Best Solution Cost: %f ~~~~~~~~~~\n', Best_Solution.totalCost); 
                 fobj(k) = Best_Solution.totalCost;
            end %for

        end %function
        %% DEBBUG GRASP
        function [Best_Solution, fobj] = debbugGRASP(MS, schedData)
            vns = VNS();
            fobj = zeros(schedData.max_k,1);
            % Define Test Parameters and File
            fprintf('Creating Test File...\n');
            c = clock;
            outputFile = strcat ('TESTE_',int2str(c(3)),'_',int2str(c(2)),'_',int2str(c(1)),'_',int2str(c(4)),'_',int2str(c(5)),'_',int2str(c(6)),'.txt');
            fileID = fopen(outputFile, 'at');
            
            % Find the best solution
            for k = 1:schedData.max_k  
                 % Create Initial Solution
                 [x_hat, garbage] = SolutionConstructor.randomizedSolution(schedData, MS,schedData.max_construct,schedData.min_construct);
                 x_hat.ID = k;
                 x_hat.writeDataToFile(fileID);

                 fprintf('\n%d Solution Construct\nCalling GeneralVNS method\n', k);
                 x_hat = vns.GeneralVNS(x_hat,schedData.max_vnd_k,schedData.max_vns_k,schedData.max_shake_operations,schedData.min_shake_operations);
                 if (x_hat.calcCost() < Best_Solution.totalCost)
                    Best_Solution = x_hat;
                 end %if
                 fprintf('~~~~~~~~~~ Best Solution Cost: %f ~~~~~~~~~~\n', Best_Solution.totalCost);  
                 fobj(k) = Best_Solution.totalCost;
            end %for
            fclose(fileID);

        end %function       
        %% Parallel GRASP
        function [Best_Solution, fobj] = parallelGRASP(MS, schedData)
            vns = VNS();
            fobj = zeros(schedData.max_k,1);
            % Preparing Parallel Programming Configurations
            fprintf('Getting ready to start parallel pool ... \n');
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool  
            end
            
            % Find the best solution
            for k = 1:schedData.max_k  
                 % Create Initial Solution
                 [x_hat, garbage] = SolutionConstructor.randomizedSolution(schedData, MS,schedData.max_construct,schedData.min_construct);
                 x_hat.ID = k;
                 x_hat.writeDataToFile(fileID);

                 fprintf('\n%d Solution Construct\nCalling GeneralVNS method\n', k);
                 x_hat = vns.GeneralVNS(x_hat,schedData.max_vnd_k,schedData.max_vns_k,schedData.max_shake_operations,schedData.min_shake_operations);
                 if (x_hat.calcCost() < Best_Solution.totalCost)
                    Best_Solution = x_hat;
                 end %if
                 fprintf('~~~~~~~~~~ Best Solution Cost: %f ~~~~~~~~~~\n', Best_Solution.totalCost);    
                 fobj(k) = Best_Solution.totalCost;
            end %for          
        end %function
        %% Reactive GRASP
        function [Best_Solution, fobj] = reactiveGRASP(MS, schedData)
            vns = VNS();
            fobj = zeros(schedData.max_k,1);
            % Find the best solution
            for k = 1:schedData.max_k  
                 % Create Initial Solution
                 [x_hat, garbage] = SolutionConstructor.randomizedSolution(schedData, MS,schedData.max_construct,schedData.min_construct);
                 x_hat.ID = k;
                 x_hat.writeDataToFile(fileID);

                 fprintf('\n%d Solution Construct\nCalling GeneralVNS method\n', k);
                 x_hat = vns.GeneralVNS(x_hat,schedData.max_vnd_k,schedData.max_vns_k,schedData.max_shake_operations,schedData.min_shake_operations);
                 if (x_hat.calcCost() < Best_Solution.totalCost)
                    Best_Solution = x_hat;
                 end %if
                 fprintf('~~~~~~~~~~ Best Solution Cost: %f ~~~~~~~~~~\n', Best_Solution.totalCost); 
                 fobj(k) = Best_Solution.totalCost;
            end %for
        end %function      
    end %methods 
    %% Private methods
    methods (Access = private)
        %% Update avgSolution
        function updateAvgSolution(obj,s)
            i = obj.alfaArray == alfa;
            obj.avgSolution(i) = (obj.avgSolution(i) + s.totalCost)/2.0;
        end %function
        %% Reactive GRASP: Update Alfa
        function updateAlfa(obj, s, alfa)
            obj.updateAvgSolution(s);
            aux = obj.avgSolution == 0;
            if(isempty(obj.avgSolution(aux)))     
                for i = 1:length(obj.avgSolution)
                    obj.alfaProb(i) = (s.totalCost/obj.avgSolution(i))*(1/(sum(s.totalCost./obj.avgSolution)));      
                end %for
            end %if
        end %function
        %% Reactive GRASP: Choose Random Alfa
        function alfa = chooseAlfa(obj)  
            alfa = randsample(obj.alfaArray,1,true,obj.alfaProb);
        end %function
    end %methods   
end %classdef