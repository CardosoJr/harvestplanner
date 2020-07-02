%% HarvestUnit Class
% Classe que descreve o talhão

% TODO: 
% - atualizar metodo showData()
classdef HarvestUnit
    %%
    properties
        woodVolume;     % Wood volume
        xCoord;         % Center x coordinates
        yCoord;         % Center y coordinates
        id;             % ID
        number;         % Number
        cutDiff;        % CutDiff. size[number of months to schedule]
        efficiency;     % Efficiency to be cut (m^3/ha)
        velocity;       % Velocity to be cut (m^3/h)
        RCL;            % Restricted Candidate List
    end %properties
    
    %%
    methods 
        
        %% Constructor
        % INPTUS
        %   data: vector which contains vital information to the construction of the object
        % OUTPUTS
        %   obj: HarvestUnit constructed object      
        function obj = HarvestUnit(data)
            % if the user choose to input the data argument
            if nargin > 0
                numMonths = length(data) - 7;
                obj.cutDiff = zeros(numMonths,1);

                obj.id = data(1);
                obj.number = data(2);
                obj.woodVolume = data(3);
                obj.xCoord = data(4);
                obj.yCoord = data(5);
                obj.efficiency = data(6);
                obj.velocity = data(7);
                obj.cutDiff(:) = data(8:end);
            % it the user didn't input the data argument    
            else  
                obj.id = 0;
                obj.woodVolume = 0;
                obj.xCoord = 0;
                obj.yCoord = 0;
                obj.cutDiff = 0;
                obj.efficiency = 0;
                obj.velocity = 0;
                
            end %if
        end
        %% Show Data Method
        % Print all Harvest Unit data in console 
        function showData(obj)
            fprintf('Printing data from HarvestUnit %d and Number: %d\n', obj.id, obj.number);
            fprintf('       Wood Volume: %d m^3', obj.woodVolume);
            fprintf('       X Coordinates: %d; Y Coordinates: %d\n', obj.xCoord, obj.yCoord);
            fprintf('       Cut Diff: ');
            numMonths = length(obj.cutDiff);
            for i=1:numMonths
                fprintf('%d |', obj.cutDiff(i));
            end
            fprintf('\n');

        end %function

        %% Calculate Cost Method
        % Calculate the total cost to harvest this specific harvest
        % unit
        function cost = calcCost(obj)
            cost = 0;
        end %function

        %% Calculate Transition Cost to another Harvest Unit Method
        % calculate the transition cost to move the machinary from this
        % harvest unit to another
        % Inputs:
        %   obj    -> "This" pointer to use this specific object
        %   obj_hu -> HarvestUnit object of the next unit to harvest on the
        %   farm schedule
        %   month  -> Month id of this particular farm schedule
        % Outputs: 
        %   cost   -> The cost of transition between harvest units
        function cost = calcTransitionCost(obj, obj_hu, coordCarb, month)
            cost = ((obj.xCoord - obj_hu.xCoord)^2 + (obj.yCoord - obj_hu.yCoord)^2)^0.5;
            cost = cost + obj_hu.cutDiff(month)*((coordCarb(1) - obj_hu.xCoord)^2 + (coordCarb(2) - obj_hu.yCoord)^2)^0.5;
        end %function
        
        %% Calculate Time To Harvest
        % Inputs: 
        %   obj    -> "This" pointer to use this specific object
        % Outputs:
        %   t -> The total time to harvest this unit
        function t = calcHarvestTime(obj)
            t = obj.woodVolume / obj.velocity;         
        end % function
        
        %% Plot HarvestUnit Coordinates
        % inputs: 
        %   h -> Figure Handle
        function plotCoordinates(obj, color, marker,size)
            p = plot(obj.xCoord,obj.yCoord);
            set(p,{'marker'},marker,'color',color, 'MarkerSize', size);
        end %function
        %% Plot Path Between HarvetUnit
        function plotPath(obj, obj_hu, color)
         %   pause
          plot([obj.xCoord obj_hu.xCoord],[obj.yCoord obj_hu.yCoord],'k--');
        end %function
        %% Write Data To File Method     
        function writeDataToFile(obj,fileID)
            fprintf(fileID,'Printing data from HarvestUnit %d and Number: %d\n', obj.id, obj.number);
            fprintf(fileID,'       Wood Volume: %d m^3', obj.woodVolume);
            fprintf(fileID,'       X Coordinates: %d; Y Coordinates: %d\n', obj.xCoord, obj.yCoord);
            fprintf(fileID,'       Cut Diff: ');
            numMonths = length(obj.cutDiff);
            for i=1:numMonths
                fprintf(fileID,'%d |', obj.cutDiff(i));
            end
            fprintf(fileID,'\n');
        end %function
    end %methods  
end %classdef