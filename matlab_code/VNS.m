%% Variable Neighborhood Search Algorithm
% 
classdef VNS
    %% Properties
    properties
        neighborhood;
    end %properties
    %% Methods
    methods
        %% Constructor
        function obj = VNS()
            obj.neighborhood = Neighborhood();
            obj.neighborhood.createOperationVector();
        end
        %% Shake Method 
        function [s] = shake(obj, x_hat,max, min)
            fprintf('**** Calling Shake Method\n')
           numOperations = randi([min max]);
           w = 0;
           for i = 1:numOperations
                while (w ~= 1)
                    op = randi([1 obj.neighborhood.numOperation-3]);
                    h = obj.neighborhood.getHandle(op);
                    [x_hat,w] = h(x_hat);
                end %while
           end %for
           s = x_hat;
        end %function
        %% VNS Method
        function s = VND(obj, x_hat, k_max)
            fprintf('**** Calling VND Method\n')
            k = 1;
            while(k < k_max)
                op = randi([1 obj.neighborhood.numOperation - 3]);
                h = obj.neighborhood.getHandle(op);
                w = 0;
                while (w ~= 1)
                    [x,w] = h(x_hat);
                end %while
                x.calcCost();
                x_hat.calcCost();
                [x_hat, k] = obj.NeighborhoodChange(x,x_hat,k);
            end %while
            s = x_hat;
        end %function
        %% NeighborhoodChange MethodS
        function [s,k] = NeighborhoodChange(~,x, x_hat,k)
            if(x.totalCost < x_hat.totalCost)
                s = x;
                k = 1;
            else
                k = k + 1;
                s  = x_hat;
            end %if
        end %function
        %% General VNS Method
        function s = GeneralVNS(obj,x_hat, l_max,k_max, max,min)
            k = 1;
            while (k < k_max)
                x = obj.shake(x_hat, max,min);
                x = obj.VND(x,l_max);
                x.calcCost();
                x_hat.calcCost();
                [x_hat,k] = obj.NeighborhoodChange(x,x_hat,k);
            end %while
        s = x_hat;
        end %function
    end %methods
end
