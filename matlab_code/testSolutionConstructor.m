%% TEST SolutionConstructor

clear
close all
clc
load('data.mat');

schedData.max_construct = 0.80;
schedData.min_construct = 1.2;

schedData.epsilon = 0.8;
schedData.alfa = 0.5;

%[x_hat, garbage] = SolutionConstructor.randomizedSolution(schedData,MS);

[x_hat, garbage] = SolutionConstructor.greedyAndRandomizedSolution1(schedData, MS);
% x_hat.ID = 1;
%x_hat.showData();
% 
% for i = 1:x_hat.numMonthsToSchedule
%    for j = 1:x_hat.listMonthSchedule(i).numFarm
%         fprintf('%d\n',x_hat.listMonthSchedule(i).listFarm(j).numHarvestUnit);
%    end
% end
% 
% fprintf('\n\n');
% 
% for i = 1:x_hat.numMonthsToSchedule
%    fprintf('%d\n', x_hat.listMonthSchedule(i).numFarm);
% end
% 

