%% Main: 
% Schedule harvest units in certain region which contains a number of farms
% and each one has a number of harvest units
% This schedule must respect demand and time restrictions

%%
clc
clear
close all
readDataFromExcel = 0;

%% Read data and construct data structures
% Read data from data base
if (readDataFromExcel == 1) 
    fprintf('Reading Data from Excel...\n')
    [MS,schedData] = loadData('data.xlsx','');
else
    fprintf('Reading Data From Matlab Variables...\n');
    load('data.mat');
end
%% Define Inputs
fprintf('Defining Algorithm Inputs...\n');

schedData.max_k = 10;
schedData.max_vns_k = 10;
schedData.max_vnd_k = 10;
schedData.max_shake_operations = 20;
schedData.min_shake_operations = 10;
schedData.max_construct = 0.85;
schedData.min_construct = 1.5;
best_solution = Solution();
best_solution.totalCost = inf;

schedData.writeDataToFile(fileID);
fprintf('\nREADY TO OPTMIZE\n');

%% Call GRASP
[s,fobj] = GRASP.debuggGrasp(MS,schedData);
%% Pos processing
k = 1:schedData.max_k;
figure;
    plot(k,fobj,'-s');
    title('Cost Function over Iterations');
    xlabel('Number of Iterations');
    ylabel('Cost Function');
s.plotPath(); 

%% Escrever planilha Excel





