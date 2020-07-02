%% Load Data 
% Ler dados do Banco de Dados

function [MS, SchedData] = loadData(fileName, sheetName) 

    %% Inputs
    % fileName: name of the file to read the input data of the optimization
    % software
    % sheetName: name of the sheet of xls file which contains all data
    %% Outputs
    % RegList
    % MonList
    %% Algorithm
    data = xlsread(fileName);
    demand = xlsread(fileName,'Demanda');
    demand = demand(2:end, 2:end);
    Max = max(max([data(:,4) data(:,5)]));
    Min = min(min([data(:,4) data(:,5)]));
    data(:,4) = 100*((data(:,4)-Min)/(Max-Min));
    data(:,5) = 100*((data(:,5)-Min)/(Max-Min));
    data(:,13) = 100*((data(:,13)-Min)/(Max-Min));
    data(:,14) = 100*((data(:,14)-Min)/(Max-Min));
    
    if data == -1 
        warning('The file could not be opened or it does not exist')
    end
    
    numFarms = data(1,10);
    numHU = zeros(numFarms,1);
    numHarvestUnits = 0;
    
    for i = 1:numFarms 
        numHarvestUnits = numHarvestUnits + data(i,11);
        numHU(i) = data(i,11);
    end
   
    numMonths = data(1,12);
    demand = demand(:,1:numMonths+6);
    demand(:,end) = 0;
    demand(:,end-1) = 0;
        
    carbCoord = data(1:4,13:14);
    
    %Inicializações
    huData = zeros(7+numMonths,1);
    sdData = zeros(numMonths, 3);
    Farms = Farm.empty(numFarms,0);
    
    % Inicialmente supoe-se que os talhoes estao agrupados
    % por fazendas
    m = 1;
    for j =1:numFarms
        
        hu = HarvestUnit.empty(data(j,11),0);
        index = 1;
        for i = m : m + data(j,11) - 1
              
              huData(1:7) = data(i,1:7);
                            
              huData(8:end) = data(i,15:numMonths+14);
              hu(index) = HarvestUnit(huData);
              index = index + 1;
        end
        m = m + data(j,11);
        if j == 1
            name = 'Pin';
        elseif j == 2
            name = 'Vgr';
        elseif j == 3
            name = 'Ext';
        elseif j ==4 
            name = 'NII'; 
        end
        Farms(j) = Farm(name,j, hu,carbCoord(j,:));
    end
    
    sdData(:,:) = data(1:numMonths,28:30);
    macVel = data(1,end);        
    SchedData = SchedulerData(sdData,numHU, macVel,demand(1:end,2:end), demand(:,1)); 
    
    MS = MonthSchedule();
    MS.listFarm = Farms;
    MS.numFarm = length(Farms);    
end
