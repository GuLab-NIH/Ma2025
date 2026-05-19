%%

clear; close all; clc

p1 = 'D:\4m_optoPaper\Histology';
cd(p1)

% load data
supp = '_ChR2';
load(['histData' supp '.mat'],'dataH')
fNames = fieldnames(dataH);
nF = length(fNames);

% generate plot array
plotData = cell(1,nF);
for ff = 1:nF
    plotData{ff} = dataH.(fNames{ff})(:);
end

% plot data
figure; hold on
barGroup(fNames,plotData,[],[],[],[],[],2)



