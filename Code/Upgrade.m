function Upgrade(deltaT)
load('cars.mat');
newData(:,4) = newData(:,4)+newData(:,3)*deltaT;
save cars.mat newData;
global PointNum
CarSet(PointNum);
Environment;
load('Environment.mat');

