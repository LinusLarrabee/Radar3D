function [PointData] = CarSet(pointcar)
% 保存CarSet.mat含十列矩阵Compen，第二行起列为cars.mat/距离与角度/CarSeries。Compen第一列为pointcar，五，六值为道数和距离
% 导入的两个mat均改为相对pointcar的值
% 这里的disthorizon表示车道上离雷达的水平距离
% line表示车道，pointcar表示第几个车作为雷达位置。
load('CarSeries.mat');
load('cars.mat');
% 第一列车型，第二列车道数，第三列速度，第四列车头距离某一节点的绝对距离
Cardata = newData;%xlsread ('Cars.xlsx');
CarSeries = carSeries;
%存储参考值；
PointData = Cardata(pointcar,:);
PointData(1) = 0;
Cardata=Cardata-PointData;
Larger = zeros(length(Cardata),4);
Uni = unique(Cardata(:,1));
for i = 1:length(Uni)
    Larger(Cardata(:,1) ==Uni(i),:) = repmat(CarSeries(i,:),...
        length(find(Cardata(:,1) ==Uni(i))),1);
end
Line = Cardata(:,2);
Disthorizon = Cardata(:,4);
%正则代表同向，负则代表逆向
Linewidth = 3.5; 
Distance = ((Linewidth*Line).^2+Disthorizon.^2).^0.5;
Angle = -atan(Linewidth*Line./Disthorizon);
Angle(Disthorizon<0) = Angle(Disthorizon<0) + pi;
    %在二三象限的值需要加pi来补，最后角度范围是-pi/2～3/2*pi
Compen = horzcat(Cardata,Distance,Angle,Larger);
Compen([1;pointcar],:)=Compen([pointcar;1],:);
Compen(1,[5,3,6]) = PointData(1,[2,3,4]);
save CarSet.mat Compen;
end
