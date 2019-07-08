% 本脚本文件需读取CarSet.mat，调用BackPower函数，生成Environment.mat
%% 参数配置（未体系化）
load('CarSet.mat');
delta = BackPower();
% 第一个是边界的反射参数，第二个是车的反射参数，默认车的反射参数相同
EnviLength = 200;
EnviWidth = 30;
MaxDistance = sqrt(EnviLength^2+EnviWidth^2);
deltaD = 1;%距离分辨率
AllLine = 4;
CarDis = 100;
CarAngle = [-pi/3,pi/3];
Linewidth = 3.5; 

%% 边界环境处理（未考虑道路中线)
Point = Compen(1,:);
[me,ne] = size(Compen);
% 7个值分别表示SigPower, Line, Velocity, Place, (LWH)
Selected = Compen((2:me),:);
Uni = unique(Selected(:,1));

Target = Selected(Selected(:,5)<=CarDis ...
    & Selected(:,6)<=CarAngle(2) & Selected(:,6)>=CarAngle(1),:);
[mt,nt] = size(Target);
deltaAngle = 0.2*pi/180;%度
allAngle = round((CarAngle(2)-CarAngle(1))/deltaAngle);
Envir = zeros(allAngle,fix(MaxDistance/deltaD));
% 第一个元素为所有的角度分量，第二个元素为距离上的分辨率

%以下为扫描的扇形区域的delta值计算，近似为三角形。
Sector = deltaAngle+CarAngle(1):deltaAngle:CarAngle(2);
if Point(5)>0
    Radius= ( AllLine - Point(5) + 1/2 ) *Linewidth ./sin(-Sector);
    MinorR = ( AllLine  + Point(5) -1/2 ) *Linewidth ./sin(Sector);
else 
    Radius= ( AllLine - abs(Point(5)-1) + 1/2 ) *Linewidth ./sin(-Sector);
    MinorR = ( AllLine  + abs(Point(5)-1) -1/2 ) *Linewidth ./sin(Sector);
end
Radius= round(min(max(Radius,MinorR),MaxDistance));
for i = 1:length(Radius)
    Envir(i,Radius(i)) = delta(1);
end

%% 车辆环境处理（未区分不同车辆）
% 竖向row，横向cal，则分别为三列的数组，row前两个元素代表y值，第三个元素代表x值
% cal则前两个元素代表x，第三个元素代表y
%切入射角为

CX = Target(:,4)-Target(:,8)*0.5;
CY = -Linewidth * Target(:,2);

%C 代表最近的点
X = CX+Target(:,8);
Y = CY+sign(CY).*Target(:,9);

% global Point

%竖线
LineC = atan(CY./CX);
LineDown = atan(Y./CX);
LineRight = atan(CY./X);
Line = horzcat(LineC,LineDown,LineRight);
Line = fix(5*(Line+pi/3)*180/pi);


for i = 1:mt
    if Line(i,1)>Line(i,2)
        CEALR = Line(i,1):-1:Line(i,2);%Car_Environment Angle List Right
        CEALD = Line(i,1):1:Line(i,3);%Car_Environment Angle List Down
    else
        CEALR = Line(i,1):1:Line(i,2);
        CEALD = Line(i,1):-1:Line(i,3);%Car_Environment Angle List Down
    end
    for j= 2: length(CEALR)
        Envir(CEALR(j),fix(CX(i)./cos(Sector(CEALR(j)))))=delta(2);
    end
    if CY(i) ~= 0 
        for j= 2: length(CEALD)
            Envir(CEALD(j),fix(abs(CY(i)./sin(Sector(CEALD(j))))))=delta(2);
        end
        midval = fix(abs(CY(i)./sin(Sector(CEALD(1)))));
    else 
        Envir(CEALD(2:length(CEALD)),CX(i))=delta(2);
        midval = CX(i);
    end
    Envir(CEALR(1),fix(midval/2+0.5*CX(i)./cos(Sector(CEALR(1))))) = delta(2);
    clearvars CEALR;
    clearvars CEALD;
end
%Envir(Envir(:,1)>=Line(:,1) & Envir(:,1)<=Line(:,2),
%round(X./cos(Sector)) = delta(2);
% 所有车型均考虑为长方体，仅雷达车为流体形，雷达位于车顶端中部
% Point 表示参考车的车道，车速和相对水平轴的位置。
% 暂时默认所有姿态下均为斜矩形板，车都在中线上行驶，车的反射参数一致
% 暂时车在三维上仅能看到最近的一台车
%% 文件保存及输出
figure(10)
surf(Envir);
save Environment.mat Envir;