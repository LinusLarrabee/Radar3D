function Environment()
% 本脚本文件需读取CarSet.mat，生成Environment.mat存放Velocity 和 Envir。
% 尚且存在的问题：暂时只允许最矮的车携带雷达，预计可使用Uni分开
%% 参数配置（应整合到Parameter文件中）
load('CarSet.mat');
Parameter();
global mete Lambda MaxDistance deltaD AllLine ...
    CarDis CarAngle Linewidth Height deltaAngle
%% 边界环境处理
% 只考虑环境尽头
Point = Compen(1,:);
[me,~] = size(Compen);
% 10个值分别表示SigPower, Line, Velocity, Place, (LWH)...
Selected = Compen((2:me),:);
Uni = unique(Selected(:,1)); %使用Uni能区分高矮车的不同能量。

Target = Selected(Selected(:,5)<=CarDis ...
    & Selected(:,6)<=CarAngle(2) & Selected(:,6)>=CarAngle(1),:);
[mt,~] = size(Target);

allAngle = round((CarAngle(2)-CarAngle(1))/deltaAngle);
Envir = zeros(allAngle,fix(MaxDistance/deltaD));
Velocity = zeros(allAngle,fix(MaxDistance/deltaD));
% 第一个元素为所有的角度分量，第二个元素为距离上的分辨率

%% 环境RCS计算
%按距离及角度分辨率在环境中划分单元，所有小单元看作是一个矩形，长为deltaD，宽为两弧线的中值
%最后按照实际面积进行修正
[~,n] = size(Envir);
queue = deltaD*(1/2:1:n-1/2);
distance =sqrt(queue.^2+Height^2);
theta = pi/2 - atan(Height./queue);

w = deltaD; %距离向长度
d = queue*deltaAngle; %方位向长度
fai = 0; %每次雷达相当于theta不定但是fai恒为0
Spa = 1; %面积矫正参数，经计算得出模拟举行和圆环的一部分面积相等。
k = 2*pi/Lambda;
deltaEnvir = pi^(-1)*(k.*d.*w).^2 .*sinc(k.*d.*sin(theta)*cos(fai)).^2 ...
    .*sinc(k*w*sin(theta)*sin(fai)).^2.*0.5.* ...
    ((1-sin(theta)*cos(fai).^2).^(0.5)+(1-sin(theta)*sin(fai).^2).^(0.5));
%% 环境反射区域计算
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

% 假设场景边界处的参数是一个材质不同的反射体（即更容易被感知到来作为最远距离）
for i = 1:length(Radius)
    Envir(i,Radius(i)) = mete(1);
    Envir(i,(1:Radius(i)-1)) =0;% deltaEnvir(1:Radius(i)-1);%0;不考虑RCS的话%
    Velocity(i,(1:Radius(i)-1)) =0;% -Point(3)*cos((i-300)*deltaAngle);%不考虑环境变量的话%
end

%% 车辆环境定位及RCS计算
% 竖向row，横向cal，则分别为三列的数组，row前两个元素代表y值，第三个元素代表x值
% cal则前两个元素代表x，第三个元素代表y

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
Line = fix((Line-CarAngle(1))/deltaAngle);


for i = 1:mt
    if Line(i,1)>Line(i,2)
        CEALR = Line(i,1):-1:Line(i,2);%Car_Environment Angle List Right
        CEALD = Line(i,1):1:Line(i,3);%Car_Environment Angle List Down
    else
        CEALR = Line(i,1):1:Line(i,2);
        CEALD = Line(i,1):-1:Line(i,3);%Car_Environment Angle List Down
    end
    for j= 2: length(CEALR)
        PorDis = CX(i)./cos(Sector(CEALR(j))); % 小面积斜距离
        theta = acos(CX(i)./sqrt(PorDis^2+0.25*Height^2));
        fai = atan(0.5*Height./sqrt(PorDis^2-CX(i)^2));
        w = Height; %距离向长度
        d = queue(fix(PorDis))*deltaAngle; %方位向长度
        Spa = 1; %面积矫正参数，这个真的懒得算了。
        k = 2*pi/Lambda;
        deltaCar1 = pi^(-1)*(k.*d.*w).^2 ...
            .*sinc(k.*d.*sin(theta)*cos(fai)).^2 ...
            .*sinc(k*w*sin(theta)*sin(fai)).^2.*0.5.* ...
            ((1-sin(theta)*cos(fai).^2).^(0.5) ...
            +(1-sin(theta)*sin(fai).^2).^(0.5));
        Envir(CEALR(j),fix(PorDis))=deltaCar1*mete(2);
        Envir(CEALR(j),(fix(PorDis)+1:n)) = 0;
        Velocity(CEALR(j),fix(PorDis))=Target(i,3) ...
            *cos((CEALR(j)-300)*deltaAngle);
        Velocity(CEALR(j),(fix(PorDis)+1:n)) = 0;
    end
    if CY(i) ~= 0 
        for j= 2: length(CEALD)
            PorDis = abs(CY(i)./sin(Sector(CEALD(j)))); % 小面积斜距离
            theta = acos(abs(CY(i))./sqrt(PorDis^2+0.25*Height^2));
            fai = atan(0.5*Height./sqrt(PorDis^2-CY(i)^2));
            w = Height; %距离向长度
            d = queue(fix(PorDis))*deltaAngle; %方位向长度
            Spa = 1; %面积矫正参数，这个真的懒得算了。
            k = 2*pi/Lambda;
            deltaCar = pi^(-1)*(k.*d.*w).^2 ...
                .*sinc(k.*d.*sin(theta)*cos(fai)).^2 ...
                .*sinc(k*w*sin(theta)*sin(fai)).^2.*0.5.* ...
                ((1-sin(theta)*cos(fai).^2).^(0.5) ...
                +(1-sin(theta)*sin(fai).^2).^(0.5));
            Envir(CEALD(j),fix(PorDis))=deltaCar*mete(2);
            Envir(CEALD(j),(fix(PorDis)+1:n))=0;
            Velocity(CEALD(j),fix(PorDis))=Target(i,3)...
                *cos((CEALD(j)-300)*deltaAngle);
            Velocity(CEALD(j),(fix(PorDis)+1:n))=0;
        end
        midval = fix(abs(CY(i)./sin(Sector(CEALD(1)))));
    else
        midval = CX(i);
    end
    Envir(CEALR(1),fix(midval/2+0.5*CX(i)./cos(Sector(CEALR(1))))) ...
        = deltaCar1*mete(2);
    Envir(CEALR(1),fix(midval/2+0.5*CX(i)/cos(Sector(CEALR(1)))+1:n))=0; 
    Velocity(CEALR(1),fix(midval/2+0.5*CX(i)./cos(Sector(CEALR(1))))) ...
        =Target(i,3)*cos((CEALR(1)-300)*deltaAngle);
    Velocity(CEALR(1),(fix(midval/2+0.5* ...
        CX(i)./cos(Sector(CEALR(1))))+1:n))=0;
    clearvars CEALR;
    clearvars CEALD;
end

% 所有车型均考虑为长方体，仅雷达车为流体形，雷达位于车顶端中部
% Point 表示参考车的车道，车速和相对水平轴的位置。
% 暂时默认所有姿态下均为斜矩形板，车都在中线上行驶，车的反射参数一致
% 暂时车在三维上仅能看到最近的一台车
%% 文件保存及输出
figure(11);
surf(Envir);
save Environment.mat Envir Velocity;
end