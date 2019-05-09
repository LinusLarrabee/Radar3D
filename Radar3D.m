%% Initial and introduction of the process

% Radar3D.m
% Created by Group 6 on 5/8/19.
close all; clear ; clc;

%% Radar Parameter
C=3.0e8;  %光速(m/s)
RF=3.140e9/2;  %雷达射频 1.57GHz
Lambda=C/RF;    %雷达工作波长
PulseNum=16;   %回波脉冲数
BandWidth=4.0e6;  %发射信号带宽 带宽B=1/τ，τ是脉冲宽度
TimeWidth=10.0e-6; %发射信号时宽
PRT=240e-6;   % 雷达发射脉冲重复周期(s),240us对应1/2*240*300=36000米最大无模糊距离
PRF=1/PRT;
Fs=1e8;  %采样频率
NoisePower=-12;%(dB);%噪声功率（目标为0dB）
Fc = 30e6;

%% 目标参数
Target = load ('test.txt');
SigPower = Target(:,1)'; %目标功率,无量纲 [5,1,100,10.250000000000000]
TargetDistance = Target(:,2)'; %目标距离,单位m  距离参数为[3000 8025 15800 8025]
TargetVelocity = Target(:,3)'; %目标径向速度 单位m/s  速度参数为[50 0 10000 100]
TargetAngle = round(Target(:,4))'; %目标角度 单位度 角度参数 [30 60 90 120]
TargetWaveNum = length(SigPower);
DelayNum = fix(Fs*2*TargetDistance/C);% 把目标距离换算成采样点（距离门） fix函数向0靠拢取整
TargetFd = 2*TargetVelocity/Lambda;%计算目标多卜勒频移2v/λ


%% 信号参数-产生线性调频信号

WaveNum=fix(Fs*TimeWidth);%回波的采样点数=脉压系数长度=暂态点数目+1
if mod(WaveNum,2)~=0
    WaveNum=WaveNum+1;
end   %把WaveNum变为偶数

Chirp = zeros(1,WaveNum);
RealChirp = zeros(1,WaveNum);
for i=-fix(WaveNum/2):fix(WaveNum/2)-1
    Ft = Fc*i/Fs+(1/2)*(BandWidth/TimeWidth)*(i/Fs)^2; %线性调频波指数幂
    Chirp(i+fix(WaveNum/2)+1)=exp(1i*2*pi*Ft); %exp(j*fi)*，产生复数矩阵Chirp
    RealChirp(i+fix(WaveNum/2)+1) =cos(2*pi*Ft);
end

coeff=conj(fliplr(Chirp)); %把Chirp矩阵翻转取共轭，产生脉压系数
% 简单解释一下就是这里取的是 h(t) = x(^*)(t_0-t) 中 t_0 = 0 的情景

%% 相关参数计算
totalAngle = 120;
startAngle = 30;
deltaAngle = 3;

AngleNum = totalAngle/deltaAngle;
Sector = startAngle + deltaAngle*(1:AngleNum);

SampleNum=fix(Fs*PRT);%计算一个脉冲周期的采样点数；
AngleCirNum = SampleNum*AngleNum; 
%先完成一个角度的一整个线性调频，再依次完成剩余39个角度，最后循环16次脉冲；
TotalNum=SampleNum*PulseNum*AngleNum;%总的采样点数；
BlindNum=fix(Fs*TimeWidth);%计算一个脉冲周期的盲区-遮挡样点数；

%% 产生目标回波串
% 产生前3个目标的回波串
SignalAll=zeros(1,TotalNum);%所有脉冲的信号,先填0
for ang = 1:AngleNum
    detect = find (TargetAngle == Sector(ang));
    if ~isempty(detect)
        for k = detect % 依次产生各个目标
            fi=2*pi/10 * fix(10*rand);
            SignalTemp=zeros(1,SampleNum);% 一个PRT
            SignalTemp(DelayNum(k)+1:DelayNum(k)+WaveNum)=...
                sqrt(SigPower(k))*exp(1i*fi)*Chirp;
            %一个脉冲的1个目标（未加多普勒速度）(DelayNum(k)+1):(DelayNum(k)+WaveNum)
            Signal = zeros(1,TotalNum);
            SignalAng = zeros(1,AngleCirNum);
            SignalAng((ang-1)*SampleNum+1:ang*SampleNum)=SignalTemp;
            for i=1:PulseNum % 16个回波脉冲
                Signal((i-1)*AngleCirNum+1:i*AngleCirNum)=SignalAng;
                %每个目标把16个SignalTemp排在一起
            end
            FreqMove=exp(1i*2*pi*TargetFd(k)*(0:TotalNum-1)/Fs);
            %目标的多普勒速度*时间=目标的多普勒相移
            Signal=Signal.*FreqMove;%加上多普勒速度后的16个脉冲1个目标
            SignalAll=SignalAll+Signal;%加上多普勒速度后的16个脉冲4个目标
        end
    end
end

figure(2);
subplot(2,1,1);plot(real(SignalAll),'r-');title('目标信号的实部');...
grid on;zoom on;
subplot(2,1,2);plot(imag(SignalAll));title('目标信号的虚部');grid on;zoom on;


%% 总的回波信号
% 产生系统噪声信号
SystemNoise=normrnd(0,10^(NoisePower/10),1,TotalNum)...
+1i*normrnd(0,10^(NoisePower/10),1,TotalNum);
%均值为0，标准差为10^(NoisePower/10)的噪声

%闭锁期无回波
EchoAll=SignalAll+SystemNoise;% +SeaClutter+TerraClutter，加噪声之后的回波
for i=1:PulseNum*AngleNum   %在接收机闭锁期,接收的回波为0
    EchoAll((i-1)*SampleNum+1:(i-1)*SampleNum+WaveNum)=0; %发射时接收为0
end
figure(3);%加噪声之后的总回波信号
subplot(2,1,1);plot(real(EchoAll),'r-');title('总回波信号的实部,闭锁期为0');
subplot(2,1,2);plot(imag(EchoAll));title('总回波信号的虚部,闭锁期为0');

%% 回波信号整形
EchoRoute = reshape(EchoAll, [SampleNum,AngleNum,PulseNum]);
for i = 1 : 10
    Echo = reshape(EchoRoute(:,i,:),1,[]);
    figure(i+3)
    subplot(3,1,1)
    plot(real(Echo),'r-');title('总回波信号的实部,闭锁期为0');
    %% 时域脉压
    pc_time0=conv(Echo,coeff);%pc_time0为Echo和coeff的卷积
    pc_time1=pc_time0(WaveNum:length(Echo)+WaveNum-1);%去掉暂态点 WaveNum-1个
    %figure(4);%时域脉压结果的幅度
    subplot(3,1,2);plot(abs(pc_time0),'r-');title('时域脉压结果的幅度,有暂态点');
    %pc_time0的模的曲线
    subplot(3,1,3);plot(abs(pc_time1));title('时域脉压结果的幅度,无暂态点');
    %pc_time1的模的曲线
    
end
%% 时域脉压
pc_time0=conv(Echo,coeff);%pc_time0为Echo和coeff的卷积
pc_time1=pc_time0(WaveNum:length(Echo)+WaveNum-1);%去掉暂态点 WaveNum-1个
figure(4);%时域脉压结果的幅度
subplot(2,1,1);plot(abs(pc_time0),'r-');title('时域脉压结果的幅度,有暂态点');
%pc_time0的模的曲线
subplot(2,1,2);plot(abs(pc_time1));title('时域脉压结果的幅度,无暂态点');
%pc_time1的模的曲线
%% 频域脉压
Echo_fft=fft(Echo,524288);
%理应进行length(Echo)+WaveNum-1点FFT,但为了提高运算速度,进行了8192点的FFT
coeff_fft=fft(coeff,524288);
pc_fft=Echo_fft.*coeff_fft;
pc_freq0=ifft(pc_fft);
figure(5);
subplot(2,1,1);plot(abs(pc_freq0(1:length(Echo)+WaveNum-1)));
title('频域脉压结果的幅度,有前暂态点');
subplot(2,1,2);
plot(abs(pc_time0(1:length(Echo)+WaveNum-1)-...
         pc_freq0(1:length(Echo)+WaveNum-1)),'r');
title('时域和频域脉压的差别');
pc_freq1=pc_freq0(WaveNum:length(Echo)+WaveNum-1);
%去掉暂态点 WaveNum-1个,后填充点若干(8192-WaveNum+1-length(Echo))
% ================按照脉冲号、距离门号重排数据=================================%
for i=1:PulseNumber
pc(i,1:SampleNumber)=pc_freq1((i-1)*SampleNumber+1:i*SampleNumber);
%每个PRT为一行，每行480个采样点的数据
end
figure(6);
plot(abs(pc(1,:)));title('频域脉压结果的幅度,没有暂态点');

% ================MTI（动目标显示）,对消静止目标和低速目标---可抑制杂波%
for i=1:PulseNumber-1  %滑动对消，少了一个脉冲
mti(i,:)=pc(i+1,:)-pc(i,:);
end
figure(7);
mesh(abs(mti));title('MTI  result');

% ================MTD（动目标检测）,区分不同速度的目标，有测速作用==%
mtd=zeros(PulseNumber,SampleNumber);
for i=1:SampleNumber
buff(1:PulseNumber)=pc(1:PulseNumber,i);
buff_fft=fft(buff);
mtd(1:PulseNumber,i)=buff_fft(1:PulseNumber);
end
figure(8);mesh(abs(mtd));title('MTD  result');

%% 虚实矩阵转换

coeff_fft_c=zeros(1,2*524288);
for i=1:8192
coeff_fft_c(2*i-1)=real(coeff_fft(i));
coeff_fft_c(2*i)=imag(coeff_fft(i));
end
echo_c=zeros(1,2*length(Echo));
for i=1:length(Echo)
echo_c(2*i-1)=real(Echo(i));
echo_c(2*i)=imag(Echo(i));
end
%% 以下是为DSP程序提供回波数据、脉压系数
[fo,message] = fopen('/Users/sunhao/Desktop/Project0.1/coeff_fft_c.dat'...
                     ,'wt');%频域脉压系数
if fo < 0
error('Failed to open coeff_fft_c.dat becauese: %s',message);
end
for i=1:2*8192
fprintf(fo,'%f,\r\n',coeff_fft_c(i));
end
fclose(fo);

[fo,message]=fopen('/Users/sunhao/Desktop/Project0.1/echo_c.dat'...
                   ,'wt');%16次回波的
if fo < 0
error('Failed to open echo_c.dat becauese: %s',message);
end
for i=1:2*length(Echo)
fprintf(fo,'%f,\r\n',echo_c(i));
end
fclose(fo);
