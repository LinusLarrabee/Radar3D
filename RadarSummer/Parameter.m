function Parameter()
global C Lambda PulseNum BandWidth ...
    TimeWidth PRT PRF Fs NoisePower Fc WaveNum
%% Radar Parameter
C=3.0e8;  %光速(m/s)
Fc=100e6;  %雷达射频 1.57GHz
Lambda=C/Fc;    %雷达工作波长
PulseNum=24;   %回波脉冲数
BandWidth=2e8;  %发射信号带宽 带宽B=1/τ，τ是脉冲宽度
TimeWidth=4.0e-8; %发射信号时宽
PRT=8e-7;   % 雷达发射脉冲重复周期(s),
%240us对应1/2*240*300=36000米最大无模糊距离
PRF=1/PRT;
Fs=2.5e10;  %采样频率
NoisePower=5;%(dB);%噪声功率（目标为0dB）

WaveNum=fix(Fs*TimeWidth);%回波的采样点数=脉压系数长度=暂态点数目+1
if mod(WaveNum,2)~=0
    WaveNum=WaveNum+1;
end   %把WaveNum变为偶数

%满足的第一个条件是：Ft在Wavenum单次的间隔不大，或者说每个2pi周期内有一定的采样点
% 若想在函数图像上看到经典的Chirp曲线

carSeries = [5,4.5,2.0,1.2;3,8,2.5,3];
save CarSeries.mat carSeries;
end
