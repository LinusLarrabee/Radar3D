function Parameter()
global C Lambda PulseNum BandWidth TimeWidth PRT PRF Fs NoisePower Fc
%% Radar Parameter
C=3.0e8;  %光速(m/s)
Fc=200e9;  %雷达射频 1.57GHz
Lambda=C/Fc;    %雷达工作波长
PulseNum=24;   %回波脉冲数
BandWidth=4.0e6;  %发射信号带宽 带宽B=1/τ，τ是脉冲宽度
TimeWidth=4.0e-8; %发射信号时宽
PRT=8e-7;   % 雷达发射脉冲重复周期(s),240us对应1/2*240*300=36000米最大无模糊距离
PRF=1/PRT;
Fs=2.5e10;  %采样频率
NoisePower=5;%(dB);%噪声功率（目标为0dB）