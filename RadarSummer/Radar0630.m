
global C PulseNum BandWidth TimeWidth PRT Fs NoisePower Fc
C
PulseNum

load('cars.mat')

Target = newData;
SigPower = Target(:,1)'; %目标功率,无量纲 [5,1,100,10.250000000000000]
TargetDistance = Target(:,2)'; %目标距离,单位m  距离参数为[3000 8025 15800 8025]
TargetVelocity = Target(:,3)'; %目标径向速度 单位m/s  速度参数为[50 0 10000 100]
TargetAngle = round(Target(:,4))'; %目标角度 单位度 角度参数 [30 60 90 120]
TargetWaveNum = length(SigPower)