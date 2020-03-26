clear;clc; close all;
%% 环境设置
tar_dis = [3000,5000,7000]; % Target Distance
RCS = [1,1,1]; % Target Power
tar_vel = 3*[1,2,3]; %Target Velocity
tar_vel_phi = [0,0.5,1,1.5]; % Target Velocity Phase 暂时不考虑MIMO的速度位置更新
tar_angle = [0.5 0.7 0.9]; %Target Angle (rad)
tar_num = length(tar_dis); 
save target.mat tar_dis RCS tar_vel tar_angle tar_num tar_vel_phi

%% 天线设置
% 假设MIMO为线性雷达，RX，TX均位于x轴
% Tx有两根，位于Rx两端，在Tx之间等间距放置Rx
group = 2;
dt = [100,500]'; % Tx之间的间距
dr = [20,100]'; % Rx之间的间距
tx1 = [0,1000]'; % 每个雷达的最左侧tx的位置
rx1 = [20,100]'; %每个雷达的最左侧rx的位置
tx_num = [2,3];
rx_num = [4,6];
Tx = dt*(1:tx_num);
Rx = dr*(0:rx_num-1);
save MIMO.mat group dt dr tx1 rx1 tx_num rx_num
