clear ;close all;
%% Radar Parameter
C = 3.0e8; %光速
Fc = 30e9; %毫米波30～300GHz频域(波长为1～10mm
Lambda = C/Fc;

PRT = 1.0e-4;
time_width=1.0e-5; %发射信号时宽
K = 1e12;
band_width=K*time_width;  %发射信号带宽 带宽B=1/τ，τ是脉冲宽度

pulse_num=16;   %回波脉冲数
Fs = 1.2 *band_width; %采样频率，建议取带宽1.2倍以上

%% Environment
load('target.mat');
% 待增加附加功能（按分辨率取整）

delay_num = fix(Fs*2*tar_dis/C);% 把目标距离换算成采样点（距离门） fix函数向0靠拢取整
phi = 0;%暂时不考虑target的更新，雷达也仅仅在一个方向
tar_fd = 2*tar_vel/Lambda;%.*abs(cos(tar_vel_phi-phi))%计算目标多卜勒频移2v/λ
tar_power = RCS./tar_dis.^4;
retar_power = ones(1,30);%tar_power/mean(tar_power);
noise_power = 0; %噪声功率，单位dB

%% 1-D 连续复信号
wave_num_half=floor(Fs*time_width/2); % WaveNumHalf：回波采样点数的一半
wave_num = wave_num_half*2; 
%回波的采样点数=脉压系数长度=暂态点数目+1，且收发分离不用考虑目标较近的盲区

A0 = -wave_num_half:wave_num_half-1; % 循环符
Ft = Fc*A0/Fs+(1/2)*K*(A0/Fs).^2; %线性调频波指数幂
Chirp=exp(1i*2*pi*Ft); %exp(j*fi)*，产生复数矩阵Chirp

fig1 = figure(1);
subplot(2,1,1);plot(real(Chirp));title('发波信号实部');
subplot(2,1,2);plot(imag(Chirp));title('发波信号虚部');
saveas(fig1, '发射信号.png');
%% MIMO
% load('MIMO.mat');
% alpha = [1,1]
% theta = pi/15;
% for i = 1:group
%     s = 
%     [at{i},ar{i}] = MIMO(dt(i),dr(i),tx_num(i),rx_num(i),Lambda,theta)
%     touch = alpha(1)*at{i}'*s;
% end

%% 2-D 回波生成

sample_num = floor(Fs*PRT);%计算一个脉冲周期的采样点数；
total_num = sample_num*pulse_num;

A1 = 1:wave_num;
A2 = 1:total_num;

deltaT = PRT;
signal_temp = zeros(tar_num,sample_num);
signal= zeros(tar_num,total_num);
for k = 1:tar_num
    phase=exp(1i*2*pi*fix(180*rand)/180); %随机相移
    signal_temp(k,delay_num(k)+A1) = sqrt(retar_power(k))*Chirp*phase;
    freq_move = exp(1i*2*pi*tar_fd(k)*A2/Fs); %多普勒频移
    signal(k,:) = repmat(signal_temp(k,:),1,pulse_num).*freq_move;
end

% 产生系统噪声信号
noise = normrnd(0,10^(noise_power/10),1,total_num)...
    +1i*normrnd(0,10^(noise_power/10),1,total_num);
%均值为0，标准差为10^(NoisePower/10)的噪声

signal_temp2 = sum(signal,1);
echo = signal_temp2+noise;

fig2 = figure(2);
subplot(2,1,1);plot(real(echo));title('回波信号实部');
subplot(2,1,2);plot(imag(echo));title('回波信号虚部');
saveas(fig2, '回波信号.png');
fig21 = figure(21);
subplot(2,1,1);plot(real(signal_temp2));title('回波信号实部(无噪声）');
subplot(2,1,2);plot(imag(signal_temp2));title('回波信号虚部（无噪声）');
saveas(fig21, '回波信号（无噪声）.png');
%% 脉冲压缩

% Generate filter coeff
coeff=conj(fliplr(Chirp)); %把Chirp矩阵翻转取共轭，产生脉压系数
% 这里取的是 h(t) = x(^*)(t_0-t) 中 t_0 = 0 的情景
%coeff=conj(fliplr(Chirp))*hamming(WaveNum);

% 时域脉冲压缩
conv_len = length(echo)+wave_num-1; %时域卷积长度
pc_time0 = conv(echo,coeff);
pc_time1 = pc_time0(wave_num:conv_len);%去掉暂态点 WaveNum-1个，结果为Echo长度

fig3 = figure(3);
subplot(2,1,1);plot(abs(pc_time0));title('时域脉压结果的幅度,有暂态点');
subplot(2,1,2);plot(abs(pc_time1));title('时域脉压结果的幅度,无暂态点');
saveas(fig3, '时域脉冲压缩.png');
% 频域脉冲压缩
echo_fft = fft(echo,conv_len);
coeff_fft = fft(coeff,conv_len);
pc_fft = echo_fft.*coeff_fft;
pc_freq0 = ifft(pc_fft);
pc_freq1 = pc_freq0(wave_num:conv_len);

fig4 = figure(4);
subplot(2,1,1);plot(abs(pc_freq1),'r-');title('频域脉压结果的幅度,无暂态点');
subplot(2,1,2);plot(abs(pc_time0-pc_freq0));title('时域和频域脉压的差别');
saveas(fig4, '时频域脉冲压缩对比.png');

pc = reshape(pc_freq1,sample_num,pulse_num)';% 奇怪的操作 reshape先列后行
fig5 = figure(5);
plot(abs(pc(1,:)));title('频域脉压结果的幅度,没有暂态点');
saveas(fig5, '频域单脉冲压缩.png');
%% MTI 及 MTD
% MTI 静目标及杂波对消
mti = pc(1:pulse_num-1,:) - pc(2:pulse_num,:); %滑动对消，少了一个脉冲

fig6 = figure(6);
mesh(abs(mti));title('MTI result');
saveas(fig6, '单延迟对消器.png');

% MTD（动目标检测）,区分不同速度的目标，有测速作用
mtd=zeros(pulse_num,sample_num);


for i=1:sample_num
    buff=pc(:,i);
    mtd(:,i)=fft(buff);
end
fig7=figure(7)
mesh(abs(mtd));title('MTD result');
saveas(fig7, 'MTD.png');
%  找到多普勒——距离地图中峰值最高的目标
abs_mtd = abs(mtd); %对正交分量的复数求模
max_target = max(max(abs_mtd));
[row,cell] = find( abs_mtd == max_target);
target_D = cell/Fs*C/2;
target_V = (fix(row - 1)/pulse_num)*Lambda/2/PRT;