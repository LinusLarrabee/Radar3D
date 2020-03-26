clear ;close all;
%% Radar Parameter
C = 3.0e8; %光速
Fc = 30e9; %毫米波30～300GHz频域(波长为1～10mm
Lambda = C/Fc;

PRT = 1.0e-4;
time_width=4.0e-5; %发射信号时宽
band_width=4.0e6;  %发射信号带宽 带宽B=1/τ，τ是脉冲宽度
mu = band_width/time_width;
pulse_num=16;   %回波脉冲数

% MIMO 多频段
M = 1; % 发射波天线数
N = 16; % 16个接收阵元
K = 160;
deltaF = K/time_width;
B = band_width + (M-1)*deltaF;
Fs = 3 *B; %采样频率，建议取带宽1.2倍以上

%% Environment
tar_dis = [3000,4000,5000]; % Target Distance
RCS = [1,1,1]; % Target Power
tar_vel = [1,1,1]; %Target Velocity
tar_vel_phi = [0,0.5,1,1.5]; % Target Velocity Phase 暂时不考虑MIMO的速度位置更新
tar_angle = [0.5 0.7 0.824]; %Target Angle (rad)
tar_num = length(tar_dis); 
save target.mat tar_dis RCS tar_vel tar_angle tar_num tar_vel_phi
%load('target.mat');
% 待增加附加功能（按分辨率取整）

delay_num = fix(Fs*2*tar_dis/C);% 把目标距离换算成采样点（距离门） fix函数向0靠拢取整
phi = 0;%暂时不考虑target的更新，雷达也仅仅在一个方向
tar_fd = 2*tar_vel/Lambda;%.*abs(cos(tar_vel_phi-phi))%计算目标多卜勒频移2v/λ
tar_power = RCS./tar_dis.^8; %幅度能量简化为RCS与距离八次方反比，考虑密集式MIMO，收发距离相近
retar_power = tar_power/mean(tar_power);
noise_power = 0; %噪声功率，单位dB

%% 1-D 连续复信号
wave_num_half=floor(Fs*time_width/2); % WaveNumHalf：回波采样点数的一半
wave_num = wave_num_half*2; 
%回波的采样点数=脉压系数长度=暂态点数目+1，且收发分离不用考虑目标较近的盲区


A0 = -wave_num_half:wave_num_half-1; % 循环符
Ft = Fc*A0/Fs+(1/2)*mu*(A0/Fs).^2; %线性调频波指数幂
Chirp=exp(1i*2*pi*Ft); %exp(j*fi)*，产生复数矩阵Chirp
fig1 = figure('visible','off');
s = zeros(M,wave_num);
for i = 1:M
    Ft = (Fc+(i-1)*deltaF)*A0/Fs+(1/2)*mu*(A0/Fs).^2; %发射正交波形
    s(i,:) = exp(1i*2*pi*Ft);
    plot(A0,fftshift(abs(fft(s(i,:)))));
    hold on;
end
title('MIMO发射信号频谱');
saveas(fig1, 'MIMO发射信号频谱.png');

fig2 = figure('visible','off');
subplot(2,1,1);plot(real(Chirp));title('发波信号实部');
subplot(2,1,2);plot(imag(Chirp));title('发波信号虚部');
saveas(fig2, '发射信号.png');

A01 = -wave_num:wave_num-1; % 两倍长度循环符
beta = -1*deltaF +mu*A01/Fs; % -1代表第一个减第二个

C12 = sin(pi*beta.*(time_width-abs(A01)/Fs)).*(time_width-abs(A01)/Fs)...
    ./(pi*beta.*(time_width-abs(A01)/Fs))/time_width;
fig3 = figure('visible','off');
plot(abs(C12));
title3 = ['K =' num2str(K) '时相邻两频段互相关函数'];
filename3 = [title3 '.png'];
title(title3);
saveas(fig3, filename3);

%% MIMO

% load('MIMO.mat');
% alpha = [1,1]
% theta = pi/15;
% for i = 1:group
%     s = 
%     [at{i},ar{i}] = MIMO(dt(i),dr(i),tx_num(i),rx_num(i),Lambda,theta)
%     touch = alpha(1)*at{i}'*s;
% end

dt = N*Lambda/2;
dr = Lambda/2; % 接收阵元间距
scan = (0:0.01:1.59);

omega_send = 2*pi*dt*sin(tar_angle)/Lambda;
a_send = exp(1i*-(0:1:M-1)'*omega_send); %发射导向矢量 默认排布符合最简
sig_t = a_send' *s;%到达目标的信号
for i = 1:tar_num
    sig_t(i,:) = retar_power(i)*sig_t(i,:); 
end
omega_back = 2*pi*dr*sin(tar_angle)/Lambda;
a_back = exp(1i*-(0:1:N-1)'*omega_back); %接收导向矢量 默认排布符合最简

sig_r = a_back*sig_t;
Rxx = sig_r*sig_r'; % 直接这样求自相关吗？

%% Bartlett
Bartlett = zeros(1,length(scan)); % Bartlett的功率谱
for i = 1:length(scan)
    omega2 = 2*pi*dr*sin(scan(i))/Lambda;
    a_send1 = exp(1i*-(0:1:M-1)'*omega2);
    a_back2 = exp(1i*-(0:1:N-1)'*omega2);
    Bartlett(i) = a_back2'*Rxx*a_back2;
end
fig11 = figure('visible','off');
plot(scan,log(abs(Bartlett)));
saveas(fig11, 'Bartlett.png');
%% 信源估计
% Rxx 信息论方法
[V,D]= eig(Rxx);
sz_Rxx = size(Rxx);
[d,~] = sort(diag(D),'descend');
[~,P] = max(d(1:N-1)./d(2:N)); %求特征值梯度最大的值，并赋值给P1（估计的目标数）
% Test1 = Rxx*V(:,1)-D(1,1)*V(:,1);

% AIC与MDL
%AIC = 2*k-2*log(L);
AIC = zeros(1,N); %k的数目按照Rxx大小来取
MDL = zeros(1,N);
L = 10;
for k = 0:N-1
    AIC_temp1 = prod(d(k+1:N)).^(1/(N-k))/sum(d(k+1:N))/(N-k);
    AIC(k+1) = -2*(N-k)*L*log(AIC_temp1)+2*k*(2*N-k); %一次实验，采样次数为1
    MDL(k+1) = -L*(N-k)*log(AIC_temp1)+k*(2*N-k)*log(L); %一次实验，采样次数为1
end
AIC_abs = abs(AIC);
fig5 = figure('visible','off');
plot(AIC_abs);
hold on;
plot(abs(MDL));
saveas(fig5, '信源数.png');

%P = 2; %极点数，应该是要用算法求的
AIC_abs=abs(AIC);
P2 = length(AIC_abs(AIC_abs>mean(AIC_abs))); %另一种对目标数的估计方法，赋值为P2
% MDC也可以同样方法赋值
%% MUSIC
v1(:,1:N-P) = V(:,1:N-P); %特征向量取列向量
v2(:,N-P+1:N) = V(:,N-P+1:N);
NewR = Rxx;
for i = 1:length(d)
    NewR  = NewR - D(i,i)*V(:,i)*V(:,i)';
end

MUSIC = zeros(1,length(scan)); % MUSIC的功率谱
for i = 1:length(scan)
    omega2 = 2*pi*dr*sin(scan(i))/Lambda;
    a_send2 = exp(1i*-(0:1:M-1)'*omega2);
    a_back2 = exp(1i*-(0:1:N-1)'*omega2);
    
    MUSIC(i) = 1/(a_back2'*v1*v1'*a_back2);
end
fig6 = figure('visible','off');
plot(scan,log(abs(MUSIC)));
title('MUSIC');
saveas(fig6, 'MUSIC.png');