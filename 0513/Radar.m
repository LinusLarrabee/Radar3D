%% Initial and introduction of the process
% Radar3D.m
% Created by Group 6 on 5/8/19.
close all; clear ; clc;

%% Radar Parameter
C=3.0e8;  %光速(m/s)
RF=30e9;  %雷达射频 1.57GHz
M = 524288;     %131126点fft
PulseNum=16;   %回波脉冲数
BandWidth=4.0e9;  %发射信号带宽 带宽B=1/τ，τ是脉冲宽度
TimeWidth=4.0e-8; %发射信号时宽
PRT=8e-7;   % 雷达发射脉冲重复周期(s),240us对应1/2*240*300=36000米最大无模糊距离
PRF=1/PRT;
Fs=2.5e10;  %采样频率
NoisePower=-20;%(dB);%噪声功率（目标为0dB）
Fc = 30e9;
Lambda=C/Fc;    %雷达工作波长

%% 目标参数
Target = xlsread ('Info.xlsx');
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

figure(1);
subplot(2,1,1),plot(RealChirp);
subplot(2,1,2);
plot((0:Fs/WaveNum:Fs-Fs/WaveNum),abs(fft(Chirp)));
figure(2);
subplot(2,1,1);
fftR = fft(RealChirp);
subplot(2,1,1);
plot((0:Fs/WaveNum:Fs-Fs/WaveNum),abs(fftR));
subplot(2,1,2);
plot((0:Fs/WaveNum:Fs/2-Fs/WaveNum),abs(fftR(1:WaveNum/2)));

coeff=conj(fliplr(Chirp)); %把Chirp矩阵翻转取共轭，产生脉压系数
% 简单解释一下就是这里取的是 h(t) = x(^*)(t_0-t) 中 t_0 = 0 的情景

%% 码内正交解调  
LocalIQNum = 0:fix(WaveNum)-1;
local_oscillator_i=cos(LocalIQNum*Fc/Fs*2*pi);%i路本振信号
local_oscillator_q=sin(LocalIQNum*Fc/Fs*2*pi);%q路本振信号
fbb_i=local_oscillator_i.*RealChirp;%i路解调   先进行一个码元的求解脉冲压缩系数
fbb_q=local_oscillator_q.*RealChirp;%q路解调
window=chebwin(51,40); %切比雪夫窗函数
[b,a]=fir1(50,2*BandWidth/Fs,window); %b a 分别是经滤波器定义的分子分母系数向量
fbb_i=[fbb_i,zeros(1,25)];  %因为该FIR滤波器又25个采样周期的延迟，为了保证所有的有效信息全部通过滤波器，因此在信号后面扩展了25个0
fbb_q=[fbb_q,zeros(1,25)];
fbb_i=filter(b,a,fbb_i);   %I 路 Q路信号经过低通滤波器
fbb_q=filter(b,a,fbb_q);
fbb_i=fbb_i(26:end);%截取有效信息
fbb_q=fbb_q(26:end);%截取有效信息  
fbb=fbb_i+1i*fbb_q;
fbb_fft_result = fft(fbb);

figure(3);subplot(2,1,1),plot(fbb_i);
xlabel('t(单位：秒)');title('雷达发射信号码内解调后I路信号');
subplot(2,1,2),plot(fbb_q);
xlabel('t(单位：秒)');title('雷达发射信号码内解调后Q路信号');
figure(4)
plot((0:Fs/WaveNum:Fs/2-Fs/WaveNum),abs(fbb_fft_result(1:WaveNum/2)));
xlabel('频率f(单位 Hz)');title('雷达发射信号码内解调信号的频谱');


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

figure(5);
subplot(2,1,1);plot(real(SignalAll),'r-');title('目标信号的实部');...
grid on;zoom on;
subplot(2,1,2);plot(imag(SignalAll));title('目标信号的虚部');grid on;zoom on;


%% 总的回波信号
% 产生系统噪声信号
SystemNoise = normrnd(0,10^(NoisePower/10),1,TotalNum)...
    +1i*normrnd(0,10^(NoisePower/10),1,TotalNum);
%均值为0，标准差为10^(NoisePower/10)的噪声

%闭锁期无回波
EchoAll=SignalAll+SystemNoise;% +SeaClutter+TerraClutter，加噪声之后的回波
for i=1:PulseNum*AngleNum   %在接收机闭锁期,接收的回波为0
    EchoAll((i-1)*SampleNum+1:(i-1)*SampleNum+WaveNum)=0; %发射时接收为0
end
figure(6);%加噪声之后的总回波信号
subplot(2,1,1);plot(real(EchoAll),'r-');title('总回波信号的实部,闭锁期为0');
subplot(2,1,2);plot(imag(EchoAll));title('总回波信号的虚部,闭锁期为0');

%% 实信号回波生成
RealSignalAll=zeros(1,TotalNum);%所有脉冲的信号,先填0
for ang = 1:AngleNum
    detect = find (TargetAngle == Sector(ang));
    if ~isempty(detect)
        for k = detect % 依次产生各个目标
            fi=2*pi/10 * fix(10*rand);
            RealSignalTemp=zeros(1,SampleNum);% 一个PRT
            RealSignalTemp(DelayNum(k)+1:DelayNum(k)+WaveNum)=...
                sqrt(SigPower(k))*cos(fi)*RealChirp;
            %一个脉冲的1个目标（未加多普勒速度）(DelayNum(k)+1):(DelayNum(k)+WaveNum)
            RealSignal = zeros(1,TotalNum);
            RealSignalAng = zeros(1,AngleCirNum);
            RealSignalAng((ang-1)*SampleNum+1:ang*SampleNum)=RealSignalTemp;
            for i=1:PulseNum % 16个回波脉冲
                RealSignal((i-1)*AngleCirNum+1:i*AngleCirNum)=RealSignalAng;
                %每个目标把16个RealSignalTemp排在一起
            end
            RealFreqMove=cos(2*pi*TargetFd(k)*(0:TotalNum-1)/Fs);
            %目标的多普勒速度*时间=目标的多普勒相移
            RealSignal=RealSignal.*RealFreqMove;%加上多普勒速度后的16个脉冲1个目标
            RealSignalAll=RealSignalAll+RealSignal;%加上多普勒速度后的16个脉冲4个目标
        end
    end
end

figure(7);plot(RealSignalAll,'r-');title('实信号回波');...
grid on;zoom on;

% 产生系统噪声信号
RealSystemNoise = normrnd(0,10^(NoisePower/10),1,TotalNum);
%均值为0，标准差为10^(NoisePower/10)的噪声

%闭锁期无回波
RealEchoAll=RealSignalAll+RealSystemNoise;% +SeaClutter+TerraClutter，加噪声之后的回波
for i=1:PulseNum*AngleNum   %在接收机闭锁期,接收的回波为0
    RealEchoAll((i-1)*SampleNum+1:(i-1)*SampleNum+WaveNum)=0; %发射时接收为0
end
figure(8);%加噪声之后的总回波信号
plot(RealEchoAll,'r-');title('带杂波的实信号回波,闭锁期为0');

%% 文件生成
Folder = {'总回波','时域脉压','时频域脉压对比','频域脉压幅度','MTI','MTD','cfar'};
for i = 1:length(Folder)
    if exist(char(Folder(i)),'dir')
        rmdir(char(Folder(i)),'s')
    end
    mkdir(char(Folder(i)))
end
DetectDV = zeros(AngleNum,6);

EchoRoute = reshape(RealEchoAll, [SampleNum,AngleNum,PulseNum]);
for argerich =1 : AngleNum
    Echo = reshape(EchoRoute(:,argerich,:),1,[]);
    DetectDV(argerich,4:6) = [Sector(argerich),TargetVelocity(argerich),TargetDistance(argerich)];
    [DetectDV(argerich,1:3)] = DSP( Echo,coeff,M);
%     hugo = figure('visible','off');
%     subplot(3,1,1)
%     plot(real(Echo),'r-');title('总回波信号的实部,闭锁期为0');
%     filename=['总回波/扫描' num2str(Sector(argerich)) '度.png'];
%     saveas(hugo,filename)
%     close(gcf)
%     %% 实信号接收调制
% 
%     EchoNum=SampleNum*PulseNum;
%     n=0:EchoNum-1;
%     receiver_oscillator_i=cos(n*Fc/Fs*pi);%I路本振信号
%     receiver_oscillator_q=cos(n*Fc/Fs*pi);%Q路本振信号
%     s_echo_i=receiver_oscillator_i.* Echo;%I路解调
%     s_echo_q=receiver_oscillator_q.* Echo;%Q路解调
%     receiverwindow=chebwin(51,40);%这是采50阶cheby窗的FIR低通滤波器
%     [Rb,Ra]=fir1(50,2*BandWidth/Fs,receiverwindow);
%     s_echo_i=[s_echo_i,zeros(1,25)];
%     s_echo_q=[s_echo_q,zeros(1,25)];
%     s_echo_i=filter(Rb,Ra,s_echo_i);
%     s_echo_q=filter(Rb,Ra,s_echo_q);
%     s_echo_i=s_echo_i(26:end);%截取有效信息
%     s_echo_q=s_echo_q(26:end);%截取有效信息
%     s_echo_mf=s_echo_i+1i*s_echo_q;
%     figure(9)
%     subplot(2,1,1),plot((0:1/Fs:PulseNum*PRT-1/Fs),s_echo_i);
%     xlabel('t(unit:s)'); title('雷达回波信号解调后的I路信号');
% 
%     subplot(2,1,2),plot((0:1/Fs:PulseNum*PRT-1/Fs),s_echo_q);
%     xlabel('t(unit:s)'); title('雷达回波信号解调后的q路信号');
% 
%     %% 时域脉压
%     pc_time0=conv(Echo,coeff);%pc_time0为Echo和coeff的卷积
%     pc_time1=pc_time0(WaveNum:length(Echo)+WaveNum-1);%去掉暂态点 WaveNum-1个
%     %figure(4);%时域脉压结果的幅度
%     subplot(3,1,2);plot(abs(pc_time0),'r-');title('时域脉压结果的幅度,有暂态点');
%     %pc_time0的模的曲线
%     subplot(3,1,3);plot(abs(pc_time1));title('时域脉压结果的幅度,无暂态点');
%     %pc_time1的模的曲线
%     filename=['时域脉压/扫描' num2str(Sector(argerich)) '度.png'];
%     saveas(hugo,filename)
%     close(gcf)
% 
%     %% 频域脉压
%     Echo_fft=fft(Echo,524288);
%     %理应进行length(Echo)+WaveNum-1点FFT,但为了提高运算速度,进行了8192点的FFT
%     coeff_fft=fft(coeff,524288);
%     pc_fft=Echo_fft.*coeff_fft;
%     pc_freq0=ifft(pc_fft);
%     hug1 = figure('visible','off');
%     subplot(2,1,1);plot(abs(pc_freq0(1:length(Echo)+WaveNum-1)));
%     title('频域脉压结果的幅度,有前暂态点');
%     subplot(2,1,2);
%     plot(abs(pc_time0(1:length(Echo)+WaveNum-1)-...
%             pc_freq0(1:length(Echo)+WaveNum-1)),'r');
%     title('时域和频域脉压的差别');
%     filename=['时频域脉压对比/扫描' num2str(Sector(argerich)) '度.png'];
%     saveas(hug1,filename)
%     close(gcf)
% 
%     pc_freq1=pc_freq0(WaveNum:length(Echo)+WaveNum-1);
%     %去掉暂态点 WaveNum-1个,后填充点若干(8192-WaveNum+1-length(Echo))
%     % ================按照脉冲号、距离门号重排数据=================================%
%     pc = zeros(PulseNum,SampleNum);
%     for i=1:PulseNum
%         pc(i,1:SampleNum)=pc_freq1((i-1)*SampleNum+1:i*SampleNum);
%         %每个PRT为一行，每行480个采样点的数据
%     end
%     hug2 = figure('visible','off');
%     plot(abs(pc(1,:)));title('频域脉压结果的幅度,没有暂态点');
%     filename=['频域脉压幅度/扫描' num2str(Sector(argerich)) '度.png'];
%     saveas(hug2,filename)
%     close(gcf)
% %     %% 脉冲压缩 
% %     coeff_fft=fft(coeff,M);
% %     for i=1:PulseNum
% %         s_echo_fft_result=fft(s_echo_mf(1,(i-1)*PRT*Fs+1:i*PRT*Fs),M);
% %         s_pc_fft=s_echo_fft_result.*coeff_fft;
% %         pc(i,:)=ifft(s_pc_fft,M);    
% %     end
% % 
% %     hug1 = figure('visible','off');
% %     plot(abs(pc(1,:)));%一个周期内三峰值点理论上为40 107 304
% % 
% %     filename=['频域脉压幅度/扫描' num2str(Sector(argerich)) '度.png'];
% %     saveas(hug1,filename)
% %     close(gcf)
% % 
% %     s_pc_result_1=pc;
% %     s_pc_result_1=reshape((s_pc_result_1)',1,PulseNum*M);   %%%%%%%%%%注意，这里由于reshape函数的算法，需要将矩阵转置才能首尾连在一起
% %     hug2 = figure('visible','off');
% %     subplot(2,1,1),plot((0:1/Fs:PulseNum*M/Fs-1/Fs),abs(s_pc_result_1)),
% %     %N_echo_frame*T_frame-ts
% %     xlabel('t(单位:s)'),title('脉冲压缩处理后结果（实部）');
% %     subplot(2,1,2),plot((0:1/Fs:PulseNum*M/Fs-1/Fs),imag(s_pc_result_1)),
% %     xlabel('t(单位:s)'),title('脉冲压缩处理后结果（虚部）');
% %     filename=['频域脉压幅度/扫描' num2str(Sector(argerich)) '度.png'];
% %     saveas(hug2,filename)
% %     close(gcf)
%     %% MTI（动目标显示）,对消静止目标和低速目标---可抑制杂波%
%     %mti = zeros(PulseNum-1,SampleNum);
%     for i=1:PulseNum-1  %滑动对消，少了一个脉冲
%         mti(i,:)=pc(i+1,:)-pc(i,:);
%     end
%     hug3 = figure('visible','off');
%     mesh(abs(mti));title('MTI  result');
%     filename=['MTI/扫描' num2str(Sector(argerich)) '度.png'];
%     saveas(hug3,filename)
%     close(gcf)
%     %% MTD（动目标检测）,区分不同速度的目标，有测速作用==%
%     mtd_PulseNum = PulseNum-1;
%     mtd=zeros(PulseNum,SampleNum);
%     for i=1:SampleNum
%         buff(1:(mtd_PulseNum))= mti(1:(mtd_PulseNum),i);
%         buff_fft=fftshift(fft(buff)); %用fftshift将零频搬移到中间 这样可以方便观察速度正负
%         mtd(1:mtd_PulseNum,i)=buff_fft(1:mtd_PulseNum)';
%     end
%     hug4 = figure('visible','off');
%     mesh(abs(mtd));title('MTD  result');
%     filename=['MTD/扫描' num2str(Sector(argerich)) '度.png'];
%     saveas(hug4,filename)
%     close(gcf);
% %     %%
% %     mtd=zeros(PulseNum,SampleNum);
% %     for i=1:SampleNum
% %        buff(1:PulseNum)=pc(1:PulseNum,i);
% %        buff_fft=fft(buff);
% %        mtd(1:PulseNum,i)=buff_fft(1:PulseNum);
% %     end
% %       figure(8);mesh(abs(mtd));title('MTD  result');
%     %% 找到多普勒——距离地图中峰值最高的目标 
%     abs_mtd = abs(mtd); %对正交分量的复数求模
%     max_target = max(max(abs_mtd));
%     [row,cell] = find( abs_mtd == max_target);
%     target_D = cell/Fs*C/2;
%     target_V = ((abs(row - 8))/PulseNum)*PRF*(Lambda/2);
%     %% 二维恒虚警 
%     ex_PulseNum = mtd_PulseNum+4;
%     ex_SampleNum = SampleNum+4;
%     ex_mtd = zeros(ex_PulseNum ,ex_SampleNum);
%     cfar_mtd = zeros(ex_PulseNum ,ex_SampleNum);
%     T = 1.0;     %恒虚警CFAR 阈值因子
% 
%     for i = 3:(ex_PulseNum-2)
%         for j = 3:(ex_SampleNum-2)         
%               ex_mtd(i,j) = abs_mtd(i-2,j-2); %构建长宽每边都扩张两格的矩阵供恒虚警窗口检测     
%         end
%     end
% 
%     for i = 3:(ex_PulseNum-2)
%         for j = 3:(ex_SampleNum-2)    %从有目标的窗口中进行检测
%             cfar_sx = ex_mtd(i-2,j-2) + ex_mtd(i-2,j-1) +ex_mtd(i-1,j-2) +ex_mtd(i,j-2)+ex_mtd(i+1,j-2)+ex_mtd(i+2,j-2) + ex_mtd(i+2,j-1);                   
%             cfar_sy = ex_mtd(i-2,j+1) + ex_mtd(i-2,j+2) +ex_mtd(i-1,j+2) +ex_mtd(i,j+2)+ex_mtd(i+1,j+2)+ex_mtd(i+2,j+1) + ex_mtd(i+2,j+2);
%             cfar_s = T*min(cfar_sx, cfar_sy);      %使用GO—CFAR二维恒虚警，
%             if (ex_mtd(i,j)>= cfar_s)
%                 cfar_mtd(i,j) = ex_mtd(i,j);
%             end
%         end
%     end
%     figure(100);mesh(abs(cfar_mtd));title('cfar  result');  
% %     hug5 = figure('visible','off');
% %     mesh(abs(mtd));title('MTD  result');
% %     filename=['cfar/扫描' num2str(Sector(argerich)) '度.png'];
% %     saveas(hug5,filename)
% %     close(gcf);
end
xlswrite('DV.xlsx',DetectDV);
