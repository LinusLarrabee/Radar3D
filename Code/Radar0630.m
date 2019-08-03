clear;
PowerK = 1000;
%% 目标参数
global C Lambda PulseNum BandWidth PRT PRF Fs NoisePower Fc ...
    WaveNum deltaD MinDis deltaAngle CarAngle
Parameter;
%% 环境参数配置及读取
global Point PointNum Swirling
Swirling = 1;
PointNum = 3;
Point = CarSet(PointNum);
Environment;

%% 信号参数-产生线性调频信号

%[Chirp,Reality] = SigGeneration();
[Chirp,Reality] = SigGeneration(1);
%% 脉冲压缩系数
if Reality ==0

    coeff=conj(fliplr(Chirp)); %把Chirp矩阵翻转取共轭，产生脉压系数
    % 简单解释一下就是这里取的是 h(t) = x(^*)(t_0-t) 中 t_0 = 0 的情景
    figure(13);
    subplot(2,1,1),plot(real(coeff));
    subplot(2,1,2),plot(imag(coeff));
elseif Reality ==1
    %实信号码?调制
    LocalIQNum = 0:fix(WaveNum)-1;
    M = 131126;     %131126点fft
    local_oscillator_i=cos(LocalIQNum*Fc/Fs*2*pi);%i路本振信号
    local_oscillator_q=sin(LocalIQNum*Fc/Fs*2*pi);%q路本振信号
    fbb_i=local_oscillator_i.*Chirp;%i路解调   先进行一个码元的求解脉冲压缩系数
    fbb_q=local_oscillator_q.*Chirp;%q路解调
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
    coeff=conj(fliplr(fbb)); %把Chirp矩阵翻转取共轭，产生脉压系数
    figure(14);subplot(2,1,1),plot(real(fbb_fft_result));
    subplot(2,1,2),plot(imag(fbb_fft_result));
    figure(3);subplot(2,1,1),plot(fbb_i);
    xlabel('t(单位：秒)');title('雷达发射信号码内解调后I路信号');
    subplot(2,1,2),plot(fbb_q);
    xlabel('t(单位：秒)');title('雷达发射信号码内解调后Q路信号');
    figure(4)
    plot((0:Fs/WaveNum:Fs/2-Fs/WaveNum),abs(fbb_fft_result(1:WaveNum/2)));
    xlabel('频率f(单位 Hz)');title('雷达发射信号码内解调信号的频谱');
end
%% 相关参数计算
load('Environment.mat');
[AngleNum,MaxDis] = size(Envir);
TargetDistance = deltaD*(0.5:1:MaxDis-0.5);
DelayNum = fix(Fs*2*TargetDistance/C);% 把目标距离换算成采样点（距离门） 
%fix函数向0靠拢取整
Sector = round((CarAngle(1) + deltaAngle*(1:AngleNum))*180/pi);
SampleNum=fix(Fs*PRT);%计算一个脉冲周期的采样点数；
AngleCirNum = SampleNum*AngleNum; 
%先完成一个角度的一整个线性调频，再依次完成剩余角度，最后循环脉冲；
TotalNum=SampleNum*PulseNum*AngleNum;%总的采样点数；
% BlindNum=fix(Fs*TimeWidth);%计算一个脉冲周期的盲区-遮挡样点数；
%% 产生目标回波串
if Reality >= 0
    Signal = zeros(1,TotalNum);%所有脉冲的信号,先填0
    for Inn = 1:PulseNum % 16个回波脉冲
        SignalAng = zeros(1,AngleCirNum);
        for ang = 1:AngleNum
            detect = find (Envir(ang,:) ~= 0);
            [~,detectN] = size(detect);
            if detectN>=MinDis
                detectAct = detect(1,MinDis:detectN);
                for k = detectAct % 依次产生各个目标
                    %fi=2*pi/10 * fix(10*rand);
                    SignalTemp=zeros(1,SampleNum);% 一个PRT
                    Power = PowerK*sqrt(Envir(ang,k)*TargetDistance(k)^(-4));
                    SignalTemp(DelayNum(k)+1:DelayNum(k)+WaveNum)=...
                         Power*Chirp;%*exp(1i*fi)
                    StatusNum =(ang-1)*SampleNum+(Inn-1)*AngleCirNum;
                    if Reality ==0
                        FreqMove=exp(1i*4*pi*Velocity(ang,k)*(StatusNum:1:StatusNum+SampleNum-1)/Fs/Lambda);
                    else 
                        FreqMove=cos(4*pi*Velocity(ang,k)*(StatusNum:1:StatusNum+SampleNum-1)/Fs/Lambda);
                    end
                    SignalTemp = SignalTemp.*FreqMove;
                    %一个脉冲的1个目标（加多普勒速度）
                    %(DelayNum(k)+1):(DelayNum(k)+WaveNum)
                    SignalAng((ang-1)*SampleNum+1:ang*SampleNum)=SignalTemp;
                end
            end
            if mod(Swirling,2) ==0
                Upgrade(PRT);
            end
        end
        Signal((Inn-1)*AngleCirNum+1:Inn*AngleCirNum)=SignalAng;
        if mod(Swirling,2) ==1
            Upgrade(PRT*AngleNum);
            % 仅仅考虑Swirling的运动分布而不考虑其概率分布
        end
    end

    figure(2);
    subplot(2,1,1);plot(real(Signal),'r-');title('目标信号的实部');...
    grid on;zoom on;
    subplot(2,1,2);plot(imag(Signal));title('目标信号的虚部');grid on;zoom on;

    %% 总的回波信号
    % 产生系统噪声信号
    SystemNoise = normrnd(0,10^(NoisePower/10),1,TotalNum)...
        +1i*normrnd(0,10^(NoisePower/10),1,TotalNum);
    %均值为0，标准差为10^(NoisePower/10)的噪声

    %闭锁期无回波
    EchoAll=Signal+SystemNoise;% +SeaClutter+TerraClutter，加噪声之后的回波
    for i=1:PulseNum*AngleNum   %在接收机闭锁期,接收的回波为0
        EchoAll((i-1)*SampleNum+1:(i-1)*SampleNum+WaveNum)=0; %发射时接收为0
    end
    f3 = figure(3);%加噪声之后的总回波信号
    subplot(2,1,1);plot(real(EchoAll),'r-');title('总回波信号的实部,闭锁期为0');
    subplot(2,1,2);plot(imag(EchoAll));title('总回波信号的虚部,闭锁期为0');
    saveas(f3,'figure3.jpg')
else
    %% 实信号回波生成
    Signal = zeros(1,TotalNum);%所有脉冲的信号,先填0
    for Inn = 1:PulseNum % 16个回波脉冲
        SignalAng = zeros(1,AngleCirNum);
        for ang = 1:AngleNum
            detect = find (Envir(ang,:) ~= 0);
            [~,detectN] = size(detect);
            if detectN>=MinDis
                detectAct = detect(1,MinDis:detectN);
                for k = detectAct % 依次产生各个目标
                    fi=2*pi/10 * fix(10*rand);
                    SignalTemp=zeros(1,SampleNum);% 一个PRT
                    Power = sqrt(Envir(ang,k)*TargetDistance(k)^(-4));
                    SignalTemp(DelayNum(k)+1:DelayNum(k)+WaveNum)=...
                        Power*cos(fi)*Chirp;
                    StatusNum =(ang-1)*SampleNum+(Inn-1)*AngleCirNum;
                    FreqMove=cos(4*pi*Velocity(ang,k)*(StatusNum:1:StatusNum+SampleNum-1)/Fs/Lambda);
                    SignalTemp = SignalTemp.*FreqMove;
                    %一个脉冲的1个目标（加多普勒速度）
                    %(DelayNum(k)+1):(DelayNum(k)+WaveNum)
                    SignalAng((ang-1)*SampleNum+1:ang*SampleNum)=SignalTemp;
                end
            end

            if mod(Swirling,2) ==0
                Upgrade(PRT);
            end
        end
        Signal((Inn-1)*AngleCirNum+1:Inn*AngleCirNum)=SignalAng;
        if mod(Swirling,2) ==1
            Upgrade(PRT*AngleNum);
            % 仅仅考虑Swirling的运动分布而不考虑其概率分布
        end
    end

    figure(7);plot(Signal,'r-');title('实信号回波');
    grid on;zoom on;

    % 产生系统噪声信号
    SystemNoise = normrnd(0,10^(NoisePower/10),1,TotalNum);
    %均值为0，标准差为10^(NoisePower/10)的噪声

    %闭锁期无回波
    EchoAll=Signal+SystemNoise;% +SeaClutter+TerraClutter，加噪声之后的回波
    for i=1:PulseNum*AngleNum   %在接收机闭锁期,接收的回波为0
        EchoAll((i-1)*SampleNum+1:(i-1)*SampleNum+WaveNum)=0; %发射时接收为0
    end
    figure(8);%加噪声之后的总回波信号
    plot(EchoAll,'r-');title('带杂波的实信号回波,闭锁期为0');
end
SignalProcess;