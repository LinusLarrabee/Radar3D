%% 回波信号整形
Folder = {'时域脉压','时频域脉压对比','频域脉压幅度','MTI','MTD','ReceiverIQ'};
for i = 1:length(Folder)
    if exist(char(Folder(i)),'dir')
        rmdir(char(Folder(i)),'s')
    end
    mkdir(char(Folder(i)))
end

EchoRoute = reshape(RealEchoAll, [SampleNum,AngleNum,PulseNum]);
LJCDis = zeros(2,40);
tag = 1;
for argerich =  1 : AngleNum
    Echo = reshape(EchoRoute(:,argerich,:),1,[]);
    %% 实信号接收调制
    EchoNum=SampleNum*PulseNum;
    n=0:EchoNum-1;
    receiver_oscillator_i=cos(n*Fc/Fs*pi);%I路本振信号
    receiver_oscillator_q=cos(n*Fc/Fs*pi);%Q路本振信号
    s_echo_i=receiver_oscillator_i.* Echo;%I路解调
    s_echo_q=receiver_oscillator_q.* Echo;%Q路解调
    receiverwindow=chebwin(51,40);%这是采50阶cheby窗的FIR低通滤波器
    [Rb,Ra]=fir1(50,2*BandWidth/Fs,receiverwindow);
    s_echo_i=[s_echo_i,zeros(1,25)];
    s_echo_q=[s_echo_q,zeros(1,25)];
    s_echo_i=filter(Rb,Ra,s_echo_i);
    s_echo_q=filter(Rb,Ra,s_echo_q);
    s_echo_i=s_echo_i(26:end);%截取有效信息
    s_echo_q=s_echo_q(26:end);%截取有效信息
    Echo=s_echo_i+1i*s_echo_q;
    hug_1 = figure('visible','off');
    subplot(2,1,1),plot((0:1/Fs:PulseNum*PRT-1/Fs),s_echo_i);
    xlabel('t(unit:s)'); title('雷达回波信号解调后的I路信号');
    subplot(2,1,2),plot((0:1/Fs:PulseNum*PRT-1/Fs),s_echo_q);
    xlabel('t(unit:s)'); title('雷达回波信号解调后的q路信号');
    filename=['ReceiverIQ/扫描' num2str(Sector(argerich)) '度.png'];
    saveas(hug_1,filename)
    %close(gcf)
     
    %% 时域脉压
    hugo = figure('visible','off');
    subplot(3,1,1)
    plot(real(Echo),'r-');title('总回波信号的实部,闭锁期为0'); 
    pc_time0=conv(Echo,coeff);%pc_time0为Echo和coeff的卷积
    pc_time1=pc_time0(WaveNum:length(Echo)+WaveNum-1);%去掉暂态点 WaveNum-1个
    %figure(4);%时域脉压结果的幅度
    subplot(3,1,2);plot(abs(pc_time0),'r-');title('时域脉压结果的幅度,有暂态点');
    %pc_time0的模的曲线
    subplot(3,1,3);plot(abs(pc_time1));title('时域脉压结果的幅度,无暂态点');
    %pc_time1的模的曲线
    filename=['时域脉压/扫描' num2str(Sector(argerich)) '度.png'];
    saveas(hugo,filename)
    %close(gcf)
    
    %% 频域脉压
    Echo_fft=fft(Echo,length(Echo)+WaveNum-1);
    %理应进行length(Echo)+WaveNum-1点FFT,但为了提高运算速度,进行了8192点的FFT
    coeff_fft=fft(coeff,length(Echo)+WaveNum-1);
    pc_fft=Echo_fft.*coeff_fft;
    pc_freq0=ifft(pc_fft);
    hug1 = figure('visible','off');
    subplot(2,1,1);plot(abs(pc_freq0(1:length(Echo)+WaveNum-1)));
    title('频域脉压结果的幅度,有前暂态点');
    subplot(2,1,2);
    plot(abs(pc_time0(1:length(Echo)+WaveNum-1)-...
            pc_freq0(1:length(Echo)+WaveNum-1)),'r');
    title('时域和频域脉压的差别');
    filename=['时频域脉压对比/扫描' num2str(Sector(argerich)) '度.png'];
    saveas(hug1,filename)
%        close(gcf)
    
    pc_freq1=pc_freq0(WaveNum:length(Echo)+WaveNum-1);
    %去掉暂态点 WaveNum-1个,后填充点若干(8192-WaveNum+1-length(Echo))
    %% 按照脉冲号、距离门号重排数据=================================%
    for i=1:PulseNum
        pc(i,1:SampleNum)=pc_freq1((i-1)*SampleNum+1:i*SampleNum);
        %每个PRT为一行，每行480个采样点的数据
    end
    hug2 = figure('visible','off');
    plot(abs(pc(1,:)));title('频域脉压结果的幅度,没有暂态点');
    filename=['频域脉压幅度/扫描' num2str(Sector(argerich)) '度.png'];
    saveas(hug2,filename)
%        close(gcf)
    %% MTI（动目标显示）,对消静止目标和低速目标---可抑制杂波%
    for i=1:PulseNum-1  %滑动对消，少了一个脉冲
        mti(i,:)=pc(i+1,:)-pc(i,:);
    end
    hug3 = figure('visible','off');
    mesh(abs(mti));title('MTI  result');
    filename=['MTI/扫描' num2str(Sector(argerich)) '度.png'];
    saveas(hug3,filename)
%        close(gcf)
    %% MTD（动目标检测）,区分不同速度的目标，有测速作用
    mtd=zeros(PulseNum,SampleNum);
    for i=1:SampleNum
        buff(1:PulseNum)=pc(1:PulseNum,i);
        buff_fft=fft(buff);
        mtd(1:PulseNum,i)=buff_fft(1:PulseNum);
    end
    hug4 = figure('visible','off');
    mesh(abs(mtd));title('MTD  result');
    filename=['MTD/扫描' num2str(Sector(argerich)) '度.png'];
    saveas(hug4,filename)
%        close(gcf)
    %%  找到多普勒——距离地图中峰值最高的目标
    abs_mtd = abs(mtd); %对正交分量的复数求模
    max_target = max(max(abs_mtd));
    [row,cell] = find( abs_mtd == max_target);
    target_D = cell/Fs*C/2
    target_V = (fix(row - 1)/PulseNum)*PRF*Lambda/2
    
         %% 二维恒虚警
    
    ex_PulseNum = PulseNum+4;
    ex_SampleNum = SampleNum+4;
    ex_mtd = zeros(ex_PulseNum ,ex_SampleNum);
    cfar_mtd = zeros(ex_PulseNum ,ex_SampleNum);
    T = 0.25;     %恒虚警CFAR 阈值因子

    for i = 3:(ex_PulseNum-2)
        for j = 3:(ex_SampleNum-2)

              ex_mtd(i,j) = abs_mtd(i-2,j-2); %构建长宽每边都扩张两格的矩阵供恒虚警窗口检测

        end
    end

    for i = 3:(ex_PulseNum-2)
        for j = 3:(ex_SampleNum-2)    %从有目标的窗口中进行检测
            cfar_sx = ex_mtd(i-2,j-2) + ex_mtd(i-2,j-1) +ex_mtd(i-1,j-2) +ex_mtd(i,j-2)+ex_mtd(i+1,j-2)+ex_mtd(i+2,j-2) + ex_mtd(i+2,j-1);                   
            cfar_sy = ex_mtd(i-2,j+1) + ex_mtd(i-2,j+2) +ex_mtd(i-1,j+2) +ex_mtd(i,j+2)+ex_mtd(i+1,j+2)+ex_mtd(i+2,j+1) + ex_mtd(i+2,j+2);
            cfar_s = T*max(cfar_sx, cfar_sy);      %使用GO?CFAR二维恒虚警，
            if (ex_mtd(i,j)>= 1200)
                cfar_mtd(i,j) = ex_mtd(i,j);
            end
        end
    end
     x=0:1:ex_SampleNum-1;
     y=-9:1:9;%通道这样设后读出的通道数乘单位值则是速度值。
    % figure(9);mesh(x,y,abs(cfar_mtd));title('cfar  result');       
        figure(argerich+10);mesh(abs(cfar_mtd));title('cfar  result'); 
     %% 恒虚警后续筛选处理
    result_target = zeros(ex_PulseNum,ex_SampleNum);
    flag = 0;
    for i = 1:ex_PulseNum
        for j = 1:ex_SampleNum

            if ((cfar_mtd(i,j) ~= 0)&&(flag == 0))
                result_target(i,j) = cfar_mtd(i,j);
                flag = 1500;
            end
            if(flag ~= 0)
                flag = flag - 1;
            end
        end
    end
    % x=0:1:ex_SampleNum-1;
    % y=-9:1:9;%通道这样设后读出的通道数乘单位值则是速度值。
    figure(argerich+50);mesh(abs(result_target));title('result');    
    %% 找到目标
        [r,c] = find(result_target>1)
        target_DJ = (c+1500)/Fs*C/2
    for index = 1:length(r)/2
        LJCDis(1,tag) = target_DJ(index)
        LJCDis(2,tag) = Sector(argerich)
        tag = tag+1;
    end
        
end

%% hmmm
%result_tar = [target_D',target_V',target_row',Sector'];
f9 = figure(900)
polarscatter(deg2rad(LJCDis(2,:)),LJCDis(1,:),'filled')
thetalim([30 150])
f10 = figure(1000)
polarscatter(deg2rad(TargetAngle),TargetDistance,'filled')
thetalim([30 150])
f11 = figure(1100)
polarscatter(deg2rad(LJCDis(2,:)),LJCDis(1,:),'filled')
thetalim([30 150])
hold on
polarscatter(deg2rad(TargetAngle),TargetDistance,'filled')

%csvwrite('result.csv',result_tar)

%saveas(f9,'cfar检测场景所有目标（figure9）.fig')
