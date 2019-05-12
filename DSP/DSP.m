function [ target_V ] = DSP( RealEchoAll,SampleNum,AngleNum,PulseNum,Fc,Fs,BandWidth,PRT,coeff,M )
    EchoRoute = reshape(RealEchoAll, [SampleNum,AngleNum,PulseNum]);
    for argerich = 10:10:40%1 : AngleNum
        Echo = reshape(EchoRoute(:,argerich,:),1,[]);
        hugo = figure('visible','off');
        subplot(3,1,1)
        plot(real(Echo),'r-');title('总回波信号的实部,闭锁期为0');
        filename=['总回波/扫描' num2str(Sector(argerich)) '度.png'];
        saveas(hug1,filename)
        close(gcf)
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
        s_echo_mf=s_echo_i+1i*s_echo_q;
        figure(9)
        subplot(2,1,1),plot((0:1/Fs:PulseNum*PRT-1/Fs),s_echo_i);
        xlabel('t(unit:s)'); title('雷达回波信号解调后的I路信号');

        subplot(2,1,2),plot((0:1/Fs:PulseNum*PRT-1/Fs),s_echo_q);
        xlabel('t(unit:s)'); title('雷达回波信号解调后的q路信号');

%         %% 时域脉压
%         pc_time0=conv(Echo,coeff);%pc_time0为Echo和coeff的卷积
%         pc_time1=pc_time0(WaveNum:length(Echo)+WaveNum-1);%去掉暂态点 WaveNum-1个
%         %figure(4);%时域脉压结果的幅度
%         subplot(3,1,2);plot(abs(pc_time0),'r-');title('时域脉压结果的幅度,有暂态点');
%         %pc_time0的模的曲线
%         subplot(3,1,3);plot(abs(pc_time1));title('时域脉压结果的幅度,无暂态点');
%         %pc_time1的模的曲线
%         filename=['时域脉压/扫描' num2str(Sector(argerich)) '度.png'];
%         saveas(hugo,filename)
%         close(gcf)
%         
%         %% 频域脉压
%         Echo_fft=fft(Echo,524288);
%         %理应进行length(Echo)+WaveNum-1点FFT,但为了提高运算速度,进行了8192点的FFT
%         coeff_fft=fft(coeff,524288);
%         pc_fft=Echo_fft.*coeff_fft;
%         pc_freq0=ifft(pc_fft);
%         hug1 = figure('visible','off');
%         subplot(2,1,1);plot(abs(pc_freq0(1:length(Echo)+WaveNum-1)));
%         title('频域脉压结果的幅度,有前暂态点');
%         subplot(2,1,2);
%         plot(abs(pc_time0(1:length(Echo)+WaveNum-1)-...
%                 pc_freq0(1:length(Echo)+WaveNum-1)),'r');
%         title('时域和频域脉压的差别');
%         filename=['时频域脉压对比/扫描' num2str(Sector(argerich)) '度.png'];
%         saveas(hug1,filename)
%         close(gcf)
%         
%         pc_freq1=pc_freq0(WaveNum:length(Echo)+WaveNum-1);
%         %去掉暂态点 WaveNum-1个,后填充点若干(8192-WaveNum+1-length(Echo))
%         % ================按照脉冲号、距离门号重排数据=================================%
%         pc = zeros(PulseNum,SampleNum);
%         for i=1:PulseNum
%             pc(i,1:SampleNum)=pc_freq1((i-1)*SampleNum+1:i*SampleNum);
%             %每个PRT为一行，每行480个采样点的数据
%         end
%         hug2 = figure('visible','off');
%         plot(abs(pc(1,:)));title('频域脉压结果的幅度,没有暂态点');
%         filename=['频域脉压幅度/扫描' num2str(Sector(argerich)) '度.png'];
%         saveas(hug2,filename)
%         close(gcf)
        %% 脉冲压缩 
        coeff_fft=fft(coeff,M);
        for i=1:PulseNum
            s_echo_fft_result=fft(s_echo_mf(1,(i-1)*PRT*Fs+1:i*PRT*Fs),M);
            s_pc_fft=s_echo_fft_result.*coeff_fft;
            pc(i,:)=ifft(s_pc_fft,M);    
        end

        hug1 = figure('visible','off');
        plot(abs(pc(1,:)));%一个周期内三峰值点理论上为40 107 304
    
        filename=['频域脉压幅度/扫描' num2str(Sector(argerich)) '度.png'];
        saveas(hug1,filename)
        close(gcf)
    
        s_pc_result_1=pc;
        s_pc_result_1=reshape((s_pc_result_1)',1,PulseNum*M);   %%%%%%%%%%注意，这里由于reshape函数的算法，需要将矩阵转置才能首尾连在一起
        hug2 = figure('visible','off');
        subplot(2,1,1),plot((0:1/Fs:PulseNum*M/Fs-1/Fs),abs(s_pc_result_1)),
        %N_echo_frame*T_frame-ts
        xlabel('t(单位:s)'),title('脉冲压缩处理后结果（实部）');
        subplot(2,1,2),plot((0:1/Fs:PulseNum*M/Fs-1/Fs),imag(s_pc_result_1)),
        xlabel('t(单位:s)'),title('脉冲压缩处理后结果（虚部）');
        filename=['频域脉压幅度/扫描' num2str(Sector(argerich)) '度.png'];
        saveas(hug2,filename)
        close(gcf)
        %% MTI（动目标显示）,对消静止目标和低速目标---可抑制杂波%
        mti = zeros(PulseNum-1,SampleNum);
        for i=1:PulseNum-1  %滑动对消，少了一个脉冲
            mti(i,:)=pc(i+1,:)-pc(i,:);
        end
        hug3 = figure('visible','off');
        mesh(abs(mti));title('MTI  result');
        filename=['MTI/扫描' num2str(Sector(argerich)) '度.png'];
        saveas(hug3,filename)
        close(gcf)
        %% MTD（动目标检测）,区分不同速度的目标，有测速作用==%
        mtd_PulseNum = PulseNum-1;
        mtd=zeros(PulseNum,SampleNum);
        for i=1:SampleNum
            buff(1:(mtd_PulseNum))= mti(1:(mtd_PulseNum),i);
            buff_fft=fftshift(fft(buff)); %用fftshift将零频搬移到中间 这样可以方便观察速度正负
            mtd(1:mtd_PulseNum,i)=buff_fft(1:mtd_PulseNum)';
        end
        hug4 = figure('visible','off');
        mesh(abs(mtd));title('MTD  result');
        filename=['MTD/扫描' num2str(Sector(argerich)) '度.png'];
        saveas(hug4,filename)
        close(gcf);
        %% 找到多普勒——距离地图中峰值最高的目标 
        abs_mtd = abs(mtd); %对正交分量的复数求模
        max_target = max(max(abs_mtd));
        [row,cell] = find( abs_mtd == max_target);
        target_D = cell/Fs*C/2;
        target_V = ((abs(row - 8))/PulseNum)*PRF*(Lambda/2);
        %% 二维恒虚警 
        ex_PulseNum = mtd_PulseNum+4;
        ex_SampleNum = SampleNum+4;
        ex_mtd = zeros(ex_PulseNum ,ex_SampleNum);
        cfar_mtd = zeros(ex_PulseNum ,ex_SampleNum);
        T = 1.0;     %恒虚警CFAR 阈值因子

        for i = 3:(ex_PulseNum-2)
            for j = 3:(ex_SampleNum-2)         
                  ex_mtd(i,j) = abs_mtd(i-2,j-2); %构建长宽每边都扩张两格的矩阵供恒虚警窗口检测     
            end
        end

        for i = 3:(ex_PulseNum-2)
            for j = 3:(ex_SampleNum-2)    %从有目标的窗口中进行检测
                cfar_sx = ex_mtd(i-2,j-2) + ex_mtd(i-2,j-1) +ex_mtd(i-1,j-2) +ex_mtd(i,j-2)+ex_mtd(i+1,j-2)+ex_mtd(i+2,j-2) + ex_mtd(i+2,j-1);                   
                cfar_sy = ex_mtd(i-2,j+1) + ex_mtd(i-2,j+2) +ex_mtd(i-1,j+2) +ex_mtd(i,j+2)+ex_mtd(i+1,j+2)+ex_mtd(i+2,j+1) + ex_mtd(i+2,j+2);
                cfar_s = T*min(cfar_sx, cfar_sy);      %使用GO—CFAR二维恒虚警，
                if (ex_mtd(i,j)>= cfar_s)
                    cfar_mtd(i,j) = ex_mtd(i,j);
                end
            end
        end
        figure(100);mesh(abs(cfar_mtd));title('cfar  result');  
    %     hug5 = figure('visible','off');
    %     mesh(abs(mtd));title('MTD  result');
    %     filename=['cfar/扫描' num2str(Sector(argerich)) '度.png'];
    %     saveas(hug5,filename)
    %     close(gcf);
    end
end
