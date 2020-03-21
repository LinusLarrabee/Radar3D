function [Chirp,Reality] = SigGeneration(Real)
global BandWidth TimeWidth Fs Fc WaveNum
%Real = 1为实信号，0为复信号，默认为实信号。

Chirp = zeros(1,WaveNum);
for i=-fix(WaveNum/2):fix(WaveNum/2)-1
    Ft = Fc*i/Fs+(1/2)*(BandWidth/TimeWidth)*(i/Fs)^2; %线性调频波指数幂
    if(nargin<1 || Real == 1)
        Chirp(i+fix(WaveNum/2)+1) =cos(2*pi*Ft);
        Reality = 1;
    else 
        Chirp(i+fix(WaveNum/2)+1)=exp(1i*2*pi*Ft); 
        %exp(j*fi)*，产生复数矩阵Chirp
        Reality = 0;
    end
end

hugo01 = figure(1);
subplot(2,1,1),plot(real(Chirp));
subplot(2,1,2);
fftR = fft(Chirp);
plot((0:Fs/WaveNum:Fs/2-Fs/WaveNum),abs(fftR(1:WaveNum/2)));
filename=['测试信号（figure1）.fig'];
saveas(hugo01,filename)